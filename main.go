// PhyloWGS-Go: Pure Go implementation of PhyloWGS MCMC cancer phylogenetics sampler
package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"phylowgs-go/cuda"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

// Global GPU engine (shared across chains)
var (
	gpuEngine *cuda.GPULikelihoodEngine
	useGPU    bool
	gpuMu     sync.Mutex
)

// ============================================================================
// Data Structures
// ============================================================================

// SSM represents a Simple Somatic Mutation
type SSM struct {
	ID              string
	Name            string
	A               []int     // variant read counts per timepoint
	D               []int     // total depth per timepoint
	MuR             float64   // reference allele probability in normal
	MuV             float64   // variant allele probability
	LogBinNormConst []float64 // precomputed log binomial coefficients
	Node            *Node     // assigned node in tree
	CNVs            []*CNVRef // CNVs affecting this SSM

	// Precomputed MH states for CNV-affected SSMs.
	// Mirrors Python's write_data_state / mh.cpp precomputation (params.py:114-171,
	// mh.hpp:80-163): per (SSM, node) the pi-independent (nr, nv) contribution
	// factors for up to 4 timing cases. At MH time the inner loop becomes a
	// plain O(K) dot product with node.Pi[tp] / node.Pi1[tp] instead of a full
	// tree walk via computeNGenomes on every iteration.
	MHStates      []SSMNodeState
	UseFourStates bool // true iff ssm co-located with its CNV
	MHStateValid  bool // cleared whenever tree topology / assignments change
}

// SSMNodeState holds the pi-independent (nr, nv) contribution of a single tree
// node to an SSM's CNV-aware likelihood, for each of the (up to) four timing
// states. The MH inner loop just needs nd.Pi[tp] * Nr_k and nd.Pi[tp] * Nv_k.
type SSMNodeState struct {
	Node               *Node
	Nr1, Nv1           float64
	Nr2, Nv2           float64
	Nr3, Nv3, Nr4, Nv4 float64
}

// SSMDataInterface implementation for cuda.SSMDataInterface
func (s *SSM) GetA() []int                   { return s.A }
func (s *SSM) GetD() []int                   { return s.D }
func (s *SSM) GetMuR() float64               { return s.MuR }
func (s *SSM) GetMuV() float64               { return s.MuV }
func (s *SSM) GetLogBinNormConst() []float64 { return s.LogBinNormConst }

// CNV data interface implementation
func (s *SSM) HasCNV() bool {
	return len(s.CNVs) > 0
}
func (s *SSM) GetMajorCN() int {
	if len(s.CNVs) == 0 {
		return 1 // diploid
	}
	return s.CNVs[0].PaternalCN
}
func (s *SSM) GetMinorCN() int {
	if len(s.CNVs) == 0 {
		return 1 // diploid
	}
	return s.CNVs[0].MaternalCN
}

// CNVRef references a CNV affecting an SSM
type CNVRef struct {
	CNV        *CNV
	MaternalCN int
	PaternalCN int
}

// PhysicalCNV holds genomic annotation for one physical CNV segment.
// Parsed from the physical_cnvs column of cnv_data.txt.
// Fields map directly to the Python cnv_logical_physical_mapping output.
type PhysicalCNV struct {
	Chrom    string `json:"chrom"`
	Start    int    `json:"start"`
	End      int    `json:"end"`
	MajorCN  int    `json:"major_cn"`
	MinorCN  int    `json:"minor_cn"`
	CellPrev string `json:"cell_prev"` // raw "0.0|0.718" string, one value per timepoint
}

// CNVSSMLink holds the SSM ID and copy numbers for one SSM affected by a CNV.
// Used when writing mutlist.json.
type CNVSSMLink struct {
	SSMID      string
	MaternalCN int
	PaternalCN int
}

// CNV represents a Copy Number Variation
type CNV struct {
	ID              string
	A               []int
	D               []int
	LogBinNormConst []float64
	Node            *Node
	AffectedSSMs    []string      // SSM IDs affected by this CNV
	SSMLinks        []CNVSSMLink  // SSM IDs with their CN values (for mutlist.json)
	PhysicalCNVs    []PhysicalCNV // physical segment annotations parsed from cnv_data.txt
}

// Node represents a clone in the phylogenetic tree
type Node struct {
	ID          int
	Parent      *Node
	Children    []*Node
	Data        []int     // indices of SSMs assigned to this node
	Params      []float64 // cellular prevalence (phi) per timepoint
	Pi          []float64 // node-specific prevalence
	Params1     []float64 // proposed params (for MH)
	Pi1         []float64 // proposed pi (for MH)
	Height      int
	Path        []*Node        // ancestors from root to this node
	AncestorSet map[*Node]bool // pre-computed ancestor lookup (for fast isAncestorOf)
}

// TSSB represents the Tree-Structured Stick-Breaking Process
type TSSB struct {
	Root          *TSSBNode
	DPAlpha       float64
	DPGamma       float64
	AlphaDecay    float64
	NumData       int
	Data          []*SSM
	CNVData       []*CNV
	Assignments   []*Node // which node each datum is assigned to
	NTPS          int     // number of timepoints
	MaxDepth      int     // maximum tree depth (default 15)
	MinDepth      int     // minimum depth before stopping (default 0)
	NodeIDCounter int     // counter for generating unique node IDs

	// GPU acceleration buffers
	phiBuf   []float64 // flat phi buffer [nSSM * nTimepoints]
	gpuReady bool      // true if static data uploaded to GPU

	// Cached mixture weights (invalidated on tree structure changes)
	cachedWeights []float64
	cachedNodes   []*Node
	weightsDirty  bool
}

// TSSBNode is a node in the TSSB structure
type TSSBNode struct {
	Node     *Node
	Main     float64   // stick length for stopping here
	Sticks   []float64 // sticks for going to children
	Children []*TSSBNode
}

// ChainResult holds results from a single MCMC chain
type ChainResult struct {
	ChainID     int
	Trees       []TreeSample
	BurninLLH   []float64
	SampleLLH   []float64
	FinalTree   *TSSB
	ElapsedTime time.Duration
}

// TreeSample holds a sampled tree state including a full population snapshot.
// The Snapshot field is a pre-encoded JSON object suitable for writing directly
// to the newline-delimited all_trees.ndjson output file.
type TreeSample struct {
	Iteration int
	LLH       float64
	NumNodes  int
	Snapshot  json.RawMessage // full tree snapshot; nil during burnin
}

// ============================================================================
// Math utilities
// ============================================================================

func logFactorial(n int) float64 {
	if n <= 1 {
		return 0
	}
	// Use math.Lgamma for accuracy and O(1) performance, matching scipy.special.gammaln
	// used by the Python original. The loop-sum alternative accumulates float rounding
	// error at large n (typical WGS read depths).
	lgamma, _ := math.Lgamma(float64(n + 1))
	return lgamma
}

func logBinCoeff(n, k int) float64 {
	if k < 0 || k > n {
		return math.Inf(-1)
	}
	return logFactorial(n) - logFactorial(k) - logFactorial(n-k)
}

func logBinomialLikelihood(a, d int, mu float64) float64 {
	if mu <= 0 {
		mu = 1e-15
	}
	if mu >= 1 {
		mu = 1 - 1e-15
	}
	logMu := math.Log(mu)
	// Use log1p for better numerical stability when mu is close to 0 or 1
	log1MinusMu := math.Log1p(-mu)
	return float64(a)*logMu + float64(d-a)*log1MinusMu
}

// logBinomialLikelihoodPrecomputed uses precomputed log values
func logBinomialLikelihoodPrecomputed(a, d int, logMu, log1MinusMu float64) float64 {
	return float64(a)*logMu + float64(d-a)*log1MinusMu
}

func logsumexp(vals []float64) float64 {
	if len(vals) == 0 {
		return math.Inf(-1)
	}
	maxVal := vals[0]
	for _, v := range vals[1:] {
		if v > maxVal {
			maxVal = v
		}
	}
	if math.IsInf(maxVal, -1) {
		return math.Inf(-1)
	}
	sum := 0.0
	for _, v := range vals {
		sum += math.Exp(v - maxVal)
	}
	return maxVal + math.Log(sum)
}

func boundBeta(a, b float64, rng *rand.Rand) float64 {
	// Simple beta distribution sampling using gamma variates
	x := gammaVariate(a, rng)
	y := gammaVariate(b, rng)
	result := x / (x + y)
	// Bound away from 0 and 1
	if result < 1e-6 {
		result = 1e-6
	}
	if result > 1-1e-6 {
		result = 1 - 1e-6
	}
	return result
}

func gammaVariate(alpha float64, rng *rand.Rand) float64 {
	// Marsaglia and Tsang's method for gamma variates
	if alpha < 1 {
		return gammaVariate(1+alpha, rng) * math.Pow(rng.Float64(), 1/alpha)
	}
	d := alpha - 1.0/3.0
	c := 1.0 / math.Sqrt(9.0*d)
	for {
		var x, v float64
		for {
			x = rng.NormFloat64()
			v = 1.0 + c*x
			if v > 0 {
				break
			}
		}
		v = v * v * v
		u := rng.Float64()
		if u < 1-0.0331*(x*x)*(x*x) {
			return d * v
		}
		if math.Log(u) < 0.5*x*x+d*(1-v+math.Log(v)) {
			return d * v
		}
	}
}

func dirichletSample(alpha []float64, rng *rand.Rand) []float64 {
	result := make([]float64, len(alpha))
	sum := 0.0
	for i, a := range alpha {
		result[i] = gammaVariate(a, rng)
		sum += result[i]
	}
	for i := range result {
		result[i] /= sum
		if result[i] < 1e-10 {
			result[i] = 1e-10
		}
	}
	// Renormalize
	sum = 0
	for _, v := range result {
		sum += v
	}
	for i := range result {
		result[i] /= sum
	}
	return result
}

func dirichletLogPDF(x, alpha []float64) float64 {
	if len(x) != len(alpha) {
		return math.Inf(-1)
	}
	logB := 0.0
	for _, a := range alpha {
		logB += lgamma(a)
	}
	sumAlpha := 0.0
	for _, a := range alpha {
		sumAlpha += a
	}
	logB -= lgamma(sumAlpha)

	result := -logB
	for i := range x {
		if x[i] <= 0 {
			return math.Inf(-1)
		}
		result += (alpha[i] - 1) * math.Log(x[i])
	}
	return result
}

func lgamma(x float64) float64 {
	lg, _ := math.Lgamma(x)
	return lg
}

func betaPDFLn(x, a, b float64) float64 {
	if x <= 0 || x >= 1 {
		return math.Inf(-1)
	}
	return (a-1)*math.Log(x) + (b-1)*math.Log(1-x) - lgamma(a) - lgamma(b) + lgamma(a+b)
}

// ============================================================================
// Data Loading
// ============================================================================

func loadSSMData(filename string) ([]*SSM, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var ssms []*SSM
	scanner := bufio.NewScanner(file)
	lineNum := 0
	for scanner.Scan() {
		lineNum++
		if lineNum == 1 {
			continue // skip header
		}
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		ssm := &SSM{
			ID:   fields[0],
			Name: fields[1],
		}

		// Parse read counts
		aStrs := strings.Split(fields[2], ",")
		dStrs := strings.Split(fields[3], ",")
		ssm.A = make([]int, len(aStrs))
		ssm.D = make([]int, len(dStrs))
		for i, s := range aStrs {
			ssm.A[i], _ = strconv.Atoi(s)
		}
		for i, s := range dStrs {
			ssm.D[i], _ = strconv.Atoi(s)
		}

		// Parse mu values
		ssm.MuR = 0.999
		ssm.MuV = 0.5
		if len(fields) > 4 {
			ssm.MuR, _ = strconv.ParseFloat(fields[4], 64)
		}
		if len(fields) > 5 {
			ssm.MuV, _ = strconv.ParseFloat(fields[5], 64)
		}

		// Precompute log binomial coefficients
		ssm.LogBinNormConst = make([]float64, len(ssm.A))
		for i := range ssm.A {
			ssm.LogBinNormConst[i] = logBinCoeff(ssm.D[i], ssm.A[i])
		}

		ssms = append(ssms, ssm)
	}

	return ssms, scanner.Err()
}

// parsePhysicalCNVs parses the physical_cnvs column from cnv_data.txt.
// Format: "chrom=1,start=141500000,end=148899999,major_cn=2,minor_cn=1,cell_prev=0.0|0.718"
// Multiple segments are separated by ";".
func parsePhysicalCNVs(raw string) []PhysicalCNV {
	var result []PhysicalCNV
	segments := strings.Split(raw, ";")
	for _, seg := range segments {
		seg = strings.TrimSpace(seg)
		if seg == "" {
			continue
		}
		var p PhysicalCNV
		for _, kv := range strings.Split(seg, ",") {
			idx := strings.Index(kv, "=")
			if idx < 0 {
				continue
			}
			key := kv[:idx]
			val := kv[idx+1:]
			switch key {
			case "chrom":
				p.Chrom = val
			case "start":
				p.Start, _ = strconv.Atoi(val)
			case "end":
				p.End, _ = strconv.Atoi(val)
			case "major_cn":
				p.MajorCN, _ = strconv.Atoi(val)
			case "minor_cn":
				p.MinorCN, _ = strconv.Atoi(val)
			case "cell_prev":
				p.CellPrev = val
			}
		}
		result = append(result, p)
	}
	return result
}

func loadCNVData(filename string, ssms []*SSM) ([]*CNV, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	// Build SSM lookup map
	ssmMap := make(map[string]*SSM)
	for _, ssm := range ssms {
		ssmMap[ssm.ID] = ssm
	}

	var cnvs []*CNV
	scanner := bufio.NewScanner(file)
	lineNum := 0
	for scanner.Scan() {
		lineNum++
		if lineNum == 1 {
			continue // skip header
		}
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			continue
		}

		cnv := &CNV{
			ID: fields[0],
		}

		// Parse read counts
		aStrs := strings.Split(fields[1], ",")
		dStrs := strings.Split(fields[2], ",")
		cnv.A = make([]int, len(aStrs))
		cnv.D = make([]int, len(dStrs))
		for i, s := range aStrs {
			cnv.A[i], _ = strconv.Atoi(s)
		}
		for i, s := range dStrs {
			cnv.D[i], _ = strconv.Atoi(s)
		}

		// Precompute log binomial coefficients
		cnv.LogBinNormConst = make([]float64, len(cnv.A))
		for i := range cnv.A {
			cnv.LogBinNormConst[i] = logBinCoeff(cnv.D[i], cnv.A[i])
		}

		// Parse SSM references
		if len(fields) > 3 && fields[3] != "" {
			ssmRefs := strings.Split(fields[3], ";")
			for _, ref := range ssmRefs {
				parts := strings.Split(ref, ",")
				if len(parts) >= 3 {
					ssmID := parts[0]
					maternalCN, _ := strconv.Atoi(parts[1])
					paternalCN, _ := strconv.Atoi(parts[2])
					cnv.AffectedSSMs = append(cnv.AffectedSSMs, ssmID)
					cnv.SSMLinks = append(cnv.SSMLinks, CNVSSMLink{
						SSMID:      ssmID,
						MaternalCN: maternalCN,
						PaternalCN: paternalCN,
					})
					if ssm, ok := ssmMap[ssmID]; ok {
						ssm.CNVs = append(ssm.CNVs, &CNVRef{
							CNV:        cnv,
							MaternalCN: maternalCN,
							PaternalCN: paternalCN,
						})
					}
				}
			}
		}

		// Parse physical_cnvs column (5th column, optional)
		if len(fields) > 4 && fields[4] != "" {
			cnv.PhysicalCNVs = parsePhysicalCNVs(fields[4])
		}

		cnvs = append(cnvs, cnv)
	}

	return cnvs, scanner.Err()
}

// ============================================================================
// Tree Operations
// ============================================================================

func newNode(parent *Node, ntps int) *Node {
	n := &Node{
		Children: make([]*Node, 0),
		Data:     make([]int, 0),
		Params:   make([]float64, ntps),
		Pi:       make([]float64, ntps),
		Params1:  make([]float64, ntps),
		Pi1:      make([]float64, ntps),
		Parent:   parent,
	}
	if parent != nil {
		parent.Children = append(parent.Children, n)
	}
	return n
}

func (n *Node) hasData() bool {
	if len(n.Data) > 0 {
		return true
	}
	for _, child := range n.Children {
		if child.hasData() {
			return true
		}
	}
	return false
}

func (n *Node) getAncestors() []*Node {
	if n.Parent == nil {
		return []*Node{n}
	}
	return append(n.Parent.getAncestors(), n)
}

func newTSSB(ssms []*SSM, cnvs []*CNV, dpAlpha, dpGamma, alphaDecay float64, rng *rand.Rand) *TSSB {
	ntps := len(ssms[0].A)

	// Create root node
	rootNode := newNode(nil, ntps)
	rootNode.ID = 0
	for i := range rootNode.Params {
		rootNode.Params[i] = 1.0
		rootNode.Pi[i] = 1.0
	}

	root := &TSSBNode{
		Node:     rootNode,
		Main:     1e-30, // Root always passes through
		Sticks:   make([]float64, 0),
		Children: make([]*TSSBNode, 0),
	}

	tssb := &TSSB{
		Root:          root,
		DPAlpha:       dpAlpha,
		DPGamma:       dpGamma,
		AlphaDecay:    alphaDecay,
		NumData:       len(ssms),
		Data:          ssms,
		CNVData:       cnvs,
		Assignments:   make([]*Node, len(ssms)),
		NTPS:          ntps,
		MaxDepth:      15,
		MinDepth:      0,
		NodeIDCounter: 1, // root is 0
	}
	rootNode.Height = 0

	// Create initial child node and assign all data to it
	childNode := newNode(rootNode, ntps)
	tssb.NodeIDCounter++
	childNode.ID = tssb.NodeIDCounter
	childNode.Height = 1
	// Initialize child params via pi conservation (Python evolve.py hack section):
	// root node spawns the initial child which takes a random fraction of root's pi.
	// Replicates alleles.__init__ with parent=root: pi = rand(1) * parent.pi.
	// rand(1) is a length-1 numpy array, so broadcasting applies ONE scalar
	// uniformly across all time-point dims — the draw must be hoisted.
	initFrac := rng.Float64()
	for i := range childNode.Params {
		childNode.Pi[i] = initFrac * rootNode.Pi[i]
		rootNode.Pi[i] -= childNode.Pi[i]
		childNode.Params[i] = childNode.Pi[i]
	}

	childTSSB := &TSSBNode{
		Node:     childNode,
		Main:     boundBeta(1.0, alphaDecay*dpAlpha, rng),
		Sticks:   make([]float64, 0),
		Children: make([]*TSSBNode, 0),
	}
	root.Children = append(root.Children, childTSSB)
	root.Sticks = append(root.Sticks, 0.999999)

	// Assign all SSMs to child node
	for i := range ssms {
		childNode.Data = append(childNode.Data, i)
		tssb.Assignments[i] = childNode
		ssms[i].Node = childNode
	}

	// Assign all CNVs to child node (CNVs are also tree data in PhyloWGS)
	for _, cnv := range cnvs {
		cnv.Node = childNode
	}

	return tssb
}

func (t *TSSB) getNodes() []*Node {
	var nodes []*Node
	var descend func(*TSSBNode)
	descend = func(root *TSSBNode) {
		nodes = append(nodes, root.Node)
		for _, child := range root.Children {
			descend(child)
		}
	}
	descend(t.Root)
	return nodes
}

func (t *TSSB) getMixture() ([]float64, []*Node) {
	var weights []float64
	var nodes []*Node

	var descend func(*TSSBNode, float64)
	descend = func(root *TSSBNode, mass float64) {
		weights = append(weights, mass*root.Main)
		nodes = append(nodes, root.Node)

		edges := sticksToEdges(root.Sticks)
		for i, child := range root.Children {
			childMass := mass * (1 - root.Main)
			if i < len(edges) {
				if i == 0 {
					childMass *= edges[i]
				} else {
					childMass *= edges[i] - edges[i-1]
				}
			}
			descend(child, childMass)
		}
	}
	descend(t.Root, 1.0)
	return weights, nodes
}

// getMixtureCached returns cached mixture weights, recomputing if dirty
func (t *TSSB) getMixtureCached() ([]float64, []*Node) {
	if t.weightsDirty || t.cachedWeights == nil {
		t.cachedWeights, t.cachedNodes = t.getMixture()
		t.weightsDirty = false
	}
	return t.cachedWeights, t.cachedNodes
}

// invalidateWeightsCache marks the mixture weight cache as dirty
func (t *TSSB) invalidateWeightsCache() {
	t.weightsDirty = true
}

// rebuildAncestorSets pre-computes ancestor sets for all nodes
// This makes isAncestorOf() O(1) instead of O(depth)
// IMPORTANT: includes the node itself, matching getAncestors() behavior.
// isAncestorOfCached(nd, nd) must return true for computeNGenomes to correctly
// identify when an SSM and its CNV are co-located on the same tree node.
func (t *TSSB) rebuildAncestorSets() {
	nodes := t.getNodes()
	for _, n := range nodes {
		n.AncestorSet = make(map[*Node]bool)
		n.AncestorSet[n] = true // include self: a node is its own ancestor
		cur := n.Parent
		for cur != nil {
			n.AncestorSet[cur] = true
			cur = cur.Parent
		}
	}
}

// isAncestorOfCached uses pre-computed ancestor sets for O(1) lookup
func isAncestorOfCached(candidate, target *Node) bool {
	if target.AncestorSet == nil {
		// Fall back to slow path if not cached
		return isAncestorOf(candidate, target)
	}
	return target.AncestorSet[candidate]
}

func sticksToEdges(sticks []float64) []float64 {
	if len(sticks) == 0 {
		return nil
	}
	edges := make([]float64, len(sticks))
	prod := 1.0
	for i, s := range sticks {
		edges[i] = 1 - prod*(1-s)
		prod *= (1 - s)
	}
	return edges
}

func (t *TSSB) setNodeHeights() {
	t.Root.Node.Height = 0
	var descend func(*TSSBNode, int)
	descend = func(root *TSSBNode, height int) {
		root.Node.Height = height
		for _, child := range root.Children {
			descend(child, height+1)
		}
	}
	for _, child := range t.Root.Children {
		descend(child, 1)
	}
}

func (t *TSSB) setNodePaths() {
	nodes := t.getNodes()
	for _, node := range nodes {
		node.Path = node.getAncestors()
	}
}

func (t *TSSB) mapDatumToNode() {
	nodes := t.getNodes()
	for _, node := range nodes {
		for _, idx := range node.Data {
			if idx < len(t.Data) {
				t.Data[idx].Node = node
			}
		}
	}
}

// ============================================================================
// Likelihood Computation
// ============================================================================

func (ssm *SSM) logLikelihood(phi []float64) float64 {
	if len(ssm.CNVs) == 0 {
		// Simple case: no CNV
		return ssm.logLikelihoodNoCNV(phi)
	}
	// CNV case - more complex
	return ssm.logLikelihoodWithCNV(phi)
}

func (ssm *SSM) logLikelihoodNoCNV(phi []float64) float64 {
	llh := 0.0
	for tp := range ssm.A {
		p := phi[tp]
		if p < 0 {
			p = 0
		}
		if p > 1 {
			p = 1
		}
		mu := (1-p)*ssm.MuR + p*ssm.MuV
		if mu < 1e-15 {
			mu = 1e-15
		}
		if mu > 1-1e-15 {
			mu = 1 - 1e-15
		}
		llh += logBinomialLikelihood(ssm.A[tp], ssm.D[tp], mu) + ssm.LogBinNormConst[tp]
	}
	return llh
}

// findMostRecentCNV finds the most recent (deepest) CNV in a node's ancestry
// that affects the given SSM. Mirrors Python's find_most_recent_cnv().
// Returns nil if no CNV affects this node for this SSM.
func findMostRecentCNV(ssm *SSM, nd *Node) *CNVRef {
	if len(ssm.CNVs) == 0 {
		return nil
	}
	// Walk up ancestors from deepest (most recent) to root
	ancestors := nd.getAncestors()
	for i := len(ancestors) - 1; i >= 0; i-- {
		anc := ancestors[i]
		// Check if any of the SSM's CNVs have their node == anc
		for _, cnvRef := range ssm.CNVs {
			if cnvRef.CNV.Node == anc {
				return cnvRef
			}
		}
	}
	return nil
}

// isAncestorOf returns true if candidate is an ancestor of target
// (i.e., candidate is on the path from root to target)
func isAncestorOf(candidate *Node, target *Node) bool {
	ancestors := target.getAncestors()
	for _, anc := range ancestors {
		if anc == candidate {
			return true
		}
	}
	return false
}

// computeNGenomes computes (nr, nv) pairs for an SSM by traversing ALL tree nodes.
// This mirrors Python's compute_n_genomes() exactly.
// Returns 2 or 4 pairs depending on SSM/CNV relationship.
func computeNGenomes(ssm *SSM, tssb *TSSB, tp int, newState bool) [][2]float64 {
	// Initialize 4 pairs (matching Python: nr1,nv1, nr2,nv2, nr3,nv3, nr4,nv4)
	var nr1, nv1, nr2, nv2, nr3, nv3, nr4, nv4 float64

	// Get SSM's assigned node (equivalent to self.node.path[-1] in Python)
	ssmNode := ssm.Node

	// Traverse ALL nodes in tree
	nodes := tssb.getNodes()
	for _, nd := range nodes {
		// Get pi for this node (mixture weight contribution)
		var pi float64
		if newState {
			pi = nd.Pi1[tp]
		} else {
			pi = nd.Pi[tp]
		}

		// Find most recent CNV affecting this node for this SSM
		mrCNV := findMostRecentCNV(ssm, nd)

		// Check if SSM node is an ancestor of current node
		// (meaning the variant is present in this node's cells)
		// Use cached lookup if available
		ssmInAncestors := isAncestorOfCached(ssmNode, nd)

		// Get CNV copy numbers if present
		var cp, cm int // paternal (major), maternal (minor)
		if mrCNV != nil {
			cp = mrCNV.PaternalCN
			cm = mrCNV.MaternalCN
		}

		if !ssmInAncestors && mrCNV == nil {
			// Case 1: No variant, no CNV - standard diploid (2 ref copies)
			nr1 += pi * 2
			nr2 += pi * 2
			nr3 += pi * 2
			nr4 += pi * 2
		} else if ssmInAncestors && mrCNV == nil {
			// Case 2: Has variant, no CNV - 1 ref, 1 var
			nr1 += pi
			nv1 += pi
			nr2 += pi
			nv2 += pi
			nr3 += pi
			nv3 += pi
			nr4 += pi
			nv4 += pi
		} else if !ssmInAncestors && mrCNV != nil {
			// Case 3: No variant, has CNV - all copies are reference
			totalCN := float64(cp + cm)
			nr1 += pi * totalCN
			nr2 += pi * totalCN
			nr3 += pi * totalCN
			nr4 += pi * totalCN
		} else if ssmInAncestors && mrCNV != nil {
			// Case 4: Has variant AND has CNV - complex timing cases
			totalCN := float64(cp + cm)

			// Cases 3 and 4 in Python: SSM occurred BEFORE CNV
			// (variant not amplified - at most 1 copy of variant)
			nr3 += pi * math.Max(0, totalCN-1)
			nv3 += pi * math.Min(1, totalCN)
			nr4 += pi * math.Max(0, totalCN-1)
			nv4 += pi * math.Min(1, totalCN)

			// Cases 1 and 2: Check timing relationship
			// If SSM node is ancestor of CNV node, SSM occurred first (variant amplified)
			// Python: nr1 = pi * mr_cnv[1] where mr_cnv[1] = tok[1] = minor_cn = MaternalCN = cm
			//         nv1 = pi * mr_cnv[2] where mr_cnv[2] = tok[2] = major_cn = PaternalCN = cp
			//         nr2 = pi * mr_cnv[2] = cp (major/paternal)
			//         nv2 = pi * mr_cnv[1] = cm (minor/maternal)
			cnvNode := mrCNV.CNV.Node
			if cnvNode != nil && isAncestorOfCached(ssmNode, cnvNode) {
				// SSM occurred before CNV - variant gets amplified/deleted by CNV
				// Case 1 (maternal): nr1 = minor_cn (MaternalCN=cm), nv1 = major_cn (PaternalCN=cp)
				nr1 += pi * float64(cm) // was cp - FIXED to match Python mr_cnv[1]=minor
				nv1 += pi * float64(cp) // was cm - FIXED to match Python mr_cnv[2]=major
				// Case 2 (paternal): nr2 = major_cn, nv2 = minor_cn
				nr2 += pi * float64(cp) // was cm - FIXED
				nv2 += pi * float64(cm) // was cp - FIXED
			} else {
				// CNV occurred before or same node as SSM - variant not amplified
				nr1 += pi * math.Max(0, totalCN-1)
				nv1 += pi * math.Min(1, totalCN)
				nr2 += pi * math.Max(0, totalCN-1)
				nv2 += pi * math.Min(1, totalCN)
			}
		}
	}

	// Decide whether to return 2 or 4 pairs
	// Return 4 pairs only if SSM is at the same node as its CNV
	returnFour := false
	if len(ssm.CNVs) == 1 && ssm.CNVs[0].CNV.Node == ssm.Node {
		returnFour = true
	}

	if returnFour {
		return [][2]float64{
			{nr1, nv1},
			{nr2, nv2},
			{nr3, nv3},
			{nr4, nv4},
		}
	}
	return [][2]float64{
		{nr1, nv1},
		{nr2, nv2},
	}
}

// logLikelihoodWithCNVTree computes CNV-aware likelihood using full tree traversal.
// This is the correct implementation matching Python's __log_complete_likelihood__.
func logLikelihoodWithCNVTree(ssm *SSM, tssb *TSSB) float64 {
	if len(ssm.CNVs) == 0 {
		// No CNV - use standard likelihood with assigned node's phi
		return ssm.logLikelihoodNoCNV(ssm.Node.Params)
	}

	llh := 0.0
	for tp := range ssm.A {
		// Compute (nr, nv) pairs for this timepoint
		possNGenomes := computeNGenomes(ssm, tssb, tp, false)

		// Filter out pairs where nv <= 0 (matching Python: nv > 0)
		var validPairs [][2]float64
		for _, ng := range possNGenomes {
			if ng[1] > 0 {
				validPairs = append(validPairs, ng)
			}
		}

		if len(validPairs) == 0 {
			// No valid configurations - return very low likelihood
			llh += math.Log(1e-99)
			continue
		}

		// Compute likelihood for each (nr, nv) pair and average via logsumexp
		prior := math.Log(1.0 / float64(len(validPairs)))
		lls := make([]float64, len(validPairs))
		for i, ng := range validPairs {
			nr, nv := ng[0], ng[1]
			total := nr + nv
			if total < 1e-15 {
				total = 1e-15
			}
			// mu = (nr * muR + nv * (1-muR)) / (nr + nv)
			mu := (nr*ssm.MuR + nv*(1-ssm.MuR)) / total
			if mu < 1e-15 {
				mu = 1e-15
			}
			if mu > 1-1e-15 {
				mu = 1 - 1e-15
			}
			lls[i] = logBinomialLikelihood(ssm.A[tp], ssm.D[tp], mu) + prior + ssm.LogBinNormConst[tp]
		}
		llh += logsumexp(lls)
	}
	return llh
}

// logLikelihoodWithCNV is a wrapper that calls the tree-traversal version
// when a TSSB is available. For backward compatibility.
func (ssm *SSM) logLikelihoodWithCNV(phi []float64) float64 {
	// This function is called without TSSB context, so we fall back to simple model
	// The tree-traversal version (logLikelihoodWithCNVTree) should be called instead
	// when TSSB is available.

	if len(ssm.CNVs) == 0 {
		return ssm.logLikelihoodNoCNV(phi)
	}

	// Simple fallback: use first CNV's copy numbers
	cnvRef := ssm.CNVs[0]
	majorCN := cnvRef.PaternalCN
	minorCN := cnvRef.MaternalCN
	totalCN := majorCN + minorCN

	if totalCN <= 0 {
		return ssm.logLikelihoodNoCNV(phi)
	}

	llh := 0.0
	for tp := range ssm.A {
		p := phi[tp]
		if p < 0 {
			p = 0
		}
		if p > 1 {
			p = 1
		}

		// Two cases: variant on major or minor allele
		var lls []float64

		// Case 1: Variant on major allele
		nr1 := (1-p)*2 + p*float64(minorCN)
		nv1 := p * float64(majorCN)
		if nv1 > 0 {
			total := nr1 + nv1
			mu := (nr1*ssm.MuR + nv1*(1-ssm.MuR)) / total
			mu = math.Max(1e-15, math.Min(1-1e-15, mu))
			lls = append(lls, logBinomialLikelihood(ssm.A[tp], ssm.D[tp], mu)+
				math.Log(0.5)+ssm.LogBinNormConst[tp])
		}

		// Case 2: Variant on minor allele
		nr2 := (1-p)*2 + p*float64(majorCN)
		nv2 := p * float64(minorCN)
		if nv2 > 0 {
			total := nr2 + nv2
			mu := (nr2*ssm.MuR + nv2*(1-ssm.MuR)) / total
			mu = math.Max(1e-15, math.Min(1-1e-15, mu))
			lls = append(lls, logBinomialLikelihood(ssm.A[tp], ssm.D[tp], mu)+
				math.Log(0.5)+ssm.LogBinNormConst[tp])
		}

		if len(lls) == 0 {
			mu := (1-p)*ssm.MuR + p*ssm.MuV
			mu = math.Max(1e-15, math.Min(1-1e-15, mu))
			llh += logBinomialLikelihood(ssm.A[tp], ssm.D[tp], mu) + ssm.LogBinNormConst[tp]
		} else {
			llh += logsumexp(lls)
		}
	}
	return llh
}

func (t *TSSB) completeDataLogLikelihood() float64 {
	weights, nodes := t.getMixture()
	llh := 0.0

	// Try GPU-accelerated path (only for non-CNV SSMs)
	if useGPU && gpuEngine != nil && t.gpuReady {
		gpuMu.Lock()
		gpuLLH, err := t.computeLikelihoodGPU(weights, nodes)
		gpuMu.Unlock()
		if err == nil {
			// GPU computed non-CNV SSMs, now add CNV SSMs via CPU tree traversal
			for _, ssm := range t.Data {
				if len(ssm.CNVs) > 0 {
					// Subtract the simple LLH that GPU computed, add correct tree LLH
					gpuLLH -= ssm.logLikelihoodNoCNV(ssm.Node.Params)
					gpuLLH += logLikelihoodWithCNVTree(ssm, t)
				}
			}
			return gpuLLH
		}
		// Fall through to CPU path on error
	}

	// Build node → weight map for CNV datum contribution
	nodeWeightMap := make(map[*Node]float64)
	for i, node := range nodes {
		nodeWeightMap[node] = weights[i]
	}

	// CPU path - use tree traversal for CNV SSMs
	for i, node := range nodes {
		if len(node.Data) > 0 {
			llh += float64(len(node.Data)) * math.Log(weights[i])
			for _, idx := range node.Data {
				ssm := t.Data[idx]
				if len(ssm.CNVs) > 0 {
					// Use tree-traversal likelihood for CNV SSMs
					llh += logLikelihoodWithCNVTree(ssm, t)
				} else {
					// Standard likelihood for non-CNV SSMs
					llh += ssm.logLikelihoodNoCNV(node.Params)
				}
			}
		}
	}

	// Include CNV datum likelihoods (matching Python: CNVs are data in the tree too)
	// Python treats CNV datums identically to SSMs without their own CNV context.
	// Each CNV datum contributes: log(weight_of_its_node) + binomial_llh(a, d, mu)
	// where mu = (1-phi)*0.999 + phi*0.5 and phi = cnv.Node.Params
	cnvNodeCounts := make(map[*Node]int)
	for _, cnv := range t.CNVData {
		if cnv.Node != nil {
			cnvNodeCounts[cnv.Node]++
		}
	}
	for nd, count := range cnvNodeCounts {
		if w, ok := nodeWeightMap[nd]; ok && w > 0 {
			llh += float64(count) * math.Log(w)
		}
	}
	for _, cnv := range t.CNVData {
		if cnv.Node == nil {
			continue
		}
		phi := cnv.Node.Params
		cnvLLH := 0.0
		for tp := range cnv.A {
			p := phi[tp]
			if p < 0 {
				p = 0
			}
			if p > 1 {
				p = 1
			}
			// CNV datums use same formula: mu = (1-phi)*muR + phi*muV
			mu := (1-p)*0.999 + p*0.5
			if mu < 1e-15 {
				mu = 1e-15
			}
			if mu > 1-1e-15 {
				mu = 1 - 1e-15
			}
			cnvLLH += logBinomialLikelihood(cnv.A[tp], cnv.D[tp], mu) + cnv.LogBinNormConst[tp]
		}
		llh += cnvLLH
	}

	return llh
}

// computeLikelihoodGPU computes log-likelihood using the GPU
// Packs phi values based on current node assignments and calls GPU kernel
func (t *TSSB) computeLikelihoodGPU(weights []float64, nodes []*Node) (float64, error) {
	nSSM := len(t.Data)
	nTP := t.NTPS

	// Ensure phi buffer is allocated
	if len(t.phiBuf) < nSSM*nTP {
		t.phiBuf = make([]float64, nSSM*nTP)
	}

	// Pack phi values for each SSM based on its assigned node
	for i, ssm := range t.Data {
		node := ssm.Node
		if node == nil {
			// Fallback - this shouldn't happen
			for tp := 0; tp < nTP; tp++ {
				t.phiBuf[i*nTP+tp] = 0.5
			}
		} else {
			for tp := 0; tp < nTP; tp++ {
				phi := node.Params[tp]
				if phi < 0 {
					phi = 0
				}
				if phi > 1 {
					phi = 1
				}
				t.phiBuf[i*nTP+tp] = phi
			}
		}
	}

	// Compute likelihoods on GPU
	ssmLLH, err := gpuEngine.ComputeBatchFast(t.phiBuf, nSSM, nTP)
	if err != nil {
		return 0, err
	}

	// Sum up total likelihood with mixture weights
	llh := 0.0

	// Build node weight map
	nodeWeights := make(map[*Node]float64)
	for i, node := range nodes {
		nodeWeights[node] = weights[i]
	}

	// Sum per-SSM likelihoods with log mixture weights
	for i, ssm := range t.Data {
		w := nodeWeights[ssm.Node]
		if w > 0 {
			llh += math.Log(w) + ssmLLH[i]
		}
	}

	return llh, nil
}

// initGPUForTSSB uploads static SSM data to GPU
func (t *TSSB) initGPUForTSSB() error {
	if !useGPU || gpuEngine == nil {
		return nil
	}

	// Convert SSM slice to interface slice
	ssms := make([]cuda.SSMDataInterface, len(t.Data))
	for i, s := range t.Data {
		ssms[i] = s
	}

	gpuMu.Lock()
	defer gpuMu.Unlock()

	err := gpuEngine.UploadStaticData(ssms)
	if err != nil {
		return err
	}

	t.gpuReady = true
	return nil
}

// ============================================================================
// MCMC Moves
// ============================================================================

func (t *TSSB) resampleAssignments(rng *rand.Rand) {
	// Pre-compute tree metadata
	t.setNodeHeights()
	t.setNodePaths()
	t.mapDatumToNode()
	t.rebuildAncestorSets() // Pre-compute ancestor lookups

	eps := 1e-15

	for n := 0; n < t.NumData; n++ {
		oldNode := t.Assignments[n]

		// Get path indices for current assignment
		ancestors := oldNode.getAncestors()
		currentIndices := make([]int, 0)
		current := t.Root
		for _, anc := range ancestors[1:] { // skip root
			for i, child := range current.Children {
				if child.Node == anc {
					currentIndices = append(currentIndices, i)
					current = child
					break
				}
			}
		}

		// Compute old likelihood. For CNV SSMs, use full tree-traversal (like Python).
		var oldLLH float64
		ssm := t.Data[n]
		if len(ssm.CNVs) > 0 {
			oldLLH = logLikelihoodWithCNVTree(ssm, t)
		} else {
			oldLLH = ssm.logLikelihoodNoCNV(oldNode.Params)
		}
		llhSlice := math.Log(rng.Float64()) + oldLLH

		maxU := 1.0
		minU := 0.0

		for iter := 0; iter < 100; iter++ { // Max iterations to prevent infinite loop
			newU := (maxU-minU)*rng.Float64() + minU

			// Use findOrCreateNode for dynamic node creation
			newNode, newPath := t.findOrCreateNode(newU, rng)

			// Skip root node (follow PhyloWGS convention: root should be empty)
			if newNode.Parent == nil && len(t.Root.Children) > 0 {
				newNode = t.Root.Children[0].Node
				newPath = []int{0}
			}

			// Compute new likelihood. For CNV SSMs, temporarily update node assignment
			// to compute correct tree-traversal likelihood (matching Python behavior).
			var newLLH float64
			if len(ssm.CNVs) > 0 {
				// Temporarily move SSM to new node for likelihood evaluation
				oldNode.Data = removeFromSlice(oldNode.Data, n)
				newNode.Data = append(newNode.Data, n)
				ssm.Node = newNode
				newLLH = logLikelihoodWithCNVTree(ssm, t)
				// Restore to old node (proposal not yet accepted)
				newNode.Data = removeFromSlice(newNode.Data, n)
				oldNode.Data = append(oldNode.Data, n)
				ssm.Node = oldNode
			} else {
				newLLH = ssm.logLikelihoodNoCNV(newNode.Params)
			}

			if newLLH > llhSlice {
				// Accept move
				oldNode.Data = removeFromSlice(oldNode.Data, n)
				newNode.Data = append(newNode.Data, n)
				t.Assignments[n] = newNode
				ssm.Node = newNode
				break
			}

			// Shrink slice based on path comparison
			if maxU-minU < eps {
				break
			}

			// Slice bracket shrinking: exact match to Python's tssb.py:
			// path_comp = path_lt(indices, new_path)
			// if path_comp < 0: min_u = new_u
			// elif path_comp >= 0: max_u = new_u
			pathComp := pathLT(currentIndices, newPath)
			if pathComp < 0 {
				minU = newU
			} else {
				maxU = newU
			}
		}
	}

	// Resample CNV datum assignments (matching Python: CNVs are data in the tree too).
	// CNV datums have simple binomial likelihood: mu = (1-phi)*0.999 + phi*0.5
	// They use the same slice sampler as SSMs without CNV context.
	for _, cnv := range t.CNVData {
		if cnv.Node == nil {
			continue
		}
		oldCNVNode := cnv.Node

		// Compute current path for slice shrinking
		cnvAncestors := oldCNVNode.getAncestors()
		cnvCurrentIndices := make([]int, 0)
		cnvCurrent := t.Root
		for _, anc := range cnvAncestors[1:] {
			for i, child := range cnvCurrent.Children {
				if child.Node == anc {
					cnvCurrentIndices = append(cnvCurrentIndices, i)
					cnvCurrent = child
					break
				}
			}
		}

		// Compute old LLH for this CNV datum
		oldCNVLLH := t.cnvDatumLLH(cnv, oldCNVNode.Params)
		llhSlice := math.Log(rng.Float64()) + oldCNVLLH

		maxU := 1.0
		minU := 0.0

		for iter := 0; iter < 100; iter++ {
			newU := (maxU-minU)*rng.Float64() + minU
			newNode, newPath := t.findOrCreateNode(newU, rng)
			if newNode.Parent == nil && len(t.Root.Children) > 0 {
				newNode = t.Root.Children[0].Node
				newPath = []int{0}
			}

			newLLH := t.cnvDatumLLH(cnv, newNode.Params)
			if newLLH > llhSlice {
				cnv.Node = newNode
				break
			}
			if maxU-minU < eps {
				break
			}
			// Path-based shrinking (same as SSMs)
			pathComp := pathLT(cnvCurrentIndices, newPath)
			if pathComp < 0 {
				minU = newU
			} else {
				maxU = newU
			}
		}
	}

	// Invalidate caches since tree structure may have changed
	t.invalidateWeightsCache()
}

// cnvDatumLLH computes the log-likelihood for a CNV datum at a given phi.
// Matches Python: CNV datums use mu = (1-phi)*0.999 + phi*0.5 (no tree traversal).
func (t *TSSB) cnvDatumLLH(cnv *CNV, phi []float64) float64 {
	llh := 0.0
	for tp := range cnv.A {
		p := phi[tp]
		if p < 0 {
			p = 0
		}
		if p > 1 {
			p = 1
		}
		mu := (1-p)*0.999 + p*0.5
		if mu < 1e-15 {
			mu = 1e-15
		}
		if mu > 1-1e-15 {
			mu = 1 - 1e-15
		}
		llh += logBinomialLikelihood(cnv.A[tp], cnv.D[tp], mu) + cnv.LogBinNormConst[tp]
	}
	return llh
}

// removeFromSlice removes the first occurrence of val from a slice of ints.
func removeFromSlice(s []int, val int) []int {
	for i, v := range s {
		if v == val {
			return append(s[:i], s[i+1:]...)
		}
	}
	return s
}

// pathLT implements Python's path_lt(path1, path2) from tssb.py.
// Encodes paths as zero-padded 3-digit strings and compares lexicographically.
// Returns: 1 if path2 > path1, -1 if path2 < path1, 0 if equal.
// (Mirrors Python: (s2 > s1) - (s2 < s1))
func pathLT(path1, path2 []int) int {
	if len(path1) == 0 && len(path2) == 0 {
		return 0
	}
	if len(path1) == 0 {
		return 1 // Python: elif len(path1)==0: return 1
	}
	if len(path2) == 0 {
		return -1 // Python: elif len(path2)==0: return -1
	}
	// Build strings: each index formatted as %03d
	s1 := ""
	for _, v := range path1 {
		s1 += fmt.Sprintf("%03d", v)
	}
	s2 := ""
	for _, v := range path2 {
		s2 += fmt.Sprintf("%03d", v)
	}
	if s2 > s1 {
		return 1
	}
	if s2 < s1 {
		return -1
	}
	return 0
}

func (t *TSSB) metropolis(iters int, std float64, rng *rand.Rand) float64 {
	// Get all nodes
	_, nodes := t.getMixture()
	if len(nodes) == 0 {
		return 0
	}

	// Build node ID map
	nodeIDMap := make(map[int]int)
	for i, node := range nodes {
		node.ID = i
		nodeIDMap[node.ID] = i
	}

	accepted := 0

	// Cache the current likelihood - it only changes when we accept a move
	cachedOldLLH := t.paramPost(nodes, false)

	for iter := 0; iter < iters; iter++ {
		// Sample new pi values for each timepoint
		for tp := 0; tp < t.NTPS; tp++ {
			// Get current pi values
			pi := make([]float64, len(nodes))
			for i, node := range nodes {
				pi[i] = node.Pi[tp]
				if pi[i] < 1e-10 {
					pi[i] = 1e-10
				}
			}

			// Sample from Dirichlet centered on current pi
			alpha := make([]float64, len(pi))
			for i := range pi {
				alpha[i] = std*pi[i] + 1
			}
			piNew := dirichletSample(alpha, rng)

			// Store proposed values
			for i, node := range nodes {
				node.Pi1[tp] = piNew[i]
			}

			// Update proposed params (cumulative prevalence)
			for i := len(nodes) - 1; i >= 0; i-- {
				node := nodes[i]
				node.Params1[tp] = node.Pi1[tp]
				for _, child := range node.Children {
					childIdx := nodeIDMap[child.ID]
					if childIdx < len(nodes) {
						node.Params1[tp] += nodes[childIdx].Params1[tp]
					}
				}
			}
		}

		// Compute acceptance ratio - use cached old LLH
		oldLLH := cachedOldLLH
		newLLH := t.paramPost(nodes, true)
		logA := newLLH - oldLLH

		// Add Dirichlet correction terms
		for tp := 0; tp < t.NTPS; tp++ {
			piOld := make([]float64, len(nodes))
			piNew := make([]float64, len(nodes))
			for i, node := range nodes {
				piOld[i] = node.Pi[tp]
				piNew[i] = node.Pi1[tp]
				if piOld[i] < 1e-10 {
					piOld[i] = 1e-10
				}
				if piNew[i] < 1e-10 {
					piNew[i] = 1e-10
				}
			}

			// Forward proposal: q(pi_new | pi_old)
			alphaFwd := make([]float64, len(piNew))
			for i := range piNew {
				alphaFwd[i] = std * piNew[i]
			}
			logA += dirichletLogPDF(piOld, alphaFwd)

			// Backward proposal: q(pi_old | pi_new)
			alphaBwd := make([]float64, len(piOld))
			for i := range piOld {
				alphaBwd[i] = std * piOld[i]
			}
			logA -= dirichletLogPDF(piNew, alphaBwd)
		}

		// Accept/reject
		if math.Log(rng.Float64()) < logA {
			accepted++
			for _, node := range nodes {
				copy(node.Params, node.Params1)
				copy(node.Pi, node.Pi1)
			}
			// Update cached LLH to new value (since we accepted)
			cachedOldLLH = newLLH
		}
	}

	return float64(accepted) / float64(iters)
}

func (t *TSSB) paramPost(nodes []*Node, useNew bool) float64 {
	// Simple serial computation - parallelization overhead too high for this hot path
	// (called 10,000+ times per MCMC iteration)
	llh := 0.0
	for _, node := range nodes {
		var params []float64
		if useNew {
			params = node.Params1
		} else {
			params = node.Params
		}
		for _, idx := range node.Data {
			ssm := t.Data[idx]
			if len(ssm.CNVs) > 0 {
				llh += logLikelihoodWithCNVTreeMH(ssm, t, useNew)
			} else {
				llh += ssm.logLikelihoodNoCNV(params)
			}
		}
	}
	// Add CNV datum likelihoods (matching Python: CNVs are data nodes too)
	for _, cnv := range t.CNVData {
		if cnv.Node == nil {
			continue
		}
		var phi []float64
		if useNew {
			phi = cnv.Node.Params1
		} else {
			phi = cnv.Node.Params
		}
		for tp := range cnv.A {
			p := phi[tp]
			if p < 0 {
				p = 0
			}
			if p > 1 {
				p = 1
			}
			mu := (1-p)*0.999 + p*0.5
			if mu < 1e-15 {
				mu = 1e-15
			}
			if mu > 1-1e-15 {
				mu = 1 - 1e-15
			}
			llh += logBinomialLikelihood(cnv.A[tp], cnv.D[tp], mu) + cnv.LogBinNormConst[tp]
		}
	}
	return llh
}

// precomputeMHStates fills in per-SSM pi-independent (nr, nv) factors for every
// CNV-affected SSM. Must be called after mapDatumToNode() and before
// metropolis(), and whenever the tree topology, SSM-to-node assignments, or
// CNV-to-node assignments change. Matches Python's params.write_data_state
// (params.py:114-171) semantically.
//
// Parallelized across SSMs with goroutines. The tree is read-only during this
// pass so no locking is required.
func (t *TSSB) precomputeMHStates() {
	// Ensure ancestor sets are current so isAncestorOfCached is valid.
	t.rebuildAncestorSets()
	nodes := t.getNodes()

	// Invalidate any non-CNV SSM precomp (they don't use it anyway, but keep flag tidy).
	for _, ssm := range t.Data {
		if len(ssm.CNVs) == 0 {
			ssm.MHStates = nil
			ssm.MHStateValid = false
			ssm.UseFourStates = false
			continue
		}
		ssm.MHStateValid = false
	}

	// Collect CNV SSMs
	var work []*SSM
	for _, ssm := range t.Data {
		if len(ssm.CNVs) > 0 {
			work = append(work, ssm)
		}
	}
	if len(work) == 0 {
		return
	}

	// Parallelize across SSMs. Each SSM is independent.
	workers := runtime.NumCPU()
	if workers > len(work) {
		workers = len(work)
	}
	if workers < 1 {
		workers = 1
	}
	var wg sync.WaitGroup
	ch := make(chan *SSM, len(work))
	for _, ssm := range work {
		ch <- ssm
	}
	close(ch)
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for ssm := range ch {
				computeSSMStates(ssm, nodes)
			}
		}()
	}
	wg.Wait()
}

// computeSSMStates fills ssm.MHStates with per-node (nr, nv) coefficients for
// each of the 4 timing cases. Matches computeNGenomes case-by-case exactly, but
// stores the pi-independent coefficient rather than aggregating across nodes.
func computeSSMStates(ssm *SSM, nodes []*Node) {
	ssmNode := ssm.Node
	if cap(ssm.MHStates) < len(nodes) {
		ssm.MHStates = make([]SSMNodeState, len(nodes))
	} else {
		ssm.MHStates = ssm.MHStates[:len(nodes)]
		for i := range ssm.MHStates {
			ssm.MHStates[i] = SSMNodeState{}
		}
	}

	for i, nd := range nodes {
		mrCNV := findMostRecentCNV(ssm, nd)
		ssmInAncestors := isAncestorOfCached(ssmNode, nd)
		var cp, cm int
		if mrCNV != nil {
			cp = mrCNV.PaternalCN
			cm = mrCNV.MaternalCN
		}

		s := SSMNodeState{Node: nd}

		switch {
		case !ssmInAncestors && mrCNV == nil:
			// Case 1: no variant, no CNV -> diploid reference
			s.Nr1 = 2
			s.Nr2 = 2
			s.Nr3 = 2
			s.Nr4 = 2
		case ssmInAncestors && mrCNV == nil:
			// Case 2: variant, no CNV -> 1 ref, 1 var
			s.Nr1, s.Nv1 = 1, 1
			s.Nr2, s.Nv2 = 1, 1
			s.Nr3, s.Nv3 = 1, 1
			s.Nr4, s.Nv4 = 1, 1
		case !ssmInAncestors && mrCNV != nil:
			// Case 3: no variant, has CNV -> all reference, total copies
			totalCN := float64(cp + cm)
			s.Nr1 = totalCN
			s.Nr2 = totalCN
			s.Nr3 = totalCN
			s.Nr4 = totalCN
		case ssmInAncestors && mrCNV != nil:
			// Case 4: variant + CNV -> timing-dependent
			totalCN := float64(cp + cm)
			// States 3,4: SSM occurred before CNV is not amplified (unamplified)
			s.Nr3 = math.Max(0, totalCN-1)
			s.Nv3 = math.Min(1, totalCN)
			s.Nr4 = math.Max(0, totalCN-1)
			s.Nv4 = math.Min(1, totalCN)

			// States 1,2: depends on SSM-CNV timing
			cnvNode := mrCNV.CNV.Node
			if cnvNode != nil && isAncestorOfCached(ssmNode, cnvNode) {
				// SSM before CNV: variant amplified
				s.Nr1 = float64(cm)
				s.Nv1 = float64(cp)
				s.Nr2 = float64(cp)
				s.Nv2 = float64(cm)
			} else {
				// CNV before or at SSM: variant not amplified
				s.Nr1 = math.Max(0, totalCN-1)
				s.Nv1 = math.Min(1, totalCN)
				s.Nr2 = math.Max(0, totalCN-1)
				s.Nv2 = math.Min(1, totalCN)
			}
		}

		ssm.MHStates[i] = s
	}

	// Decide whether to return 2 or 4 pairs at runtime. Matches computeNGenomes.
	ssm.UseFourStates = len(ssm.CNVs) == 1 && ssm.CNVs[0].CNV.Node == ssm.Node
	ssm.MHStateValid = true
}

// logLikelihoodWithCNVTreeMHPrecomputed is the O(K) fast path used during MH
// when ssm.MHStateValid is true. Mirrors Python's mh.cpp log_complete_ll inner
// loop: a plain dot product of precomputed (nr, nv) factors with node.Pi[tp]
// (or node.Pi1[tp] when newState=true).
func logLikelihoodWithCNVTreeMHPrecomputed(ssm *SSM, newState bool) float64 {
	llh := 0.0
	states := ssm.MHStates
	useFour := ssm.UseFourStates

	for tp := range ssm.A {
		var nr1, nv1, nr2, nv2, nr3, nv3, nr4, nv4 float64
		if newState {
			for k := range states {
				s := &states[k]
				pi := s.Node.Pi1[tp]
				nr1 += pi * s.Nr1
				nv1 += pi * s.Nv1
				nr2 += pi * s.Nr2
				nv2 += pi * s.Nv2
				if useFour {
					nr3 += pi * s.Nr3
					nv3 += pi * s.Nv3
					nr4 += pi * s.Nr4
					nv4 += pi * s.Nv4
				}
			}
		} else {
			for k := range states {
				s := &states[k]
				pi := s.Node.Pi[tp]
				nr1 += pi * s.Nr1
				nv1 += pi * s.Nv1
				nr2 += pi * s.Nr2
				nv2 += pi * s.Nv2
				if useFour {
					nr3 += pi * s.Nr3
					nv3 += pi * s.Nv3
					nr4 += pi * s.Nr4
					nv4 += pi * s.Nv4
				}
			}
		}

		// Build pairs and apply the same nv > 0 filter as computeNGenomes.
		var pairs [4][2]float64
		n := 0
		if nv1 > 0 {
			pairs[n] = [2]float64{nr1, nv1}
			n++
		}
		if nv2 > 0 {
			pairs[n] = [2]float64{nr2, nv2}
			n++
		}
		if useFour {
			if nv3 > 0 {
				pairs[n] = [2]float64{nr3, nv3}
				n++
			}
			if nv4 > 0 {
				pairs[n] = [2]float64{nr4, nv4}
				n++
			}
		}

		if n == 0 {
			llh += math.Log(1e-99)
			continue
		}

		prior := math.Log(1.0 / float64(n))
		var lls [4]float64
		for i := 0; i < n; i++ {
			nr, nv := pairs[i][0], pairs[i][1]
			total := nr + nv
			if total < 1e-15 {
				total = 1e-15
			}
			mu := (nr*ssm.MuR + nv*(1-ssm.MuR)) / total
			if mu < 1e-15 {
				mu = 1e-15
			}
			if mu > 1-1e-15 {
				mu = 1 - 1e-15
			}
			lls[i] = logBinomialLikelihood(ssm.A[tp], ssm.D[tp], mu) + prior + ssm.LogBinNormConst[tp]
		}
		llh += logsumexp(lls[:n])
	}
	return llh
}

// logLikelihoodWithCNVTreeMH is like logLikelihoodWithCNVTree but supports MH state
func logLikelihoodWithCNVTreeMH(ssm *SSM, tssb *TSSB, newState bool) float64 {
	if len(ssm.CNVs) == 0 {
		var params []float64
		if newState {
			params = ssm.Node.Params1
		} else {
			params = ssm.Node.Params
		}
		return ssm.logLikelihoodNoCNV(params)
	}

	// Fast path: precomputed per-node (nr, nv) factors (mirrors Python mh.cpp).
	if ssm.MHStateValid {
		return logLikelihoodWithCNVTreeMHPrecomputed(ssm, newState)
	}

	llh := 0.0
	for tp := range ssm.A {
		// Compute (nr, nv) pairs for this timepoint using newState flag
		possNGenomes := computeNGenomes(ssm, tssb, tp, newState)

		// Filter out pairs where nv <= 0
		var validPairs [][2]float64
		for _, ng := range possNGenomes {
			if ng[1] > 0 {
				validPairs = append(validPairs, ng)
			}
		}

		if len(validPairs) == 0 {
			llh += math.Log(1e-99)
			continue
		}

		prior := math.Log(1.0 / float64(len(validPairs)))
		lls := make([]float64, len(validPairs))
		for i, ng := range validPairs {
			nr, nv := ng[0], ng[1]
			total := nr + nv
			if total < 1e-15 {
				total = 1e-15
			}
			mu := (nr*ssm.MuR + nv*(1-ssm.MuR)) / total
			if mu < 1e-15 {
				mu = 1e-15
			}
			if mu > 1-1e-15 {
				mu = 1 - 1e-15
			}
			lls[i] = logBinomialLikelihood(ssm.A[tp], ssm.D[tp], mu) + prior + ssm.LogBinNormConst[tp]
		}
		llh += logsumexp(lls)
	}
	return llh
}

func (t *TSSB) resampleSticks(rng *rand.Rand) {
	var descend func(*TSSBNode, int)
	descend = func(root *TSSBNode, depth int) {
		dataDown := 0
		// Process children in reverse order
		for i := len(root.Children) - 1; i >= 0; i-- {
			child := root.Children[i]
			childData := countData(child)
			descend(child, depth+1)

			postAlpha := 1.0 + float64(childData)
			postBeta := t.DPGamma + float64(dataDown)
			if depth != 0 {
				root.Sticks[i] = boundBeta(postAlpha, postBeta, rng)
			} else {
				root.Sticks[i] = 0.999999
			}
			dataDown += childData
		}

		dataHere := len(root.Node.Data)
		postAlpha := 1.0 + float64(dataHere)
		postBeta := math.Pow(t.AlphaDecay, float64(depth))*t.DPAlpha + float64(dataDown)
		if depth >= 1 { // min_depth = 0 typically
			root.Main = boundBeta(postAlpha, postBeta, rng)
		} else {
			root.Main = 1e-30
		}
	}
	descend(t.Root, 0)
	t.invalidateWeightsCache()
}

func countData(node *TSSBNode) int {
	count := len(node.Node.Data)
	for _, child := range node.Children {
		count += countData(child)
	}
	return count
}

// killNode removes a child from its parent, returning the child's pi to the parent.
// Matches Python's alleles.kill().
func killNode(child *Node, parent *Node) {
	// Return pi to parent
	for i := range parent.Pi {
		parent.Pi[i] += child.Pi[i]
	}
	// Remove from parent's Node.Children
	for i, c := range parent.Children {
		if c == child {
			parent.Children = append(parent.Children[:i], parent.Children[i+1:]...)
			break
		}
	}
}

func (t *TSSB) cullTree() {
	var descend func(*TSSBNode) int
	descend = func(root *TSSBNode) int {
		counts := make([]int, len(root.Children))
		for i, child := range root.Children {
			counts[i] = descend(child)
		}

		// Keep only children with data
		keep := 0
		for i := len(counts) - 1; i >= 0; i-- {
			if counts[i] > 0 {
				keep = i + 1
				break
			}
		}

		// Remove empty trailing children -- return pi to parent (like Python's kill())
		if keep < len(root.Children) {
			for i := keep; i < len(root.Children); i++ {
				killNode(root.Children[i].Node, root.Node)
			}
			root.Children = root.Children[:keep]
			if keep < len(root.Sticks) {
				root.Sticks = root.Sticks[:keep]
			}
		}

		total := len(root.Node.Data)
		for _, c := range counts[:keep] {
			total += c
		}
		return total
	}
	descend(t.Root)
	t.invalidateWeightsCache()
}

// stickOrderSubWeights builds the normalized sub-weight vector for one
// iteration of resample_stick_orders. Matches Python (tssb.py:208-210):
//
//	sub_weights = hstack([all_weights[sub_indices], 1.0 - sum(all_weights)])
//	sub_weights = sub_weights / sum(sub_weights)
//
// The "unallocated region" bucket is 1 - sum(ALL stick weights), i.e. the
// true unallocated stick mass. Go previously (incorrectly) computed the
// bucket as 1 - sum(NOT-yet-ordered weights), conflating already-ordered
// children's mass with genuinely unallocated mass and inflating P(spawn).
//
// Returns a slice of length (len(allWeights) - len(newOrder) + 1), with the
// remaining children first (in original index order) and the leftover bucket
// last. All entries are normalized to sum to 1.
func stickOrderSubWeights(allWeights []float64, newOrder []int) []float64 {
	// Build subIndices: indices NOT yet in newOrder, in original order.
	inOrder := make(map[int]bool, len(newOrder))
	for _, k := range newOrder {
		inOrder[k] = true
	}
	subIndices := make([]int, 0, len(allWeights)-len(newOrder))
	for i := 0; i < len(allWeights); i++ {
		if !inOrder[i] {
			subIndices = append(subIndices, i)
		}
	}

	subWeights := make([]float64, len(subIndices)+1)
	for i, idx := range subIndices {
		subWeights[i] = allWeights[idx]
	}

	// Leftover bucket = 1 - sum(ALL stick weights) (Python's formula).
	totalAll := 0.0
	for _, w := range allWeights {
		totalAll += w
	}
	leftover := 1.0 - totalAll
	if leftover < 0 {
		leftover = 0
	}
	subWeights[len(subIndices)] = leftover

	// Normalize.
	sumW := 0.0
	for _, w := range subWeights {
		sumW += w
	}
	if sumW > 0 {
		for i := range subWeights {
			subWeights[i] /= sumW
		}
	}
	return subWeights
}

// spawnChild creates a new child TSSBNode from a parent TSSBNode
// Equivalent to Python's root['node'].spawn() pattern
func (t *TSSB) spawnChild(parent *TSSBNode, depth int, rng *rand.Rand) *TSSBNode {
	// Create new Node
	childNode := newNode(parent.Node, t.NTPS)
	t.NodeIDCounter++
	childNode.ID = t.NodeIDCounter
	childNode.Height = depth + 1

	// Pi conservation: child takes a random fraction of parent's pi,
	// parent gives up that fraction. Matches Python alleles.__init__:
	//   self.pi = rand(1) * parent.pi
	//   parent.pi = parent.pi - self.pi
	//   self.params = self.pi
	// NOTE: rand(1) returns a length-1 numpy array; broadcasting (1,) * (ntps,)
	// applies ONE scalar uniformly across all time-point dims. The draw MUST
	// be hoisted outside the per-dim loop to preserve this proportionality.
	frac := rng.Float64()
	for i := range childNode.Pi {
		childNode.Pi[i] = frac * parent.Node.Pi[i]
		parent.Node.Pi[i] -= childNode.Pi[i]
		childNode.Params[i] = childNode.Pi[i]
	}

	// Compute main stick break
	var main float64
	if t.MinDepth <= depth+1 {
		main = boundBeta(1.0, math.Pow(t.AlphaDecay, float64(depth+1))*t.DPAlpha, rng)
	} else {
		main = 0.0
	}

	return &TSSBNode{
		Node:     childNode,
		Main:     main,
		Sticks:   make([]float64, 0),
		Children: make([]*TSSBNode, 0),
	}
}

// resampleStickOrders reorders children of each node by stick lengths
// and prunes branches with no assigned data.
// Equivalent to original tssb.resample_stick_orders()
func (t *TSSB) resampleStickOrders(rng *rand.Rand) {
	var descend func(*TSSBNode, int)
	descend = func(root *TSSBNode, depth int) {
		if len(root.Children) == 0 {
			return
		}

		// Find children that have data (represented set)
		represented := make(map[int]bool)
		for i, child := range root.Children {
			if child.Node.hasData() {
				represented[i] = true
			}
		}

		// If no children have data, skip reordering (they'll be pruned anyway)
		if len(represented) == 0 {
			return
		}

		// Compute weights from sticks
		edges := sticksToEdges(root.Sticks)
		allWeights := make([]float64, len(edges))
		if len(edges) > 0 {
			allWeights[0] = edges[0]
			for i := 1; i < len(edges); i++ {
				allWeights[i] = edges[i] - edges[i-1]
			}
		}

		newOrder := make([]int, 0)

		// Sample children in random order, weighted by stick lengths
		for len(represented) > 0 {
			u := rng.Float64()

			// Build sub-weights for children not yet in new_order
			// Safeguard: maximum number of child creations per sampling
			maxCreations := 10
			creations := 0

			for {
				// Rebuild subIndices (indices not yet in newOrder, in order).
				subIndices := make([]int, 0)
				{
					inOrder := make(map[int]bool, len(newOrder))
					for _, k := range newOrder {
						inOrder[k] = true
					}
					for i := 0; i < len(root.Sticks); i++ {
						if !inOrder[i] {
							subIndices = append(subIndices, i)
						}
					}
				}

				// If no indices left to sample, break
				if len(subIndices) == 0 {
					break
				}

				// Compute sub-weights via Python-faithful helper (leftover
				// bucket = 1 - sum(ALL weights), not 1 - sum(remaining)).
				subWeights := stickOrderSubWeights(allWeights, newOrder)

				// If all sub-weights are zero (degenerate), pick first remaining.
				sumW := 0.0
				for _, w := range subWeights {
					sumW += w
				}
				if sumW == 0 {
					if len(subIndices) > 0 {
						newOrder = append(newOrder, subIndices[0])
						delete(represented, subIndices[0])
					}
					break
				}

				// Sample index from cumulative weights
				cumSum := 0.0
				index := len(subIndices) // default to leftover
				for i, w := range subWeights {
					cumSum += w
					if u < cumSum {
						index = i
						break
					}
				}

				if index == len(subIndices) {
					// Fell into uncreated region - create new child
					creations++
					if creations > maxCreations {
						// Too many creations - force pick from existing
						if len(subIndices) > 0 {
							// Pick the first represented one, or first available
							picked := -1
							for _, si := range subIndices {
								if represented[si] {
									picked = si
									break
								}
							}
							if picked == -1 && len(subIndices) > 0 {
								picked = subIndices[0]
							}
							if picked >= 0 {
								newOrder = append(newOrder, picked)
								delete(represented, picked)
							}
						}
						break
					}

					newStick := boundBeta(1, t.DPGamma, rng)
					root.Sticks = append(root.Sticks, newStick)
					newChild := t.spawnChild(root, depth, rng)
					root.Children = append(root.Children, newChild)

					// Recompute weights
					edges = sticksToEdges(root.Sticks)
					allWeights = make([]float64, len(edges))
					if len(edges) > 0 {
						allWeights[0] = edges[0]
						for i := 1; i < len(edges); i++ {
							allWeights[i] = edges[i] - edges[i-1]
						}
					}
				} else {
					index = subIndices[index]
					newOrder = append(newOrder, index)
					delete(represented, index)
					break
				}
			}
		}

		// Build new children list and recurse
		newChildren := make([]*TSSBNode, len(newOrder))
		for i, k := range newOrder {
			newChildren[i] = root.Children[k]
			descend(root.Children[k], depth+1)
		}

		// Kill children not in new order (no data) -- return pi to parent
		inNewOrder := make(map[int]bool)
		for _, k := range newOrder {
			inNewOrder[k] = true
		}
		for k := 0; k < len(root.Children); k++ {
			if !inNewOrder[k] {
				// Kill node: return pi to parent and remove from Node.Children
				killNode(root.Children[k].Node, root.Node)
			}
		}

		root.Children = newChildren
		root.Sticks = make([]float64, len(newChildren))
	}

	descend(t.Root, 0)

	// Immediately resample sticks after reordering
	t.resampleSticks(rng)
}

// findOrCreateNode traverses the stick-breaking tree for position u,
// creating new nodes as needed (lazy instantiation).
// Equivalent to original tssb.find_node()
func (t *TSSB) findOrCreateNode(u float64, rng *rand.Rand) (*Node, []int) {
	var descend func(*TSSBNode, float64, int) (*Node, []int)
	descend = func(root *TSSBNode, u float64, depth int) (*Node, []int) {
		// Check max depth
		if depth >= t.MaxDepth {
			return root.Node, []int{}
		}

		// Stop at this node?
		if u < root.Main {
			return root.Node, []int{}
		}

		// Rescale u to remaining interval
		u = (u - root.Main) / (1.0 - root.Main)

		if depth > 0 {
			// Create children as needed until u falls within existing stick space
			// Safeguard: limit number of children to prevent infinite loops
			maxChildCreations := 20
			creations := 0
			for {
				if len(root.Children) == 0 {
					// No children yet - create first one
					newStick := boundBeta(1, t.DPGamma, rng)
					root.Sticks = append(root.Sticks, newStick)
					newChild := t.spawnChild(root, depth, rng)
					root.Children = append(root.Children, newChild)
					creations++
				}

				// Compute edge = 1 - prod(1 - sticks)
				prod := 1.0
				for _, s := range root.Sticks {
					prod *= (1.0 - s)
				}
				edge := 1.0 - prod

				if u < edge || creations >= maxChildCreations {
					break // u falls within existing sticks OR we've created enough
				}

				// u falls beyond existing sticks - create new child
				newStick := boundBeta(1, t.DPGamma, rng)
				root.Sticks = append(root.Sticks, newStick)
				newChild := t.spawnChild(root, depth, rng)
				root.Children = append(root.Children, newChild)
				creations++
			}

			// Find which child u falls into
			edges := make([]float64, len(root.Sticks))
			prod := 1.0
			for i, s := range root.Sticks {
				prod *= (1.0 - s)
				edges[i] = 1.0 - prod
			}

			// Find index where u <= edges[index]
			index := 0
			for i, e := range edges {
				if u <= e {
					index = i
					break
				}
				index = i
			}

			// Rescale u for the child
			edgeLower := 0.0
			if index > 0 {
				edgeLower = edges[index-1]
			}
			edgeUpper := edges[index]
			if edgeUpper > edgeLower {
				u = (u - edgeLower) / (edgeUpper - edgeLower)
			} else {
				u = 0
			}

			node, path := descend(root.Children[index], u, depth+1)
			return node, append([]int{index}, path...)
		} else {
			// At root (depth 0) - always go to first child.
			// Python's find_node unconditionally does path.insert(0, index) after
			// both the depth>0 and depth==0 branches. So root-level index IS included.
			if len(root.Children) == 0 {
				// Should not happen in normal use, but handle gracefully
				return root.Node, []int{}
			}
			index := 0
			node, path := descend(root.Children[index], u, depth+1)
			return node, append([]int{index}, path...)
		}
	}

	return descend(t.Root, u, 0)
}

// dpAlphaLLH computes the log-likelihood of the DP-alpha hyperparameter over
// the current tree. Matches Python's TSSB.dp_alpha_llh (tssb.py:247-255):
//
//	if self.min_depth <= depth:
//	    llh += betapdfln(root['main'], 1.0, (self.alpha_decay**depth)*dp_alpha)
//
// With the default min_depth=0, the root's main-stick beta term IS included;
// it is NOT a constant that cancels in the slice sampler's LLH comparison
// because it depends on alpha at depth 0.
func (t *TSSB) dpAlphaLLH(alpha float64) float64 {
	var descend func(*TSSBNode, int) float64
	descend = func(root *TSSBNode, depth int) float64 {
		llh := 0.0
		if t.MinDepth <= depth {
			llh = betaPDFLn(root.Main, 1.0, math.Pow(t.AlphaDecay, float64(depth))*alpha)
		}
		for _, child := range root.Children {
			llh += descend(child, depth+1)
		}
		return llh
	}
	return descend(t.Root, 0)
}

func (t *TSSB) resampleHypers(rng *rand.Rand) {
	// Resample dp_alpha using slice sampler
	minDPAlpha := 1.0
	maxDPAlpha := 50.0

	// Slice sample for dp_alpha
	llhSlice := math.Log(rng.Float64()) + t.dpAlphaLLH(t.DPAlpha)
	lower := minDPAlpha
	upper := maxDPAlpha
	for iter := 0; iter < 100; iter++ {
		newAlpha := (upper-lower)*rng.Float64() + lower
		newLLH := t.dpAlphaLLH(newAlpha)
		if newLLH > llhSlice {
			t.DPAlpha = newAlpha
			break
		}
		if newAlpha < t.DPAlpha {
			lower = newAlpha
		} else {
			upper = newAlpha
		}
	}

	// Resample alpha_decay
	minDecay := 0.05
	maxDecay := 0.80
	llhSlice = math.Log(rng.Float64()) + t.dpAlphaLLH(t.DPAlpha)
	lower = minDecay
	upper = maxDecay
	for iter := 0; iter < 100; iter++ {
		newDecay := (upper-lower)*rng.Float64() + lower
		oldDecay := t.AlphaDecay
		t.AlphaDecay = newDecay
		newLLH := t.dpAlphaLLH(t.DPAlpha)
		t.AlphaDecay = oldDecay
		if newLLH > llhSlice {
			t.AlphaDecay = newDecay
			break
		}
		if newDecay < t.AlphaDecay {
			lower = newDecay
		} else {
			upper = newDecay
		}
	}

	// Resample dp_gamma (stick concentration for children)
	// Matches Python's tssb.resample_hypers(dp_gamma=True)
	minDPGamma := 1.0
	maxDPGamma := 10.0
	dpGammaLLH := func(gamma float64) float64 {
		var descend func(*TSSBNode) float64
		descend = func(root *TSSBNode) float64 {
			llh := 0.0
			for i, child := range root.Children {
				if i < len(root.Sticks) {
					llh += betaPDFLn(root.Sticks[i], 1.0, gamma)
				}
				llh += descend(child)
			}
			return llh
		}
		return descend(t.Root)
	}
	llhSlice = math.Log(rng.Float64()) + dpGammaLLH(t.DPGamma)
	lower = minDPGamma
	upper = maxDPGamma
	for iter := 0; iter < 100; iter++ {
		newGamma := (upper-lower)*rng.Float64() + lower
		newLLH := dpGammaLLH(newGamma)
		if newLLH > llhSlice {
			t.DPGamma = newGamma
			break
		}
		if newGamma < t.DPGamma {
			lower = newGamma
		} else {
			upper = newGamma
		}
	}
}

// ============================================================================
// MCMC Chain
// ============================================================================

func runChain(chainID int, ssms []*SSM, cnvs []*CNV, burnin, samples, mhIters int, seed int64) ChainResult {
	start := time.Now()
	rng := rand.New(rand.NewSource(seed))

	// Deep copy CNVs for this chain.
	// CNV.Node is written on every MCMC iteration (resampleAssignments). All
	// chains run concurrently, so each chain must own its own *CNV objects to
	// avoid data races and cross-chain corruption of node assignments.
	chainCNVs := make([]*CNV, len(cnvs))
	origToChainCNV := make(map[*CNV]*CNV, len(cnvs))
	for i, cnv := range cnvs {
		chainCNV := &CNV{
			ID:              cnv.ID,
			A:               append([]int{}, cnv.A...),
			D:               append([]int{}, cnv.D...),
			LogBinNormConst: append([]float64{}, cnv.LogBinNormConst...),
			AffectedSSMs:    cnv.AffectedSSMs, // read-only during MCMC
			SSMLinks:        cnv.SSMLinks,      // read-only during MCMC
			PhysicalCNVs:    cnv.PhysicalCNVs,  // read-only during MCMC
			// Node is assigned by newTSSB and updated by resampleAssignments
		}
		chainCNVs[i] = chainCNV
		origToChainCNV[cnv] = chainCNV
	}

	// Deep copy SSMs for this chain.
	// Each SSM gets its own CNVRef slice pointing to chain-local *CNV copies so
	// that findMostRecentCNV / computeNGenomes see the correct per-chain node
	// assignment for each CNV.
	ssmsCopy := make([]*SSM, len(ssms))
	for i, ssm := range ssms {
		chainRefs := make([]*CNVRef, len(ssm.CNVs))
		for j, ref := range ssm.CNVs {
			chainRefs[j] = &CNVRef{
				CNV:        origToChainCNV[ref.CNV],
				MaternalCN: ref.MaternalCN,
				PaternalCN: ref.PaternalCN,
			}
		}
		ssmsCopy[i] = &SSM{
			ID:              ssm.ID,
			Name:            ssm.Name,
			A:               append([]int{}, ssm.A...),
			D:               append([]int{}, ssm.D...),
			MuR:             ssm.MuR,
			MuV:             ssm.MuV,
			LogBinNormConst: append([]float64{}, ssm.LogBinNormConst...),
			CNVs:            chainRefs,
		}
	}

	// Initialize TSSB with chain-local SSMs and CNVs
	tssb := newTSSB(ssmsCopy, chainCNVs, 25.0, 1.0, 0.25, rng)

	// Initialize GPU for this TSSB (uploads static SSM data)
	if err := tssb.initGPUForTSSB(); err != nil {
		log.Printf("Chain %d: GPU init failed, using CPU: %v", chainID, err)
	}

	mhStd := 100.0
	var burninLLH []float64
	var sampleLLH []float64
	var trees []TreeSample

	totalIters := burnin + samples

	for iter := -burnin; iter < samples; iter++ {
		// MCMC iteration
		tssb.resampleAssignments(rng)
		tssb.cullTree()

		// Pre-compute tree metadata for MH
		tssb.setNodeHeights()
		tssb.setNodePaths()
		tssb.mapDatumToNode()

		// Precompute per-(SSM, node) (nr, nv) factors once per MH loop so the
		// metropolis inner loop is an O(K) dot product instead of a full
		// tree walk via computeNGenomes on every iteration. Mirrors Python's
		// params.write_data_state() + mh.cpp precomputation.
		tssb.precomputeMHStates()

		// MH sampling for params
		mhAcc := tssb.metropolis(mhIters, mhStd, rng)

		// Invalidate precomputed states after MH: subsequent operations
		// (resampleSticks, resampleStickOrders, resampleAssignments in the
		// next iteration) may change tree topology or SSM-to-node
		// assignments, so the cached (nr, nv) factors are no longer valid.
		for _, ssm := range tssb.Data {
			ssm.MHStateValid = false
		}

		// Adapt MH step size
		if mhAcc < 0.08 && mhStd < 10000 {
			mhStd *= 2.0
		}
		if mhAcc > 0.5 && mhAcc < 0.99 {
			mhStd /= 2.0
		}

		// Resample tree structure
		tssb.resampleSticks(rng)
		tssb.resampleStickOrders(rng) // Reorder children and prune empty branches
		tssb.resampleHypers(rng)

		// Compute likelihood
		llh := tssb.completeDataLogLikelihood()

		if iter < 0 {
			burninLLH = append(burninLLH, llh)
		} else {
			sampleLLH = append(sampleLLH, llh)
			_, nodes := tssb.getMixture()
			// Capture a full tree snapshot for the posterior output file.
			// snapshotTree does not call removeEmptyNodes so the live tree is
			// unaffected; cullTree + resampleStickOrders already maintain a
			// reasonably clean tree during the run.
			snap := snapshotTree(tssb, llh, iter, chainID)
			trees = append(trees, TreeSample{
				Iteration: iter,
				LLH:       llh,
				NumNodes:  len(nodes),
				Snapshot:  snap,
			})
		}

		// Progress logging
		if (iter+burnin+1)%100 == 0 {
			elapsed := time.Since(start)
			progress := float64(iter+burnin+1) / float64(totalIters) * 100
			log.Printf("Chain %d: iter=%d/%d (%.1f%%) llh=%.2f elapsed=%v",
				chainID, iter+burnin+1, totalIters, progress, llh, elapsed.Round(time.Second))
		}
	}

	return ChainResult{
		ChainID:     chainID,
		Trees:       trees,
		BurninLLH:   burninLLH,
		SampleLLH:   sampleLLH,
		FinalTree:   tssb,
		ElapsedTime: time.Since(start),
	}
}

// ============================================================================
// Output
// ============================================================================

// writeMutList writes mutlist.json to outDir, matching the Python
// write_results.py / result_generator.py output format.
// ssms and cnvs come from the best chain's FinalTree.
func writeMutList(outDir string, ssms []*SSM, cnvs []*CNV) error {
	// --- SSM section ---
	type ssmEntry struct {
		RefReads             []int   `json:"ref_reads"`
		TotalReads           []int   `json:"total_reads"`
		ExpectedRefInRef     float64 `json:"expected_ref_in_ref"`
		ExpectedRefInVariant float64 `json:"expected_ref_in_variant"`
		Name                 string  `json:"name,omitempty"`
	}
	ssmsOut := make(map[string]ssmEntry, len(ssms))
	for _, s := range ssms {
		ssmsOut[s.ID] = ssmEntry{
			RefReads:             s.A,
			TotalReads:           s.D,
			ExpectedRefInRef:     s.MuR,
			ExpectedRefInVariant: s.MuV,
			Name:                 s.Name,
		}
	}

	// --- CNV section ---
	type cnvSSMEntry struct {
		SSMID      string `json:"ssm_id"`
		MaternalCN int    `json:"maternal_cn"`
		PaternalCN int    `json:"paternal_cn"`
	}
	type cnvEntry struct {
		RefReads     []int         `json:"ref_reads"`
		TotalReads   []int         `json:"total_reads"`
		PhysicalCNVs []PhysicalCNV `json:"physical_cnvs"`
		SSMs         []cnvSSMEntry `json:"ssms"`
	}
	cnvsOut := make(map[string]cnvEntry, len(cnvs))
	for _, c := range cnvs {
		ssmLinks := make([]cnvSSMEntry, 0, len(c.SSMLinks))
		for _, link := range c.SSMLinks {
			ssmLinks = append(ssmLinks, cnvSSMEntry{
				SSMID:      link.SSMID,
				MaternalCN: link.MaternalCN,
				PaternalCN: link.PaternalCN,
			})
		}
		physCNVs := c.PhysicalCNVs
		if physCNVs == nil {
			physCNVs = []PhysicalCNV{} // emit [] not null
		}
		cnvsOut[c.ID] = cnvEntry{
			RefReads:     c.A,
			TotalReads:   c.D,
			PhysicalCNVs: physCNVs,
			SSMs:         ssmLinks,
		}
	}

	out := struct {
		SSMs map[string]ssmEntry `json:"ssms"`
		CNVs map[string]cnvEntry `json:"cnvs"`
	}{
		SSMs: ssmsOut,
		CNVs: cnvsOut,
	}

	data, err := json.MarshalIndent(out, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(filepath.Join(outDir, "mutlist.json"), data, 0644)
}

func writeResults(outDir string, results []ChainResult) error {
	if err := os.MkdirAll(outDir, 0755); err != nil {
		return err
	}

	// Write summary
	summaryPath := filepath.Join(outDir, "summary.json")
	summary := struct {
		NumChains  int      `json:"num_chains"`
		TotalTime  string   `json:"total_time"`
		ChainTimes []string `json:"chain_times"`
		BestLLH    float64  `json:"best_llh"`
		TreeCounts []int    `json:"trees_per_chain"`
	}{
		NumChains:  len(results),
		TreeCounts: make([]int, len(results)),
		ChainTimes: make([]string, len(results)),
	}

	bestLLH := math.Inf(-1)
	var maxTime time.Duration
	for i, r := range results {
		summary.TreeCounts[i] = len(r.Trees)
		summary.ChainTimes[i] = r.ElapsedTime.String()
		if r.ElapsedTime > maxTime {
			maxTime = r.ElapsedTime
		}
		for _, t := range r.Trees {
			if t.LLH > bestLLH {
				bestLLH = t.LLH
			}
		}
	}
	summary.TotalTime = maxTime.String()
	summary.BestLLH = bestLLH

	data, _ := json.MarshalIndent(summary, "", "  ")
	if err := os.WriteFile(summaryPath, data, 0644); err != nil {
		return err
	}

	// Write per-chain LLH trace (backward-compatible TSV format)
	for _, r := range results {
		chainPath := filepath.Join(outDir, fmt.Sprintf("chain_%d_samples.txt", r.ChainID))
		f, err := os.Create(chainPath)
		if err != nil {
			return err
		}
		fmt.Fprintln(f, "Iteration\tLLH\tNumNodes")
		for _, t := range r.Trees {
			fmt.Fprintf(f, "%d\t%.6f\t%d\n", t.Iteration, t.LLH, t.NumNodes)
		}
		f.Close()
	}

	// Write per-chain full posterior tree snapshots as newline-delimited JSON.
	// Each line is a self-contained JSON object with fields: iteration, chain_id,
	// llh, num_populations, populations, structure, mut_assignments.
	// This file replaces the Python original's trees.zip and can be consumed by
	// downstream post-processing tools (e.g. write_results.py with a Go-aware
	// adapter, or a future Go-native posterior analysis pipeline).
	for _, r := range results {
		ndjsonPath := filepath.Join(outDir, fmt.Sprintf("chain_%d_trees.ndjson", r.ChainID))
		f, err := os.Create(ndjsonPath)
		if err != nil {
			return err
		}
		w := bufio.NewWriter(f)
		for _, t := range r.Trees {
			if t.Snapshot != nil {
				w.Write(t.Snapshot)
				w.WriteByte('\n')
			}
		}
		if err := w.Flush(); err != nil {
			f.Close()
			return err
		}
		f.Close()
	}

	// Write best_tree.json from the chain with best final LLH
	bestChainIdx := 0
	bestOverallLLH := math.Inf(-1)
	for i, r := range results {
		for _, llhVal := range r.SampleLLH {
			if llhVal > bestOverallLLH {
				bestOverallLLH = llhVal
				bestChainIdx = i
			}
		}
	}
	if tssb := results[bestChainIdx].FinalTree; tssb != nil {
		removeEmptyNodes(tssb.Root, nil)
		treeSummary := summarizePops(tssb, bestOverallLLH, results[bestChainIdx].ChainID)
		treeData, _ := json.MarshalIndent(treeSummary, "", "  ")
		treePath := filepath.Join(outDir, "best_tree.json")
		if err := os.WriteFile(treePath, treeData, 0644); err != nil {
			return err
		}
		if err := writeMutList(outDir, tssb.Data, tssb.CNVData); err != nil {
			return err
		}
	}

	return nil
}

// removeEmptyNodes strips empty nodes from the TSSB tree.
// Matches Python's remove_empty_nodes (util2.py:128-155):
//   - Empty leaf: remove from parent
//   - Empty internal node: reparent children to grandparent
//   - Root is never removed
func removeEmptyNodes(root *TSSBNode, parent *TSSBNode) {
	// Process children depth-first (copy slice since we modify it)
	children := make([]*TSSBNode, len(root.Children))
	copy(children, root.Children)
	for _, child := range children {
		removeEmptyNodes(child, root)
	}

	if len(root.Node.Data) > 0 {
		return // Node has data, keep it
	}

	if parent == nil {
		return // Never remove root
	}

	// Find this node's index in parent
	idx := -1
	for i, c := range parent.Children {
		if c == root {
			idx = i
			break
		}
	}
	if idx < 0 {
		return // Already removed
	}

	if len(root.Children) == 0 {
		// Empty leaf: remove from parent
		parent.Children = append(parent.Children[:idx], parent.Children[idx+1:]...)
		if idx < len(parent.Sticks) {
			parent.Sticks = append(parent.Sticks[:idx], parent.Sticks[idx+1:]...)
		}
		// Remove from Node.Children
		for i, c := range parent.Node.Children {
			if c == root.Node {
				parent.Node.Children = append(parent.Node.Children[:i], parent.Node.Children[i+1:]...)
				break
			}
		}
	} else {
		// Empty internal node: reparent children to grandparent
		for i, child := range root.Children {
			parent.Children = append(parent.Children, child)
			if i < len(root.Sticks) {
				parent.Sticks = append(parent.Sticks, root.Sticks[i])
			}
			child.Node.Parent = parent.Node
			parent.Node.Children = append(parent.Node.Children, child.Node)
		}
		// Remove this node from parent
		parent.Children = append(parent.Children[:idx], parent.Children[idx+1:]...)
		if idx < len(parent.Sticks) {
			parent.Sticks = append(parent.Sticks[:idx], parent.Sticks[idx+1:]...)
		}
		for i, c := range parent.Node.Children {
			if c == root.Node {
				parent.Node.Children = append(parent.Node.Children[:i], parent.Node.Children[i+1:]...)
				break
			}
		}
	}
}

// summarizePops traverses the cleaned tree and extracts population summaries.
// Matches Python's ResultGenerator._summarize_pops (result_generator.py:43-88):
//   - Population 0 = root (normal cells)
//   - Children sorted by decreasing mean phi
//   - Preorder traversal assigns population IDs
func summarizePops(tssb *TSSB, llh float64, chainID int) map[string]interface{} {
	type popSummary struct {
		CellPrev []float64 `json:"cellular_prevalence"`
		NumSSMs  int       `json:"num_ssms"`
		NumCNVs  int       `json:"num_cnvs"`
	}
	type mutAssignment struct {
		SSMs []string `json:"ssms"`
		CNVs []string `json:"cnvs"`
	}

	pops := make(map[string]popSummary)
	structure := make(map[string][]int)
	mutAss := make(map[string]mutAssignment)
	popIdx := 0

	// Traverse TSSBNode tree (not Node.Children) since removeEmptyNodes
	// operates on TSSBNode.Children. This matches Python's result_generator.py
	// which traverses the TSSB dict structure (root['children']).
	var traverse func(tn *TSSBNode, parentIdx int)
	traverse = func(tn *TSSBNode, parentIdx int) {
		currentIdx := popIdx
		idStr := strconv.Itoa(currentIdx)
		node := tn.Node

		// Extract SSMs and CNVs
		ssms := []string{}
		cnvs := []string{}
		for _, ssmIdx := range node.Data {
			id := tssb.Data[ssmIdx].ID
			if len(id) > 0 && id[0] == 'c' {
				cnvs = append(cnvs, id)
			} else {
				ssms = append(ssms, id)
			}
		}

		// Sanitize phi values for JSON (replace Inf/NaN with 0)
		phi := make([]float64, len(node.Params))
		for k, v := range node.Params {
			if math.IsInf(v, 0) || math.IsNaN(v) {
				phi[k] = 0
			} else {
				phi[k] = v
			}
		}
		pops[idStr] = popSummary{
			CellPrev: phi,
			NumSSMs:  len(ssms),
			NumCNVs:  len(cnvs),
		}
		if len(ssms) > 0 || len(cnvs) > 0 {
			mutAss[idStr] = mutAssignment{SSMs: ssms, CNVs: cnvs}
		}

		// Sort children by decreasing mean phi (matching Python result_generator.py:81)
		sortedChildren := make([]*TSSBNode, len(tn.Children))
		copy(sortedChildren, tn.Children)
		sort.Slice(sortedChildren, func(i, j int) bool {
			meanI := 0.0
			for _, v := range sortedChildren[i].Node.Params {
				meanI += v
			}
			if len(sortedChildren[i].Node.Params) > 0 {
				meanI /= float64(len(sortedChildren[i].Node.Params))
			}
			meanJ := 0.0
			for _, v := range sortedChildren[j].Node.Params {
				meanJ += v
			}
			if len(sortedChildren[j].Node.Params) > 0 {
				meanJ /= float64(len(sortedChildren[j].Node.Params))
			}
			return meanI > meanJ
		})

		for _, child := range sortedChildren {
			popIdx++
			structure[idStr] = append(structure[idStr], popIdx)
			traverse(child, currentIdx)
		}
	}

	traverse(tssb.Root, -1)

	// Count only non-empty populations (with SSMs or CNVs), matching Python
	nonEmpty := 0
	for _, p := range pops {
		if p.NumSSMs > 0 || p.NumCNVs > 0 {
			nonEmpty++
		}
	}

	return map[string]interface{}{
		"chain_id":        chainID,
		"llh":             llh,
		"num_populations": nonEmpty, // non-empty pops only, matching Python convention
		"populations":     pops,
		"structure":       structure,
		"mut_assignments": mutAss,
	}
}

// snapshotTree encodes the current TSSB state as a JSON object suitable for
// one line of an all_trees.ndjson file.  It calls summarizePops directly (no
// removeEmptyNodes) because cullTree + resampleStickOrders have already run and
// the tree is sufficiently clean mid-run.  removeEmptyNodes is a post-processing
// step; skipping it here preserves the live tree for the next MCMC iteration.
func snapshotTree(tssb *TSSB, llh float64, iteration, chainID int) json.RawMessage {
	summary := summarizePops(tssb, llh, chainID)
	// Stamp iteration into the snapshot so each line is self-contained
	summary["iteration"] = iteration
	data, err := json.Marshal(summary)
	if err != nil {
		// This should not happen; summarizePops only uses serializable types.
		// Return nil so the caller skips this snapshot rather than crashing.
		return nil
	}
	return data
}

// ============================================================================
// Main
// ============================================================================

func main() {
	// Parse command line flags
	samples := flag.Int("s", 2500, "Number of MCMC samples")
	burnin := flag.Int("B", 1000, "Burn-in samples")
	chains := flag.Int("j", 6, "Number of parallel chains")
	outDir := flag.String("O", "output", "Output directory")
	mhIters := flag.Int("i", 5000, "MH iterations per MCMC iteration")
	seed := flag.Int64("r", 0, "Random seed (0 = use time)")
	noGPU := flag.Bool("no-gpu", false, "Disable GPU acceleration")
	cpuProfile := flag.String("cpuprofile", "", "Write CPU profile to file")

	flag.Parse()

	// CPU profiling
	if *cpuProfile != "" {
		f, err := os.Create(*cpuProfile)
		if err != nil {
			log.Fatalf("Could not create CPU profile: %v", err)
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatalf("Could not start CPU profile: %v", err)
		}
		defer pprof.StopCPUProfile()
	}
	args := flag.Args()

	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "Usage: phylowgs-go [options] ssm_data.txt cnv_data.txt\n")
		flag.PrintDefaults()
		os.Exit(1)
	}

	ssmFile := args[0]
	cnvFile := args[1]

	if *seed == 0 {
		*seed = time.Now().UnixNano()
	}

	log.Printf("PhyloWGS-Go starting")
	log.Printf("  SSM file: %s", ssmFile)
	log.Printf("  CNV file: %s", cnvFile)
	log.Printf("  Burn-in: %d, Samples: %d, Chains: %d", *burnin, *samples, *chains)
	log.Printf("  MH iterations: %d", *mhIters)
	log.Printf("  Random seed: %d", *seed)
	log.Printf("  CPUs: %d", runtime.NumCPU())

	// Load data
	ssms, err := loadSSMData(ssmFile)
	if err != nil {
		log.Fatalf("Failed to load SSM data: %v", err)
	}
	log.Printf("Loaded %d SSMs with %d timepoints", len(ssms), len(ssms[0].A))

	cnvs, err := loadCNVData(cnvFile, ssms)
	if err != nil {
		log.Fatalf("Failed to load CNV data: %v", err)
	}
	log.Printf("Loaded %d CNVs", len(cnvs))

	// Initialize GPU if available and not disabled
	if !*noGPU && cuda.Available() {
		nTimepoints := len(ssms[0].A)
		engine, err := cuda.NewGPUEngine(len(ssms), nTimepoints)
		if err != nil {
			log.Printf("Warning: GPU initialization failed: %v", err)
			log.Printf("Continuing with CPU-only mode")
		} else {
			gpuEngine = engine
			useGPU = true
			defer gpuEngine.Close()
			log.Printf("GPU acceleration enabled")
		}
	} else if *noGPU {
		log.Printf("GPU acceleration disabled by --no-gpu flag")
	} else {
		log.Printf("CUDA not available, using CPU-only mode")
	}

	// Run chains in parallel
	start := time.Now()
	results := make([]ChainResult, *chains)
	var wg sync.WaitGroup

	for i := 0; i < *chains; i++ {
		wg.Add(1)
		go func(chainID int) {
			defer wg.Done()
			chainSeed := *seed + int64(chainID)*1000
			results[chainID] = runChain(chainID, ssms, cnvs, *burnin, *samples, *mhIters, chainSeed)
		}(i)
	}

	wg.Wait()
	totalTime := time.Since(start)

	// Compute statistics
	totalSamples := 0
	var allLLH []float64
	for _, r := range results {
		totalSamples += len(r.Trees)
		allLLH = append(allLLH, r.SampleLLH...)
	}

	// Sort LLH to find best
	sort.Float64s(allLLH)
	bestLLH := allLLH[len(allLLH)-1]
	medianLLH := allLLH[len(allLLH)/2]

	log.Printf("\n=== Results ===")
	log.Printf("Total time: %v", totalTime.Round(time.Millisecond))
	log.Printf("Total samples: %d", totalSamples)
	log.Printf("Throughput: %.2f samples/sec/chain", float64(totalSamples)/totalTime.Seconds()/float64(*chains))
	log.Printf("Best LLH: %.2f", bestLLH)
	log.Printf("Median LLH: %.2f", medianLLH)

	// Write results
	if err := writeResults(*outDir, results); err != nil {
		log.Printf("Warning: failed to write results: %v", err)
	}
	log.Printf("Results written to %s", *outDir)
}
