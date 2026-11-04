package main

// Go port of pwgsresults/result_munger.py, part 2: mutation reassignment
// after node removal, remove_superclones and remove_multiprimary_trees.
// Python's progress print()s are not reproduced (stdout only, not results).

import (
	"fmt"
	"math"
)

// reassignMuts ports _reassign_muts: walk the idx map ascending, collecting
// removed nodes' MutAss and moving renumbered ones down, then send the
// collected mutations to their destination.
func (rs *ResultSet) reassignMuts(tidx int, idxMap map[int]int, dest mutDest) error {
	mutass := rs.MutAss[tidx]
	var deleted []*MutAss
	for _, s := range sortedKeys(idxMap) {
		n := idxMap[s]
		if n != nodeRemoved && n >= s {
			return fmt.Errorf("tree %d: renumbering %d to %d does not move down", tidx, s, n)
		}
		m, ok := mutass[s]
		if !ok {
			return fmt.Errorf("tree %d: node %d has no mut_assignments entry", tidx, s)
		}
		if n == nodeRemoved {
			deleted = append(deleted, m)
		} else {
			mutass[n] = m
		}
		delete(mutass, s)
	}
	for _, dm := range deleted {
		var err error
		if dest == destBest {
			err = rs.moveMutsToBestNode(tidx, dm)
		} else {
			err = rs.moveMutsToClonalNode(tidx, dm)
		}
		if err != nil {
			return err
		}
	}
	return nil
}

// mutReads returns (ref_reads, total_reads) for an SSM or CNV from Muts.
func (rs *ResultSet) mutReads(isCNV bool, id string) ([]int, []int, error) {
	if isCNV {
		r, ok := rs.Muts.CNVReads[id]
		if !ok {
			return nil, nil, fmt.Errorf("cnv %s not in mutlist", id)
		}
		return r[0], r[1], nil
	}
	s, ok := rs.Muts.SSMs[id]
	if !ok {
		return nil, nil, fmt.Errorf("ssm %s not in mutlist", id)
	}
	return s.RefReads, s.TotalReads, nil
}

// bestNode ports the inner loop of _move_muts_to_best_node: ascending node
// order, strict `<` (ties keep the lower index), lowest starts at 1, root
// skipped. Returns -1 when no node is within phi 1.
func bestNode(pops map[int]*Population, impliedPhi float64) int {
	lowest, best := 1.0, -1
	for _, p := range sortedKeys(pops) {
		delta := math.Abs(pyMean(pops[p].CellPrev) - impliedPhi)
		if delta < lowest && p != 0 {
			lowest, best = delta, p
		}
	}
	return best
}

// moveMutsToBestNode ports _move_muts_to_best_node (SSMs first, then CNVs).
func (rs *ResultSet) moveMutsToBestNode(tidx int, dm *MutAss) error {
	mutass, pops := rs.MutAss[tidx], rs.Trees[tidx].Populations
	for _, isCNV := range []bool{false, true} {
		ids := dm.SSMs
		if isCNV {
			ids = dm.CNVs
		}
		for _, id := range ids {
			ref, total, err := rs.mutReads(isCNV, id)
			if err != nil {
				return fmt.Errorf("tree %d: %w", tidx, err)
			}
			r, t := pyMean(ref), pyMean(total)
			phi := 2 * (t - r) / t
			if 1.0 < phi { // Python min(phi, 1.0): keeps phi unless 1.0 < phi (so NaN stays NaN)
				phi = 1.0
			}
			best := bestNode(pops, phi)
			if best < 0 {
				return fmt.Errorf("tree %d: no node within phi 1 of mutation %s", tidx, id)
			}
			target, ok := mutass[best]
			if !ok { // Python KeyError on mutass[best_node]
				return fmt.Errorf("tree %d: best node %d for %s has no mut_assignments entry", tidx, best, id)
			}
			if isCNV {
				target.CNVs = appendOwned(target.CNVs, id)
			} else {
				target.SSMs = appendOwned(target.SSMs, id)
			}
		}
	}
	return nil
}

// moveMutsToClonalNode ports _move_muts_to_clonal_node.
func (rs *ResultSet) moveMutsToClonalNode(tidx int, dm *MutAss) error {
	rootCh, ok := rs.Trees[tidx].Structure[0]
	if !ok || !containsInt(rootCh, 1) {
		return fmt.Errorf("tree %d: clonal node 1 is not a child of the root", tidx)
	}
	clonal, ok := rs.MutAss[tidx][1]
	if !ok {
		return fmt.Errorf("tree %d: clonal node 1 has no mut_assignments entry", tidx)
	}
	clonal.SSMs = appendOwned(clonal.SSMs, dm.SSMs...)
	clonal.CNVs = appendOwned(clonal.CNVs, dm.CNVs...)
	return nil
}

// RemoveSuperclones ports ResultMunger.remove_superclones.
func (rs *ResultSet) RemoveSuperclones() error {
	for _, tidx := range sortedKeys(rs.Trees) {
		if err := rs.removeSuperclone(tidx); err != nil {
			return err
		}
	}
	return nil
}

func (rs *ResultSet) removeSuperclone(tidx int) error {
	ts := rs.Trees[tidx]
	rootCh := ts.Structure[0]
	if len(rootCh) == 0 {
		return fmt.Errorf("tree %d: root has no children", tidx)
	}
	if len(rootCh) != 1 {
		return nil // multiprimary: ignored
	}
	if rootCh[0] != 1 {
		return fmt.Errorf("tree %d: clonal node is %d, expected 1", tidx, rootCh[0])
	}
	if ch, ok := ts.Structure[1]; !ok || len(ch) != 1 {
		return nil // no subclones, or several
	}
	if c := ts.Structure[1][0]; c != 2 {
		return fmt.Errorf("tree %d: clonal node's child is %d, expected 2", tidx, c)
	}
	clone, child := ts.Populations[1], ts.Populations[2]
	if clone == nil || child == nil {
		return fmt.Errorf("tree %d: missing population 1 or 2", tidx)
	}
	if child.NumSSMs == 0 {
		return nil
	}
	if !(float64(clone.NumSSMs)/float64(child.NumSSMs) <= 0.33) {
		return nil
	}
	if len(clone.CellPrev) != len(child.CellPrev) {
		return fmt.Errorf("tree %d: cellular_prevalence lengths differ (%d vs %d)", tidx, len(clone.CellPrev), len(child.CellPrev))
	}
	diff := make([]float64, len(clone.CellPrev))
	merged := make([]float64, len(clone.CellPrev))
	n1, n2 := float64(clone.NumSSMs), float64(child.NumSSMs)
	for i := range diff {
		diff[i] = math.Abs(clone.CellPrev[i] - child.CellPrev[i])
		merged[i] = (clone.CellPrev[i]*n1 + child.CellPrev[i]*n2) / float64(clone.NumSSMs+child.NumSSMs)
	}
	if !(pyMean(diff) <= 0.1) {
		return nil
	}
	clone.CellPrev = merged
	return rs.removeNodes([]int{2}, tidx, destClonal)
}

// RemoveMultiprimaryTrees ports ResultMunger.remove_multiprimary_trees.
//
// Quirk kept: Python deletes the multiprimary trees *before* its renumbering
// loop, so `tidx in polyidxs` is never true inside it and the loop is a
// no-op. Surviving trees therefore keep their original indices (e.g. 0,2,4),
// and this port does no renumbering.
func (rs *ResultSet) RemoveMultiprimaryTrees(maxProp float64) error {
	var poly []int
	for _, tidx := range sortedKeys(rs.Trees) {
		rootCh := rs.Trees[tidx].Structure[0]
		if len(rootCh) == 0 {
			return fmt.Errorf("tree %d: root has no children", tidx)
		}
		if len(rootCh) != 1 {
			poly = append(poly, tidx)
		}
	}
	total := len(rs.Trees)
	if total == 0 {
		return fmt.Errorf("no trees: float division by zero")
	}
	frac := float64(len(poly)) / float64(total)
	if frac >= maxProp {
		return fmt.Errorf("%d%% of trees are multiprimary (%d of %d), so not enough to report good posterior.", int(100*frac), len(poly), total)
	}
	for _, p := range poly {
		if _, ok := rs.MutAss[p]; !ok {
			return fmt.Errorf("tree %d: no mut_assignments to delete", p)
		}
		delete(rs.Trees, p)
		delete(rs.MutAss, p)
	}
	if len(rs.Trees) != len(rs.MutAss) {
		return fmt.Errorf("tree and mut_assignments indices differ after multiprimary removal")
	}
	for tidx := range rs.Trees {
		if _, ok := rs.MutAss[tidx]; !ok {
			return fmt.Errorf("tree %d has no mut_assignments after multiprimary removal", tidx)
		}
	}
	return nil
}
