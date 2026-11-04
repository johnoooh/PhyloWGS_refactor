package main

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"testing"
)

const maxCompareDiffs = 20

// isFloatLiteral reports whether a JSON number literal is a Python float
// (has '.', 'e', 'E', or is NaN/Infinity); otherwise it is a Python int.
func isFloatLiteral(lit string) bool {
	return strings.ContainsAny(lit, ".eE") || strings.Contains(lit, "NaN") || strings.Contains(lit, "Infinity")
}

func decodeUseNumber(b []byte) (interface{}, error) {
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.UseNumber()
	var v interface{}
	if err := dec.Decode(&v); err != nil {
		return nil, err
	}
	if dec.More() {
		return nil, fmt.Errorf("trailing data after JSON value")
	}
	return v, nil
}

type jsonDiffer struct {
	diffs []string
}

func (d *jsonDiffer) add(path []string, format string, args ...interface{}) {
	if len(d.diffs) < maxCompareDiffs {
		d.diffs = append(d.diffs, "/"+strings.Join(path, "/")+": "+fmt.Sprintf(format, args...))
	}
}

func jsonKind(v interface{}) string {
	switch v.(type) {
	case map[string]interface{}:
		return "object"
	case []interface{}:
		return "array"
	case string:
		return "string"
	case json.Number:
		return "number"
	case bool:
		return "bool"
	case nil:
		return "null"
	}
	return fmt.Sprintf("%T", v)
}

// isMutList reports whether path is .../mut_assignments/<node>/ssms|cnvs.
func isMutList(path []string) bool {
	n := len(path)
	return n >= 3 && path[n-3] == "mut_assignments" && (path[n-1] == "ssms" || path[n-1] == "cnvs")
}

func sortedStrings(a []interface{}) []interface{} {
	out := append([]interface{}(nil), a...)
	sort.SliceStable(out, func(i, j int) bool { return fmt.Sprint(out[i]) < fmt.Sprint(out[j]) })
	return out
}

func (d *jsonDiffer) number(path []string, got, want json.Number) {
	gl, wl := got.String(), want.String()
	if isFloatLiteral(gl) != isFloatLiteral(wl) {
		d.add(path, "number class differs: got %s, want %s", gl, wl)
		return
	}
	if !isFloatLiteral(gl) {
		if gl != wl {
			d.add(path, "int differs: got %s, want %s", gl, wl)
		}
		return
	}
	g, err1 := strconv.ParseFloat(gl, 64)
	w, err2 := strconv.ParseFloat(wl, 64)
	if err1 != nil || err2 != nil {
		d.add(path, "unparseable float: got %s, want %s", gl, wl)
		return
	}
	if g == w {
		return
	}
	if len(path) > 0 && path[0] == "tree_densities" {
		if math.Abs(g-w) <= 1e-9*math.Max(math.Abs(g), math.Abs(w)) {
			return
		}
	}
	d.add(path, "float differs: got %s, want %s", gl, wl)
}

func (d *jsonDiffer) walk(path []string, got, want interface{}) {
	if jsonKind(got) != jsonKind(want) {
		d.add(path, "type differs: got %s, want %s", jsonKind(got), jsonKind(want))
		return
	}
	switch w := want.(type) {
	case map[string]interface{}:
		g := got.(map[string]interface{})
		keys := make([]string, 0, len(w)+len(g))
		for k := range w {
			keys = append(keys, k)
		}
		for k := range g {
			if _, ok := w[k]; !ok {
				keys = append(keys, k)
			}
		}
		sort.Strings(keys)
		for _, k := range keys {
			sub := append(append([]string(nil), path...), k)
			gv, gok := g[k]
			wv, wok := w[k]
			switch {
			case !gok:
				d.add(sub, "missing key (present in want)")
			case !wok:
				d.add(sub, "extra key (absent in want)")
			default:
				d.walk(sub, gv, wv)
			}
		}
	case []interface{}:
		g := got.([]interface{})
		if isMutList(path) {
			g, w = sortedStrings(g), sortedStrings(w)
		}
		if len(g) != len(w) {
			d.add(path, "array length differs: got %d, want %d", len(g), len(w))
			return
		}
		for i := range w {
			d.walk(append(append([]string(nil), path...), strconv.Itoa(i)), g[i], w[i])
		}
	case json.Number:
		d.number(path, got.(json.Number), w)
	default:
		if got != want {
			d.add(path, "value differs: got %v, want %v", got, want)
		}
	}
}

// diffPythonJSON returns up to maxCompareDiffs differences between two JSON
// documents, treating int vs float literals as different (Python number classes).
func diffPythonJSON(got, want []byte) []string {
	g, err := decodeUseNumber(got)
	if err != nil {
		return []string{"got is not valid JSON: " + err.Error()}
	}
	w, err := decodeUseNumber(want)
	if err != nil {
		return []string{"want is not valid JSON: " + err.Error()}
	}
	d := &jsonDiffer{}
	d.walk(nil, g, w)
	return d.diffs
}

// comparePythonJSON fails t with the first differences between got and want.
func comparePythonJSON(t *testing.T, label string, got, want []byte) {
	t.Helper()
	if diffs := diffPythonJSON(got, want); len(diffs) > 0 {
		t.Errorf("%s: %d difference(s):\n  %s", label, len(diffs), strings.Join(diffs, "\n  "))
	}
}

func TestComparePythonJSON_Kills(t *testing.T) {
	cases := []struct {
		name, got, want string
		wantDiff        string // substring of the first diff; "" means no diff expected
	}{
		{"int vs float class", `{"a":1}`, `{"a":1.0}`, "/a: number class differs"},
		{"float vs int class", `{"a":1.0}`, `{"a":1}`, "/a: number class differs"},
		{"exponent is float", `{"a":1e0}`, `{"a":1}`, "number class differs"},
		{"int value", `{"a":1}`, `{"a":2}`, "/a: int differs"},
		{"float exact", `{"a":0.1}`, `{"a":0.1000001}`, "/a: float differs"},
		{"densities within tol", `{"tree_densities":{"0":1.0000000001}}`, `{"tree_densities":{"0":1.0}}`, ""},
		{"densities outside tol", `{"tree_densities":{"0":1.001}}`, `{"tree_densities":{"0":1.0}}`, "/tree_densities/0"},
		{"same pair outside densities", `{"x":{"0":1.0000000001}}`, `{"x":{"0":1.0}}`, "/x/0"},
		{"ssms order ignored", `{"mut_assignments":{"1":{"ssms":["s1","s0"]}}}`, `{"mut_assignments":{"1":{"ssms":["s0","s1"]}}}`, ""},
		{"cnvs order ignored", `{"trees":{"0":{"mut_assignments":{"1":{"cnvs":["c1","c0"]}}}}}`, `{"trees":{"0":{"mut_assignments":{"1":{"cnvs":["c0","c1"]}}}}}`, ""},
		{"other list order matters", `{"structure":{"0":[2,1]}}`, `{"structure":{"0":[1,2]}}`, "/structure/0/0"},
		{"ssms content differs", `{"mut_assignments":{"1":{"ssms":["s1","s2"]}}}`, `{"mut_assignments":{"1":{"ssms":["s0","s1"]}}}`, "ssms"},
		{"missing key", `{"a":1}`, `{"a":1,"b":2}`, "/b: missing key"},
		{"extra key", `{"a":1,"b":2}`, `{"a":1}`, "/b: extra key"},
		{"type differs", `{"a":"1"}`, `{"a":1}`, "/a: type differs"},
		{"list length", `{"a":[1]}`, `{"a":[1,2]}`, "/a: array length"},
		{"identical", `{"a":[1,2.5,{"b":null,"c":true,"d":"x"}]}`, `{"a":[1,2.5,{"b":null,"c":true,"d":"x"}]}`, ""},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			diffs := diffPythonJSON([]byte(c.got), []byte(c.want))
			if c.wantDiff == "" {
				if len(diffs) != 0 {
					t.Fatalf("expected no diff, got %v", diffs)
				}
				return
			}
			if len(diffs) == 0 || !strings.Contains(diffs[0], c.wantDiff) {
				t.Fatalf("expected first diff containing %q, got %v", c.wantDiff, diffs)
			}
		})
	}
}

func TestComparePythonJSON_CapsAndRejectsGarbage(t *testing.T) {
	var got, want []string
	for i := 0; i < 50; i++ {
		got = append(got, fmt.Sprintf(`"k%d":1`, i))
		want = append(want, fmt.Sprintf(`"k%d":1.0`, i))
	}
	diffs := diffPythonJSON([]byte("{"+strings.Join(got, ",")+"}"), []byte("{"+strings.Join(want, ",")+"}"))
	if len(diffs) != maxCompareDiffs {
		t.Fatalf("want %d diffs (cap), got %d", maxCompareDiffs, len(diffs))
	}
	if d := diffPythonJSON([]byte(`{"a":`), []byte(`{}`)); len(d) != 1 || !strings.Contains(d[0], "not valid JSON") {
		t.Fatalf("garbage input must be reported, got %v", d)
	}
}

// errPrefixRE matches an exception-class prefix that error.txt must not carry.
var errPrefixRE = regexp.MustCompile(`^\w[\w.]*: `)

// TestNativeResultsFixtures checks every case x flag set has the expected
// files, that each expected JSON file is self-consistent under the comparator,
// and that the known error cases carry the reference's message.
func TestNativeResultsFixtures(t *testing.T) {
	root := "testdata/native_results"
	cases := map[string][]string{
		"linear": {"default", "prod"}, "identical": {"default", "prod"},
		"single": {"default", "prod"}, "superclone": {"default", "prod"},
		"cnvonly": {"default", "prod"}, "multiprimary": {"default", "prod"},
		"sim": {"default", "prod", "keepsc", "mincount", "nonames"},
	}
	errCases := map[string]string{
		"multiprimary/default": "of trees are multiprimary",
		"single/default":       "`dataset` input should have multiple elements.",
		"single/prod":          "`dataset` input should have multiple elements.",
	}
	for c, sets := range cases {
		for _, fs := range sets {
			dir := root + "/" + c + "/" + fs + "/expected/"
			if want, isErr := errCases[c+"/"+fs]; isErr {
				b, err := os.ReadFile(dir + "error.txt")
				if err != nil || !strings.Contains(string(b), want) || errPrefixRE.Match(b) {
					t.Errorf("%s/%s: bad error.txt %q (%v)", c, fs, b, err)
				}
				continue
			}
			for _, f := range []string{"summ.json", "muts.json", "mutass.json"} {
				b, err := os.ReadFile(dir + f)
				if err != nil || len(b) == 0 {
					t.Errorf("%s/%s/%s: missing or empty (%v)", c, fs, f, err)
					continue
				}
				comparePythonJSON(t, c+"/"+fs+"/"+f, b, b)
			}
		}
	}
}
