package main

import (
	"math"
	"strconv"
	"strings"
)

// PyFloat marshals like Python 2.7's json.dumps(float) (repr): shortest
// round-trip digits, always a decimal point or exponent, so a value that is
// whole still reads back as a Python float (1.0, not 1).
//
// NaN and +/-Inf make MarshalJSON return Python's tokens (NaN, Infinity,
// -Infinity), but encoding/json validates MarshalJSON output and rejects
// them, so any NaN/Inf PyFloat makes the enclosing json.Marshal return an
// error. This is intended: a NaN reaching output is an upstream bug and must
// fail loudly.
type PyFloat float64

// MarshalJSON implements json.Marshaler.
func (f PyFloat) MarshalJSON() ([]byte, error) {
	v := float64(f)
	switch {
	case math.IsNaN(v):
		return []byte("NaN"), nil
	case math.IsInf(v, 1):
		return []byte("Infinity"), nil
	case math.IsInf(v, -1):
		return []byte("-Infinity"), nil
	}
	e := strconv.FormatFloat(v, 'e', -1, 64) // d.ddde±XX
	if i := strings.IndexByte(e, 'e'); i >= 0 {
		if exp, err := strconv.Atoi(e[i+1:]); err == nil && (exp < -4 || exp >= 16) {
			return []byte(e), nil
		}
	}
	s := strconv.FormatFloat(v, 'f', -1, 64)
	if !strings.Contains(s, ".") {
		s += ".0"
	}
	return []byte(s), nil
}

// pyFloats converts a float64 slice to PyFloats (nil stays an empty slice so
// it marshals as [] rather than null).
func pyFloats(xs []float64) []PyFloat {
	out := make([]PyFloat, len(xs))
	for i, x := range xs {
		out[i] = PyFloat(x)
	}
	return out
}

// retypeSummary wraps llh and every populations[*].cellular_prevalence element
// of a JSON-decoded summary in PyFloat, in place, so re-marshalling keeps the
// Python float class. Already-typed values are left alone.
func retypeSummary(generic map[string]interface{}) {
	if v, ok := generic["llh"].(float64); ok {
		generic["llh"] = PyFloat(v)
	}
	pops, _ := generic["populations"].(map[string]interface{})
	for _, p := range pops {
		pop, ok := p.(map[string]interface{})
		if !ok {
			continue
		}
		cp, ok := pop["cellular_prevalence"].([]interface{})
		if !ok {
			continue
		}
		for i, x := range cp {
			if v, ok := x.(float64); ok {
				cp[i] = PyFloat(v)
			}
		}
	}
}
