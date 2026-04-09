# Local trace baseline report

Generated from `PhyloWGS_refactor/sim_validation/local_trace_runs`

Single-page calibration of Go-port spawn dynamics and K distributions
on a 4-fixture local slice, run BEFORE HPC results land. Used to
interpret tomorrow's HPC numbers.

## Python reference K distributions (same fixtures)

```
K3_S1_T200_M30_C2_rep0 N= 1000  range=4-16  median= 13.0  mean=11.95  mode=(13, 279)  best K=5 llh=-315.80
K3_S1_T200_M30_C0_rep0 N= 1000  range=3- 5  median=  3.0  mean=3.19  mode=(3, 820)  best K=3 llh=-108.04
K5_S1_T200_M50_C2_rep0 N= 1000  range=4- 6  median=  4.0  mean=4.13  mode=(4, 880)  best K=4 llh=-6649.55
K5_S1_T200_M50_C0_rep0 N= 1000  range=4-11  median=  6.0  mean=6.20  mode=(7, 273)  best K=4 llh=-212.87
```

## expA_baseline — HPC-parity baseline (B=500, S=1000, 1 chain)

### K3_S1_T200_M30_C2_rep0

**trace: `trace.ndjson`** total_iters=1500  burnin=500  sample=1000

```
  [burnin]
    K_before_cull   {'n': 500, 'min': 15, 'max': 208, 'mean': 83.254, 'median': 76.0, 'stdev': 34.648}
    K_after_cull    {'n': 500, 'min': 8, 'max': 103, 'mean': 41.122, 'median': 39.5, 'stdev': 14.81}
    spawns_find     {'n': 500, 'min': 3, 'max': 151, 'mean': 42.726, 'median': 38.0, 'stdev': 24.859}
    spawns_stick    {'n': 500, 'min': 0, 'max': 45, 'mean': 13.64, 'median': 13.0, 'stdev': 8.294}
    kills_in_cull   {'n': 500, 'min': 3, 'max': 142, 'mean': 42.132, 'median': 37.0, 'stdev': 24.267}
    cap_hits_find   0
    dp_alpha        {'n': 500, 'min': 3.5141111600835058, 'max': 49.981678963890204, 'mean': 25.371, 'median': 24.021915094887888, 'stdev': 12.171}
    dp_gamma        {'n': 500, 'min': 1.0001034692134678, 'max': 1.3966908797258486, 'mean': 1.074, 'median': 1.0513263271167137, 'stdev': 0.068}
    alpha_decay     {'n': 500, 'min': 0.22827609706134153, 'max': 0.786266318416499, 'mean': 0.475, 'median': 0.4882761468339276, 'stdev': 0.118}
    best_llh        -337.15
  [sample]
    K_before_cull   {'n': 1000, 'min': 18, 'max': 452, 'mean': 119.929, 'median': 106.0, 'stdev': 61.394}
    K_after_cull    {'n': 1000, 'min': 14, 'max': 137, 'mean': 53.383, 'median': 50.0, 'stdev': 20.51}
    spawns_find     {'n': 1000, 'min': 1, 'max': 340, 'mean': 67.039, 'median': 56.0, 'stdev': 45.234}
    spawns_stick    {'n': 1000, 'min': 0, 'max': 89, 'mean': 18.868, 'median': 17.0, 'stdev': 11.909}
    kills_in_cull   {'n': 1000, 'min': 4, 'max': 344, 'mean': 66.546, 'median': 56.0, 'stdev': 45.51}
    cap_hits_find   0
    dp_alpha        {'n': 1000, 'min': 3.3239574655220703, 'max': 49.92083861911006, 'mean': 23.513, 'median': 20.79325088813018, 'stdev': 12.112}
    dp_gamma        {'n': 1000, 'min': 1.000038267828546, 'max': 1.5180959117559074, 'mean': 1.079, 'median': 1.0531526985744488, 'stdev': 0.082}
    alpha_decay     {'n': 1000, 'min': 0.22660384463637867, 'max': 0.7997933919715973, 'mean': 0.558, 'median': 0.5722985365352633, 'stdev': 0.138}
    best_llh        -347.40
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=4-16  median= 13.0  mean=11.95  mode=(13, 279)  best K=5 llh=-315.80
Go local           N= 1000  range=6-15  median= 12.0  mean=11.55  mode=(12, 332)  best K=7 llh=-347.40
```

### K3_S1_T200_M30_C0_rep0

**trace: `trace.ndjson`** total_iters=1500  burnin=500  sample=1000

```
  [burnin]
    K_before_cull   {'n': 500, 'min': 2, 'max': 42, 'mean': 6.086, 'median': 5.0, 'stdev': 3.763}
    K_after_cull    {'n': 500, 'min': 2, 'max': 7, 'mean': 3.186, 'median': 3.0, 'stdev': 0.623}
    spawns_find     {'n': 500, 'min': 0, 'max': 40, 'mean': 2.862, 'median': 2.0, 'stdev': 3.645}
    spawns_stick    {'n': 500, 'min': 0, 'max': 7, 'mean': 0.214, 'median': 0.0, 'stdev': 0.857}
    kills_in_cull   {'n': 500, 'min': 0, 'max': 40, 'mean': 2.9, 'median': 2.0, 'stdev': 3.671}
    cap_hits_find   0
    dp_alpha        {'n': 500, 'min': 1.3207086600414633, 'max': 49.94040061029808, 'mean': 30.129, 'median': 32.43174867340839, 'stdev': 12.704}
    dp_gamma        {'n': 500, 'min': 1.0000243929743038, 'max': 1.3804899988511186, 'mean': 1.071, 'median': 1.0533087784859236, 'stdev': 0.064}
    alpha_decay     {'n': 500, 'min': 0.05005476682712496, 'max': 0.625763016019053, 'mean': 0.116, 'median': 0.09567905382504903, 'stdev': 0.072}
    best_llh        -107.74
  [sample]
    K_before_cull   {'n': 1000, 'min': 3, 'max': 26, 'mean': 5.578, 'median': 4.0, 'stdev': 3.212}
    K_after_cull    {'n': 1000, 'min': 3, 'max': 10, 'mean': 3.235, 'median': 3.0, 'stdev': 0.66}
    spawns_find     {'n': 1000, 'min': 0, 'max': 21, 'mean': 2.306, 'median': 1.0, 'stdev': 3.02}
    spawns_stick    {'n': 1000, 'min': 0, 'max': 8, 'mean': 0.158, 'median': 0.0, 'stdev': 0.71}
    kills_in_cull   {'n': 1000, 'min': 0, 'max': 21, 'mean': 2.343, 'median': 1.0, 'stdev': 3.044}
    cap_hits_find   0
    dp_alpha        {'n': 1000, 'min': 1.4680853465783528, 'max': 49.977949543571754, 'mean': 31.237, 'median': 32.049535471632524, 'stdev': 12.052}
    dp_gamma        {'n': 1000, 'min': 1.0000412986422078, 'max': 1.46393605922378, 'mean': 1.066, 'median': 1.0441472576467614, 'stdev': 0.069}
    alpha_decay     {'n': 1000, 'min': 0.050211262769863936, 'max': 0.4117078624905886, 'mean': 0.104, 'median': 0.09036205770685495, 'stdev': 0.05}
    best_llh        -107.73
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=3- 5  median=  3.0  mean=3.19  mode=(3, 820)  best K=3 llh=-108.04
Go local           N= 1000  range=2- 4  median=  2.0  mean=2.14  mode=(2, 871)  best K=2 llh=-107.73
```

### K5_S1_T200_M50_C2_rep0

**trace: `trace.ndjson`** total_iters=1500  burnin=500  sample=1000

```
  [burnin]
    K_before_cull   {'n': 500, 'min': 4, 'max': 48, 'mean': 8.266, 'median': 7.0, 'stdev': 4.273}
    K_after_cull    {'n': 500, 'min': 3, 'max': 10, 'mean': 4.468, 'median': 4.0, 'stdev': 0.882}
    spawns_find     {'n': 500, 'min': 0, 'max': 38, 'mean': 3.77, 'median': 3.0, 'stdev': 3.977}
    spawns_stick    {'n': 500, 'min': 0, 'max': 11, 'mean': 0.602, 'median': 0.0, 'stdev': 1.427}
    kills_in_cull   {'n': 500, 'min': 0, 'max': 38, 'mean': 3.798, 'median': 3.0, 'stdev': 3.909}
    cap_hits_find   0
    dp_alpha        {'n': 500, 'min': 1.3355423815057048, 'max': 49.45613219779591, 'mean': 21.57, 'median': 19.426259938436424, 'stdev': 12.459}
    dp_gamma        {'n': 500, 'min': 1.0000405356878534, 'max': 1.380711822826332, 'mean': 1.068, 'median': 1.0478274984405682, 'stdev': 0.068}
    alpha_decay     {'n': 500, 'min': 0.05002816122247102, 'max': 0.6189720446264673, 'mean': 0.119, 'median': 0.09820683427995033, 'stdev': 0.077}
    best_llh        -6648.68
  [sample]
    K_before_cull   {'n': 1000, 'min': 4, 'max': 36, 'mean': 8.786, 'median': 7.5, 'stdev': 4.695}
    K_after_cull    {'n': 1000, 'min': 4, 'max': 11, 'mean': 4.492, 'median': 4.0, 'stdev': 0.879}
    spawns_find     {'n': 1000, 'min': 0, 'max': 31, 'mean': 4.258, 'median': 3.0, 'stdev': 4.369}
    spawns_stick    {'n': 1000, 'min': 0, 'max': 9, 'mean': 0.483, 'median': 0.0, 'stdev': 1.193}
    kills_in_cull   {'n': 1000, 'min': 0, 'max': 31, 'mean': 4.294, 'median': 3.0, 'stdev': 4.402}
    cap_hits_find   0
    dp_alpha        {'n': 1000, 'min': 1.1196799992208613, 'max': 49.87148925459602, 'mean': 22.818, 'median': 21.50719062188626, 'stdev': 13.552}
    dp_gamma        {'n': 1000, 'min': 1.0000296680015028, 'max': 1.446454669801181, 'mean': 1.065, 'median': 1.0454378279648733, 'stdev': 0.065}
    alpha_decay     {'n': 1000, 'min': 0.05003484329775853, 'max': 0.718396451038849, 'mean': 0.129, 'median': 0.10004737005977861, 'stdev': 0.086}
    best_llh        -6648.56
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=4- 6  median=  4.0  mean=4.13  mode=(4, 880)  best K=4 llh=-6649.55
Go local           N= 1000  range=3- 5  median=  3.0  mean=3.17  mode=(3, 845)  best K=3 llh=-6648.56
```

### K5_S1_T200_M50_C0_rep0

**trace: `trace.ndjson`** total_iters=1500  burnin=500  sample=1000

```
  [burnin]
    K_before_cull   {'n': 500, 'min': 4, 'max': 124, 'mean': 30.426, 'median': 29.0, 'stdev': 15.046}
    K_after_cull    {'n': 500, 'min': 4, 'max': 50, 'mean': 15.736, 'median': 15.0, 'stdev': 6.738}
    spawns_find     {'n': 500, 'min': 0, 'max': 90, 'mean': 14.982, 'median': 13.0, 'stdev': 10.818}
    spawns_stick    {'n': 500, 'min': 0, 'max': 21, 'mean': 4.51, 'median': 4.0, 'stdev': 4.145}
    kills_in_cull   {'n': 500, 'min': 0, 'max': 95, 'mean': 14.69, 'median': 12.0, 'stdev': 11.215}
    cap_hits_find   0
    dp_alpha        {'n': 500, 'min': 1.052992488035015, 'max': 49.23969773534828, 'mean': 8.444, 'median': 5.411237296217706, 'stdev': 7.925}
    dp_gamma        {'n': 500, 'min': 1.0000912406301123, 'max': 1.4742559606076342, 'mean': 1.066, 'median': 1.0452487070545344, 'stdev': 0.07}
    alpha_decay     {'n': 500, 'min': 0.06883110385454057, 'max': 0.7946942139261521, 'mean': 0.441, 'median': 0.42975663696669275, 'stdev': 0.182}
    best_llh        -214.13
  [sample]
    K_before_cull   {'n': 1000, 'min': 4, 'max': 131, 'mean': 31.817, 'median': 28.0, 'stdev': 15.278}
    K_after_cull    {'n': 1000, 'min': 4, 'max': 49, 'mean': 16.594, 'median': 16.0, 'stdev': 6.361}
    spawns_find     {'n': 1000, 'min': 0, 'max': 84, 'mean': 15.749, 'median': 13.0, 'stdev': 11.443}
    spawns_stick    {'n': 1000, 'min': 0, 'max': 30, 'mean': 4.29, 'median': 3.0, 'stdev': 4.277}
    kills_in_cull   {'n': 1000, 'min': 0, 'max': 103, 'mean': 15.223, 'median': 13.0, 'stdev': 11.51}
    cap_hits_find   0
    dp_alpha        {'n': 1000, 'min': 1.2444777095599036, 'max': 49.98254432953915, 'mean': 21.802, 'median': 20.49164725436067, 'stdev': 14.313}
    dp_gamma        {'n': 1000, 'min': 1.0000600413533323, 'max': 1.4258258629000313, 'mean': 1.069, 'median': 1.0464332954524278, 'stdev': 0.067}
    alpha_decay     {'n': 1000, 'min': 0.07026322252433453, 'max': 0.7985676361033074, 'mean': 0.331, 'median': 0.2920163773092378, 'stdev': 0.157}
    best_llh        -215.05
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=4-11  median=  6.0  mean=6.20  mode=(7, 273)  best K=4 llh=-212.87
Go local           N= 1000  range=3-11  median=  7.0  mean=6.74  mode=(7, 353)  best K=5 llh=-215.05
```

## expB_long — Long convergence (B=2000, S=5000, 1 chain)

### K3_S1_T200_M30_C2_rep0

**trace: `trace.ndjson`** total_iters=7000  burnin=2000  sample=5000

```
  [burnin]
    K_before_cull   {'n': 2000, 'min': 8, 'max': 415, 'mean': 115.671, 'median': 100.0, 'stdev': 59.979}
    K_after_cull    {'n': 2000, 'min': 6, 'max': 163, 'mean': 51.934, 'median': 48.0, 'stdev': 21.132}
    spawns_find     {'n': 2000, 'min': 1, 'max': 320, 'mean': 64.203, 'median': 52.0, 'stdev': 43.163}
    spawns_stick    {'n': 2000, 'min': 0, 'max': 87, 'mean': 18.018, 'median': 16.0, 'stdev': 11.193}
    kills_in_cull   {'n': 2000, 'min': 0, 'max': 301, 'mean': 63.736, 'median': 52.0, 'stdev': 43.638}
    cap_hits_find   0
    dp_alpha        {'n': 2000, 'min': 1.1966584092433816, 'max': 49.94011938408166, 'mean': 23.103, 'median': 21.829438738831293, 'stdev': 12.188}
    dp_gamma        {'n': 2000, 'min': 1.0000761474822688, 'max': 1.4761204487394333, 'mean': 1.074, 'median': 1.049525784746184, 'stdev': 0.073}
    alpha_decay     {'n': 2000, 'min': 0.18510337863577198, 'max': 0.7997803146792861, 'mean': 0.558, 'median': 0.5582652275695001, 'stdev': 0.138}
    best_llh        -315.88
  [sample]
    K_before_cull   {'n': 5000, 'min': 11, 'max': 498, 'mean': 113.776, 'median': 104.0, 'stdev': 56.318}
    K_after_cull    {'n': 5000, 'min': 6, 'max': 159, 'mean': 51.724, 'median': 49.0, 'stdev': 20.166}
    spawns_find     {'n': 5000, 'min': 2, 'max': 371, 'mean': 62.594, 'median': 54.0, 'stdev': 40.574}
    spawns_stick    {'n': 5000, 'min': 0, 'max': 80, 'mean': 17.934, 'median': 16.0, 'stdev': 11.118}
    kills_in_cull   {'n': 5000, 'min': 0, 'max': 367, 'mean': 62.051, 'median': 53.0, 'stdev': 41.08}
    cap_hits_find   0
    dp_alpha        {'n': 5000, 'min': 2.3605488533716867, 'max': 49.9878137567823, 'mean': 22.379, 'median': 20.107868977751423, 'stdev': 12.068}
    dp_gamma        {'n': 5000, 'min': 1.0000097355835191, 'max': 1.5432642674611339, 'mean': 1.073, 'median': 1.0503756804099154, 'stdev': 0.074}
    alpha_decay     {'n': 5000, 'min': 0.14093464869397643, 'max': 0.7999527678653265, 'mean': 0.559, 'median': 0.565714299261292, 'stdev': 0.135}
    best_llh        -325.01
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=4-16  median= 13.0  mean=11.95  mode=(13, 279)  best K=5 llh=-315.80
Go local           N= 5000  range=3-15  median= 12.0  mean=11.42  mode=(12, 1663)  best K=5 llh=-325.01
```

### K3_S1_T200_M30_C0_rep0

**trace: `trace.ndjson`** total_iters=7000  burnin=2000  sample=5000

```
  [burnin]
    K_before_cull   {'n': 2000, 'min': 2, 'max': 35, 'mean': 5.821, 'median': 5.0, 'stdev': 3.256}
    K_after_cull    {'n': 2000, 'min': 2, 'max': 8, 'mean': 3.258, 'median': 3.0, 'stdev': 0.649}
    spawns_find     {'n': 2000, 'min': 0, 'max': 33, 'mean': 2.526, 'median': 2.0, 'stdev': 3.121}
    spawns_stick    {'n': 2000, 'min': 0, 'max': 11, 'mean': 0.234, 'median': 0.0, 'stdev': 0.944}
    kills_in_cull   {'n': 2000, 'min': 0, 'max': 33, 'mean': 2.563, 'median': 2.0, 'stdev': 3.142}
    cap_hits_find   0
    dp_alpha        {'n': 2000, 'min': 1.4174458967736467, 'max': 49.97036803782061, 'mean': 31.272, 'median': 32.82074083131809, 'stdev': 12.267}
    dp_gamma        {'n': 2000, 'min': 1.0000024643099588, 'max': 1.519261160350235, 'mean': 1.064, 'median': 1.0441644040062599, 'stdev': 0.064}
    alpha_decay     {'n': 2000, 'min': 0.05000540116593236, 'max': 0.3996965988052771, 'mean': 0.109, 'median': 0.09458281650015782, 'stdev': 0.052}
    best_llh        -107.75
  [sample]
    K_before_cull   {'n': 5000, 'min': 3, 'max': 27, 'mean': 5.943, 'median': 5.0, 'stdev': 3.354}
    K_after_cull    {'n': 5000, 'min': 3, 'max': 10, 'mean': 3.27, 'median': 3.0, 'stdev': 0.636}
    spawns_find     {'n': 5000, 'min': 0, 'max': 24, 'mean': 2.629, 'median': 2.0, 'stdev': 3.224}
    spawns_stick    {'n': 5000, 'min': 0, 'max': 12, 'mean': 0.21, 'median': 0.0, 'stdev': 0.856}
    kills_in_cull   {'n': 5000, 'min': 0, 'max': 24, 'mean': 2.673, 'median': 2.0, 'stdev': 3.247}
    cap_hits_find   0
    dp_alpha        {'n': 5000, 'min': 1.0522464730806382, 'max': 49.99079655324331, 'mean': 30.297, 'median': 31.234800616134258, 'stdev': 12.381}
    dp_gamma        {'n': 5000, 'min': 1.000009362937358, 'max': 1.4659447213445786, 'mean': 1.064, 'median': 1.0442328411753783, 'stdev': 0.063}
    alpha_decay     {'n': 5000, 'min': 0.050016255366911615, 'max': 0.774552417152469, 'mean': 0.113, 'median': 0.09770879839146085, 'stdev': 0.061}
    best_llh        -107.75
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=3- 5  median=  3.0  mean=3.19  mode=(3, 820)  best K=3 llh=-108.04
Go local           N= 5000  range=2- 4  median=  2.0  mean=2.17  mode=(2, 4171)  best K=2 llh=-107.75
```

### K5_S1_T200_M50_C2_rep0

**trace: `trace.ndjson`** total_iters=7000  burnin=2000  sample=5000

```
  [burnin]
    K_before_cull   {'n': 2000, 'min': 3, 'max': 40, 'mean': 8.9, 'median': 8.0, 'stdev': 4.429}
    K_after_cull    {'n': 2000, 'min': 3, 'max': 11, 'mean': 4.487, 'median': 4.0, 'stdev': 0.885}
    spawns_find     {'n': 2000, 'min': 0, 'max': 35, 'mean': 4.367, 'median': 3.0, 'stdev': 4.127}
    spawns_stick    {'n': 2000, 'min': 0, 'max': 11, 'mean': 0.546, 'median': 0.0, 'stdev': 1.37}
    kills_in_cull   {'n': 2000, 'min': 0, 'max': 36, 'mean': 4.413, 'median': 3.0, 'stdev': 4.185}
    cap_hits_find   0
    dp_alpha        {'n': 2000, 'min': 1.0249761781091744, 'max': 49.90784769400419, 'mean': 23.27, 'median': 21.908939924625685, 'stdev': 12.998}
    dp_gamma        {'n': 2000, 'min': 1.0000393437057284, 'max': 1.5015939589610345, 'mean': 1.067, 'median': 1.0477389929377008, 'stdev': 0.066}
    alpha_decay     {'n': 2000, 'min': 0.050047430947452315, 'max': 0.7031882353646978, 'mean': 0.128, 'median': 0.10132355080717292, 'stdev': 0.087}
    best_llh        -6648.68
  [sample]
    K_before_cull   {'n': 5000, 'min': 4, 'max': 51, 'mean': 8.587, 'median': 7.0, 'stdev': 4.414}
    K_after_cull    {'n': 5000, 'min': 4, 'max': 13, 'mean': 4.445, 'median': 4.0, 'stdev': 0.832}
    spawns_find     {'n': 5000, 'min': 0, 'max': 44, 'mean': 4.111, 'median': 3.0, 'stdev': 4.14}
    spawns_stick    {'n': 5000, 'min': 0, 'max': 16, 'mean': 0.513, 'median': 0.0, 'stdev': 1.325}
    kills_in_cull   {'n': 5000, 'min': 0, 'max': 44, 'mean': 4.142, 'median': 3.0, 'stdev': 4.169}
    cap_hits_find   0
    dp_alpha        {'n': 5000, 'min': 1.070835935940536, 'max': 49.97393825632721, 'mean': 24.297, 'median': 23.24794030801978, 'stdev': 12.907}
    dp_gamma        {'n': 5000, 'min': 1.0000027180152962, 'max': 1.600380280011489, 'mean': 1.069, 'median': 1.0485681383424956, 'stdev': 0.068}
    alpha_decay     {'n': 5000, 'min': 0.05002348595548358, 'max': 0.7985012282652016, 'mean': 0.114, 'median': 0.08787690679820226, 'stdev': 0.081}
    best_llh        -6648.57
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=4- 6  median=  4.0  mean=4.13  mode=(4, 880)  best K=4 llh=-6649.55
Go local           N= 5000  range=3- 5  median=  3.0  mean=3.14  mode=(3, 4371)  best K=3 llh=-6648.57
```

### K5_S1_T200_M50_C0_rep0

**trace: `trace.ndjson`** total_iters=7000  burnin=2000  sample=5000

```
  [burnin]
    K_before_cull   {'n': 2000, 'min': 4, 'max': 223, 'mean': 28.968, 'median': 26.0, 'stdev': 16.42}
    K_after_cull    {'n': 2000, 'min': 4, 'max': 48, 'mean': 15.026, 'median': 14.0, 'stdev': 6.854}
    spawns_find     {'n': 2000, 'min': 0, 'max': 183, 'mean': 14.372, 'median': 12.0, 'stdev': 12.243}
    spawns_stick    {'n': 2000, 'min': 0, 'max': 31, 'mean': 3.944, 'median': 3.0, 'stdev': 4.073}
    kills_in_cull   {'n': 2000, 'min': 0, 'max': 178, 'mean': 13.941, 'median': 11.0, 'stdev': 12.214}
    cap_hits_find   0
    dp_alpha        {'n': 2000, 'min': 1.023810235390486, 'max': 49.997178999782676, 'mean': 10.585, 'median': 6.466633932591202, 'stdev': 10.201}
    dp_gamma        {'n': 2000, 'min': 1.0000261394063104, 'max': 1.5871675181337284, 'mean': 1.08, 'median': 1.0556987361616634, 'stdev': 0.08}
    alpha_decay     {'n': 2000, 'min': 0.05064135393972493, 'max': 0.7999538974658993, 'mean': 0.398, 'median': 0.38508778017463696, 'stdev': 0.195}
    best_llh        -214.03
  [sample]
    K_before_cull   {'n': 5000, 'min': 4, 'max': 206, 'mean': 31.684, 'median': 28.0, 'stdev': 18.016}
    K_after_cull    {'n': 5000, 'min': 4, 'max': 54, 'mean': 15.996, 'median': 15.0, 'stdev': 7.102}
    spawns_find     {'n': 5000, 'min': 0, 'max': 168, 'mean': 16.051, 'median': 13.0, 'stdev': 13.6}
    spawns_stick    {'n': 5000, 'min': 0, 'max': 39, 'mean': 4.483, 'median': 3.0, 'stdev': 4.583}
    kills_in_cull   {'n': 5000, 'min': 0, 'max': 160, 'mean': 15.688, 'median': 12.0, 'stdev': 13.622}
    cap_hits_find   0
    dp_alpha        {'n': 5000, 'min': 1.0046354230814627, 'max': 49.69527850515377, 'mean': 8.788, 'median': 5.902059185851227, 'stdev': 8.099}
    dp_gamma        {'n': 5000, 'min': 1.0000331654496273, 'max': 1.8019100665654446, 'mean': 1.077, 'median': 1.0539315908478941, 'stdev': 0.075}
    alpha_decay     {'n': 5000, 'min': 0.050372754032151065, 'max': 0.7995596351469814, 'mean': 0.438, 'median': 0.42666600161092183, 'stdev': 0.185}
    best_llh        -212.95
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=4-11  median=  6.0  mean=6.20  mode=(7, 273)  best K=4 llh=-212.87
Go local           N= 5000  range=3-10  median=  6.0  mean=6.06  mode=(6, 1666)  best K=3 llh=-212.95
```

## expC_multichain — Multi-chain variance K3_C2 (B=500, S=1000, 4 chains)

### K3_S1_T200_M30_C2_rep0

**trace: `trace.ndjson.chain0`** total_iters=1500  burnin=500  sample=1000

```
  [burnin]
    K_before_cull   {'n': 500, 'min': 7, 'max': 301, 'mean': 85.908, 'median': 76.0, 'stdev': 49.347}
    K_after_cull    {'n': 500, 'min': 5, 'max': 137, 'mean': 41.414, 'median': 40.0, 'stdev': 19.888}
    spawns_find     {'n': 500, 'min': 1, 'max': 202, 'mean': 45.036, 'median': 38.0, 'stdev': 33.079}
    spawns_stick    {'n': 500, 'min': 0, 'max': 63, 'mean': 13.91, 'median': 12.0, 'stdev': 10.041}
    kills_in_cull   {'n': 500, 'min': 0, 'max': 210, 'mean': 44.494, 'median': 36.0, 'stdev': 33.522}
    cap_hits_find   0
    dp_alpha        {'n': 500, 'min': 1.655247782301983, 'max': 49.88073076409877, 'mean': 27.087, 'median': 26.243306373815216, 'stdev': 13.513}
    dp_gamma        {'n': 500, 'min': 1.000079046560237, 'max': 1.3280133522754998, 'mean': 1.075, 'median': 1.0565362195462469, 'stdev': 0.067}
    alpha_decay     {'n': 500, 'min': 0.16896162105415524, 'max': 0.7997766718562576, 'mean': 0.467, 'median': 0.4301138890438815, 'stdev': 0.161}
    best_llh        -324.91
  [sample]
    K_before_cull   {'n': 1000, 'min': 14, 'max': 413, 'mean': 110.383, 'median': 103.0, 'stdev': 48.939}
    K_after_cull    {'n': 1000, 'min': 8, 'max': 114, 'mean': 50.991, 'median': 49.0, 'stdev': 17.968}
    spawns_find     {'n': 1000, 'min': 1, 'max': 311, 'mean': 59.621, 'median': 52.0, 'stdev': 35.758}
    spawns_stick    {'n': 1000, 'min': 0, 'max': 68, 'mean': 18.053, 'median': 17.0, 'stdev': 10.19}
    kills_in_cull   {'n': 1000, 'min': 4, 'max': 312, 'mean': 59.392, 'median': 52.5, 'stdev': 35.913}
    cap_hits_find   0
    dp_alpha        {'n': 1000, 'min': 2.0828042883119315, 'max': 49.93824757182972, 'mean': 23.627, 'median': 21.53262705566292, 'stdev': 12.349}
    dp_gamma        {'n': 1000, 'min': 1.0000667998631967, 'max': 1.4074816738207523, 'mean': 1.072, 'median': 1.052498265023619, 'stdev': 0.068}
    alpha_decay     {'n': 1000, 'min': 0.23033396730914493, 'max': 0.797487906860112, 'mean': 0.546, 'median': 0.5704948425226439, 'stdev': 0.133}
    best_llh        -334.85
```

**trace: `trace.ndjson.chain1`** total_iters=1500  burnin=500  sample=1000

```
  [burnin]
    K_before_cull   {'n': 500, 'min': 24, 'max': 418, 'mean': 127.028, 'median': 116.0, 'stdev': 54.938}
    K_after_cull    {'n': 500, 'min': 16, 'max': 141, 'mean': 56.34, 'median': 54.0, 'stdev': 19.385}
    spawns_find     {'n': 500, 'min': 11, 'max': 307, 'mean': 71.496, 'median': 62.0, 'stdev': 40.997}
    spawns_stick    {'n': 500, 'min': 0, 'max': 77, 'mean': 19.662, 'median': 18.0, 'stdev': 11.188}
    kills_in_cull   {'n': 500, 'min': 5, 'max': 328, 'mean': 70.688, 'median': 62.0, 'stdev': 41.468}
    cap_hits_find   0
    dp_alpha        {'n': 500, 'min': 3.1743092623851306, 'max': 49.97168353252337, 'mean': 20.336, 'median': 17.971781878154644, 'stdev': 10.9}
    dp_gamma        {'n': 500, 'min': 1.000193089689036, 'max': 1.5362187337423774, 'mean': 1.079, 'median': 1.0597758749168857, 'stdev': 0.077}
    alpha_decay     {'n': 500, 'min': 0.31658741050753225, 'max': 0.799705327389071, 'mean': 0.6, 'median': 0.618793088710232, 'stdev': 0.126}
    best_llh        -354.26
  [sample]
    K_before_cull   {'n': 1000, 'min': 13, 'max': 248, 'mean': 82.928, 'median': 78.0, 'stdev': 33.48}
    K_after_cull    {'n': 1000, 'min': 6, 'max': 96, 'mean': 41.07, 'median': 40.0, 'stdev': 14.718}
    spawns_find     {'n': 1000, 'min': 1, 'max': 155, 'mean': 42.103, 'median': 38.0, 'stdev': 23.226}
    spawns_stick    {'n': 1000, 'min': 0, 'max': 52, 'mean': 13.638, 'median': 12.0, 'stdev': 8.612}
    kills_in_cull   {'n': 1000, 'min': 1, 'max': 180, 'mean': 41.858, 'median': 38.0, 'stdev': 23.548}
    cap_hits_find   0
    dp_alpha        {'n': 1000, 'min': 3.6689780642550662, 'max': 49.99832161280508, 'mean': 27.927, 'median': 27.975391749426727, 'stdev': 12.457}
    dp_gamma        {'n': 1000, 'min': 1.0000202370381017, 'max': 1.4691112658366121, 'mean': 1.068, 'median': 1.0487993367379553, 'stdev': 0.063}
    alpha_decay     {'n': 1000, 'min': 0.20290473475961773, 'max': 0.7982964812822851, 'mean': 0.455, 'median': 0.44807123321778863, 'stdev': 0.108}
    best_llh        -330.56
```

**trace: `trace.ndjson.chain2`** total_iters=1500  burnin=500  sample=1000

```
  [burnin]
    K_before_cull   {'n': 500, 'min': 19, 'max': 324, 'mean': 123.108, 'median': 111.0, 'stdev': 54.027}
    K_after_cull    {'n': 500, 'min': 11, 'max': 132, 'mean': 55.588, 'median': 54.0, 'stdev': 19.887}
    spawns_find     {'n': 500, 'min': 5, 'max': 233, 'mean': 68.402, 'median': 58.0, 'stdev': 39.535}
    spawns_stick    {'n': 500, 'min': 0, 'max': 70, 'mean': 19.646, 'median': 18.0, 'stdev': 11.465}
    kills_in_cull   {'n': 500, 'min': 1, 'max': 217, 'mean': 67.52, 'median': 59.0, 'stdev': 39.705}
    cap_hits_find   0
    dp_alpha        {'n': 500, 'min': 2.580611948061443, 'max': 49.74355219684444, 'mean': 22.969, 'median': 20.900394798082438, 'stdev': 10.959}
    dp_gamma        {'n': 500, 'min': 1.0000642771065766, 'max': 1.5086480919709142, 'mean': 1.077, 'median': 1.0507912326170759, 'stdev': 0.077}
    alpha_decay     {'n': 500, 'min': 0.1954395485669789, 'max': 0.7993055699175549, 'mean': 0.566, 'median': 0.5680319059106727, 'stdev': 0.126}
    best_llh        -348.47
  [sample]
    K_before_cull   {'n': 1000, 'min': 26, 'max': 535, 'mean': 116.714, 'median': 102.0, 'stdev': 59.089}
    K_after_cull    {'n': 1000, 'min': 9, 'max': 156, 'mean': 52.988, 'median': 49.0, 'stdev': 20.89}
    spawns_find     {'n': 1000, 'min': 7, 'max': 370, 'mean': 64.376, 'median': 54.0, 'stdev': 42.461}
    spawns_stick    {'n': 1000, 'min': 0, 'max': 76, 'mean': 18.138, 'median': 16.0, 'stdev': 11.428}
    kills_in_cull   {'n': 1000, 'min': 4, 'max': 422, 'mean': 63.726, 'median': 53.0, 'stdev': 43.335}
    cap_hits_find   0
    dp_alpha        {'n': 1000, 'min': 3.7885070794737117, 'max': 49.62377511542914, 'mean': 24.433, 'median': 23.56349172542359, 'stdev': 11.704}
    dp_gamma        {'n': 1000, 'min': 1.0000544452971845, 'max': 1.6034555886293202, 'mean': 1.077, 'median': 1.0551388267681483, 'stdev': 0.074}
    alpha_decay     {'n': 1000, 'min': 0.22334747743930494, 'max': 0.7990730424344877, 'mean': 0.551, 'median': 0.5587442389210964, 'stdev': 0.129}
    best_llh        -330.41
```

**trace: `trace.ndjson.chain3`** total_iters=1500  burnin=500  sample=1000

```
  [burnin]
    K_before_cull   {'n': 500, 'min': 6, 'max': 379, 'mean': 120.044, 'median': 109.0, 'stdev': 62.907}
    K_after_cull    {'n': 500, 'min': 5, 'max': 132, 'mean': 53.974, 'median': 52.0, 'stdev': 22.528}
    spawns_find     {'n': 500, 'min': 0, 'max': 291, 'mean': 66.826, 'median': 55.0, 'stdev': 45.036}
    spawns_stick    {'n': 500, 'min': 0, 'max': 63, 'mean': 18.69, 'median': 17.0, 'stdev': 12.133}
    kills_in_cull   {'n': 500, 'min': 0, 'max': 247, 'mean': 66.07, 'median': 57.5, 'stdev': 45.719}
    cap_hits_find   0
    dp_alpha        {'n': 500, 'min': 5.216731927292817, 'max': 49.57446374163151, 'mean': 25.593, 'median': 24.4470925919399, 'stdev': 11.354}
    dp_gamma        {'n': 500, 'min': 1.000030151779088, 'max': 1.3779576454794211, 'mean': 1.067, 'median': 1.048805244750425, 'stdev': 0.064}
    alpha_decay     {'n': 500, 'min': 0.16286633074747883, 'max': 0.7998476975809081, 'mean': 0.552, 'median': 0.548433553637925, 'stdev': 0.148}
    best_llh        -325.84
  [sample]
    K_before_cull   {'n': 1000, 'min': 21, 'max': 362, 'mean': 116.547, 'median': 108.0, 'stdev': 49.453}
    K_after_cull    {'n': 1000, 'min': 7, 'max': 221, 'mean': 53.111, 'median': 51.0, 'stdev': 19.075}
    spawns_find     {'n': 1000, 'min': 3, 'max': 249, 'mean': 64.369, 'median': 55.0, 'stdev': 36.265}
    spawns_stick    {'n': 1000, 'min': 0, 'max': 79, 'mean': 18.205, 'median': 16.0, 'stdev': 10.559}
    kills_in_cull   {'n': 1000, 'min': 4, 'max': 272, 'mean': 63.436, 'median': 56.0, 'stdev': 36.254}
    cap_hits_find   9
    dp_alpha        {'n': 1000, 'min': 3.6875114468259063, 'max': 49.94454047867749, 'mean': 21.249, 'median': 18.363141019329504, 'stdev': 11.957}
    dp_gamma        {'n': 1000, 'min': 1.0000493966055943, 'max': 1.9225411063866076, 'mean': 1.068, 'median': 1.0455776902610903, 'stdev': 0.076}
    alpha_decay     {'n': 1000, 'min': 0.2705047662753809, 'max': 0.7994039511138272, 'mean': 0.583, 'median': 0.5989177235981407, 'stdev': 0.116}
    best_llh        -323.20
```

**posterior K distribution (vs Python reference)**

```
Python (orig)      N= 1000  range=4-16  median= 13.0  mean=11.95  mode=(13, 279)  best K=5 llh=-315.80
Go local           N= 4000  range=3-16  median= 12.0  mean=11.40  mode=(12, 1325)  best K=4 llh=-323.20
```
