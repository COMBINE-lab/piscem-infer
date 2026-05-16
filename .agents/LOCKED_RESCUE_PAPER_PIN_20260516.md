# Locked Rescue Paper Pin, 2026-05-16

This note records the implementation and parameters used for the current paper
method, plus the near-identical EB-min comparison point.

## Code State

Safety tag to create after commit:

```text
paper-locked-rescue-floor-smooth-ebmin-20260516
```

The intended default paper implementation is condition-rescue locked-mass phase
2 with structural transcript selection and residual EC collapsing after
subtracting locked allocations. Condition-rescued transcripts are included in
phase 2 through locked mass, and their final estimates are not overwritten with
raw phase-1 estimates when locked rescue is enabled.

## Paper Method

Default command shape:

```bash
target/release/piscem-infer consensus-quant \
  -m <manifest.csv> \
  -l <library_type> \
  --filter-mode support \
  --txp-selection \
  --lock-condition-rescue-allocations \
  --phase2-max-iter 1500 \
  --phase2-convergence-thresh 0.0005
```

The `--condition-rescue-lock-mode` default is
`floor-smooth-confidence` when `--lock-condition-rescue-allocations` is set.

Effective locked-rescue parameters:

```text
condition_rescue_lock_mode = floor-smooth-confidence
condition_rescue_min_lock_threshold = 0.1
condition_rescue_full_lock_threshold = 0.5
condition_rescue_lock_fraction = 0.75
phase2_max_iter = 1500
phase2_convergence_thresh = 0.0005
filter_mode = support
txp_selection = true
```

Library types used in benchmark reruns:

```text
simulation: IU
SEQC:       IU
airway:     ISR
```

## EB-min Comparison

EB-min is the same implementation and command as the paper method, except:

```text
condition_rescue_min_lock_threshold = 0.0584
```

Command addition:

```bash
--condition-rescue-min-lock-threshold 0.0584
```

This value came from the empirical-Bayes audit calculation where the inferred
signal posterior crossed the conservative operating point near rescued
posterior fraction 0.0584. It was close to, but not consistently better than,
the 0.1 paper default.

## Benchmark Outputs Used For Comparison

Simulation outputs:

```text
sim_data_gencode/rerun_pims_rescue_floor_smooth_current_default_20260515
sim_data_gencode/rerun_pims_rescue_floor_smooth_ebmin058_txpsel_20260515
```

External EB-min outputs:

```text
airway_benchmark/quant/putative_locked_rescue_floor_smooth_ebmin058_20260516
seqc_benchmark/quant_putative_locked_rescue_floor_smooth_ebmin058_20260516
```

External paper-method outputs are the corresponding floor-smooth confidence
locked-rescue runs with the 0.1 minimum threshold and structural transcript
selection enabled.

## Count-Space Simulation Summary

Exact-duplicate-collapsed count-space metrics:

```text
paper floor-smooth 0.1:
  TP=4809.5 FP=176.2 FN=189.5
  precision=0.964665 recall=0.962092 F1=0.963377
  logRMSE=0.276072 Pearson=0.979929

EB-min 0.0584:
  TP=4811.7 FP=183.8 FN=187.3
  precision=0.963200 recall=0.962526 F1=0.962863
  logRMSE=0.279231 Pearson=0.979477
```

## External Benchmark Summary

Airway untreated replicate CV:

```text
paper floor-smooth 0.1:
  union median CV=0.1823107402208193
  intersection median CV=0.2739505621502545

EB-min 0.0584:
  union median CV=0.18652006060522353
  intersection median CV=0.2739574774617739
```

SEQC titration/CV/qPCR:

```text
paper floor-smooth 0.1:
  Pe(C)=0.9527664317373598
  Pe(D)=0.9518787974468762
  Sp(C)=0.9492675320506603
  Sp(D)=0.9482890730323267
  FC slope=0.9408802679940751
  Det(A)=52567
  union median CV=0.07565492328544615
  intersection median CV=0.2099004883980233
  qPCR Pearson(A)=0.5678610270421885
  qPCR Spearman(A)=0.6927673051367089
  qPCR FC Pearson=0.9056935697887065

EB-min 0.0584:
  Pe(C)=0.9533918223103796
  Pe(D)=0.9521887141504309
  Sp(C)=0.9499460458216682
  Sp(D)=0.9486172768994311
  FC slope=0.9406504823742544
  Det(A)=53110
  union median CV=0.08146107799335907
  intersection median CV=0.20940192077033384
  qPCR Pearson(A)=0.5678019222369249
  qPCR Spearman(A)=0.6928208308032315
  qPCR FC Pearson=0.9060498260982729
```
