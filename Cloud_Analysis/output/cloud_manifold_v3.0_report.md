# Cloud Manifold Pipeline — Run Report

**Version:** v3.0  
**Generated:** 2026-04-14 12:26:39  

---

## Surface Fit

| Model | RMSE | MAE | R² |
|-------|------|-----|----|
| parametric | 2.657e-04 | 1.914e-04 | 0.9370 |
| gp | 4.113e-04 | 3.048e-04 | 0.8491 |

**Parametric model:**  R_max(T,C) = (0.8744·C + 0.0166)·exp(−0.0044·T) / σ·T_K⁴

---

## Probabilistic Module

Grid: 40 × 20  |  Valid cells: 639  |  Mean H: 0.329 bits

---

## Edge Asymmetry

| Metric | Value |
|--------|-------|
| Events detected | 700,269 |
| Median rise time | 4.00 steps |
| Median fall time | 5.00 steps |
| Asymmetry index AI | +0.1111 |
| Interpretation | slow build / fast collapse |

---

## Full Log

```
=================================================================
CLOUD MANIFOLD PIPELINE  v3.0
Run started: 2026-04-14 12:19:18
=================================================================
Loading cached Parquet: output\era5_cloud_processed.parquet
  Loaded 16,346,138 rows | years 1980–2025

[Section 5] Surface fit module

Envelope grid: 639 cells with ≥ 30 observations

Parametric fit:  a = 0.87437   b = 0.01661   alpha = 0.00444
  R_max(0°C, C=0.5)  = 0.001438
  R_max(20°C, C=0.5) = 0.000992
  RMSE=2.66e-04   MAE=1.91e-04   R²=0.9370

GP fit:  l=10.000   sigma=0.010   noise=1.0e-04
  RMSE=4.11e-04   MAE=3.05e-04   R²=0.8491

── Envelope fit comparison ──────────────────────────────
Model                 RMSE         MAE        R²
──────────────  ──────────  ──────────  ────────
parametric       2.657e-04   1.914e-04    0.9370
gp               4.113e-04   3.048e-04    0.8491
──────────────────────────────────────────────────

[Section 6] Regime classification (ρ-based)

Cloud regime distribution:
  clear         :  5,345,432  (32.7 %)
  convective    :  2,418,221  (14.8 %)
  stratiform    :          0  (0.0 %)
  mixed         :  6,348,707  (38.8 %)
  transitional  :  2,233,778  (13.7 %)

[Section 7] Probabilistic module

Probabilistic module:
  Grid: 40 × 20 cells | valid cells: 639
  Mean entropy (valid cells): 0.329 bits
  Max  entropy (valid cells): 1.000 bits

── Probabilistic module summary ─────────────────────────
  clear         : mean=0.205  max=1.000  cells_dominant=131
  convective    : mean=0.090  max=0.584  cells_dominant=2
  stratiform    : mean=0.000  max=0.000  cells_dominant=0
  mixed         : mean=0.436  max=1.000  cells_dominant=315
  transitional  : mean=0.269  max=1.000  cells_dominant=191
──────────────────────────────────────────────────────

Example query  P(regime | T=15.0°C, C=0.7):
  clear         : 0.000  
  convective    : 0.128  ███
  stratiform    : 0.000  
  mixed         : 0.872  █████████████████
  transitional  : 0.000  

[Section 8–10 prep] Building ρ time series ...

[Section 8] Asymmetry analysis at ρ = 1

Edge event analysis:
  Complete cycles detected : 700,269
  Median rise time  τ_rise : 4.00 steps
  Median fall time  τ_fall : 5.00 steps
  Asymmetry index   AI     : +0.1111  (slow build / fast collapse)
  Saved: output\20260414_122324_edge_asymmetry.png

[Section 9] Conditional drift field

Drift field:
  Grid: 20 × 30 | valid cells: 434
  Fraction of cells with positive drift (push toward ρ=1): 0.255
  Drift zero-crossing (equilibrium): ρ ≈ 0.750
  Saved: output\20260414_122359_drift_field.png

[Section 10] Probabilistic state estimator

State estimator:
  Transition kernel: 15 × 20 × 15 | valid (ρ, T) cells: 213
  Regime transition matrix built from 16,315,710 observed transitions

Regime transition matrix  P(regime_{t+1} | regime_t):
                   clear  convec  strati   mixed  transi
  ──────────────────────────────────────────────────────
  clear            0.483   0.090   0.000   0.290   0.136
  convective       0.198   0.325   0.000   0.383   0.094
  stratiform       0.200   0.200   0.200   0.200   0.200
  mixed            0.246   0.143   0.000   0.467   0.143
  transitional     0.324   0.105   0.000   0.406   0.165

State estimator queries:

  State: low ρ, mild  (ρ=0.3, T=10.0°C, C=0.3)
    k=1: dominant=clear          clea=0.91  conv=0.00  stra=0.00  mixe=0.06  tran=0.03
    k=3: dominant=clear          clea=0.94  conv=0.03  stra=0.00  mixe=0.02  tran=0.01
    k=6: dominant=clear          clea=0.95  conv=0.04  stra=0.00  mixe=0.01  tran=0.00

  State: mid ρ, cool  (ρ=0.7, T=5.0°C, C=0.7)
    k=1: dominant=clear          clea=0.58  conv=0.10  stra=0.00  mixe=0.22  tran=0.09
    k=3: dominant=clear          clea=0.73  conv=0.18  stra=0.00  mixe=0.07  tran=0.02
    k=6: dominant=clear          clea=0.76  conv=0.22  stra=0.00  mixe=0.02  tran=0.00

  State: near-envelope, cold  (ρ=0.92, T=-5.0°C, C=0.9)
    k=1: dominant=mixed          clea=0.10  conv=0.18  stra=0.00  mixe=0.63  tran=0.10
    k=3: dominant=clear          clea=0.51  conv=0.13  stra=0.00  mixe=0.30  tran=0.06
    k=6: dominant=clear          clea=0.87  conv=0.04  stra=0.00  mixe=0.07  tran=0.01

  State: at edge, warm  (ρ=0.98, T=15.0°C, C=0.85)
    k=1: dominant=convective     clea=0.03  conv=0.56  stra=0.00  mixe=0.32  tran=0.09
    k=3: dominant=convective     clea=0.03  conv=0.84  stra=0.00  mixe=0.10  tran=0.03
    k=6: dominant=convective     clea=0.05  conv=0.91  stra=0.00  mixe=0.03  tran=0.01

[Section 11] Generating all plots ...
  Saved: output\20260414_122540_surface_comparison.png
  Saved: output\20260414_122540_surface_residuals.png
  Saved: output\20260414_122546_rho_scatter.png
  Saved: output\20260414_122549_R_3D_regimes.png
  Saved: output\20260414_122552_regime_probability_maps.png
  Saved: output\20260414_122552_entropy_map.png
  Saved: output\20260414_122630_rho_distribution_by_regime.png
  Saved: output\20260414_122633_transition_kernel.png
  Saved: output\20260414_122635_regime_transition_matrix.png
  Saved: output\20260414_122637_state_estimator_examples.png
```
