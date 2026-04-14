# Cloud Energetic Classification — ClimateManifoldDynamics/cloud

**Author:** László Tatai  
**ORCID:** 0009-0007-5153-6306  
**Version:** v3.0 (April 2026)  
**Part of:** [ClimateManifoldDynamics](https://github.com/LesliePi/ClimateManifoldDynamics)  
**Theory DOI:** https://doi.org/10.5281/zenodo.19568175  
**License:** Apache 2.0 WITH Commons Clause v1.0

---

## What this is

This module introduces a cloud-resolved energetic descriptor **R** — the ratio of dynamic
to radiative energy flux — and characterises how ERA5 atmospheric observations are
constrained to a curved surface in (T, C, R) state space.

It is a Supplementary Note to the Thermodynamic Manifold Framework (TMF), providing a
measurable bridge between observable cloud fields and the effective dissipation κ(μ).

---

## Key results (ERA5, 1980–2025, n = 16.3M observations)

- **Physical envelope:** R_max(T,C) = (0.874·C + 0.017)·exp(−0.004·T) / σ·T_K⁴   (R² = 0.937)
- **Asymmetry at ρ = 1:** AI = +0.111 (slow build-up, fast collapse), n = 700,269 events
- **Spontaneous organisation threshold:** ρ ≈ 0.75 (drift field zero-crossing)
- **State estimator:** P(regime_{t+k} | ρ_t, T_t, C_t) from historical ERA5 statistics

---

## Files

```
cloud/
├── cloud_manifold_v3_0.py     # Main pipeline (ERA5 load → fit → classify → dynamics → plots)
├── Supplementary_Note_S1.md   # Full theoretical and empirical description
└── README.md                  # This file
```

---

## Requirements

```
python >= 3.10
numpy, pandas, matplotlib, xarray, scipy
```

Install:
```bash
pip install numpy pandas matplotlib xarray scipy pyarrow
```

---

## Usage

```bash
# Edit DATA_PATH in the config section to point to your ERA5 .nc files
python cloud_manifold_v3_0.py
```

Output directory: `output/`  
Parquet cache: `output/era5_cloud_processed.parquet` (auto-generated on first run)

### ERA5 variables needed

Download from the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/):

- `2m_temperature` (t2m)
- `total_cloud_cover` (tcc)
- `2m_dewpoint_temperature` (d2m)  ← preferred for RH
- `surface_pressure` (sp)  ← optional, used with specific humidity

---

## Pipeline sections

| Section | Description |
|---------|-------------|
| 1 | Config and constants |
| 2 | Utilities |
| 3 | Physics (R, saturation, RH) |
| 4 | ERA5 load + Parquet cache |
| 5 | Surface fit: parametric + GP, ρ = R/R_max |
| 6 | ρ-based regime classification |
| 7 | P(regime\|T,C) probabilistic module + entropy map |
| 8 | Asymmetry analysis at ρ = 1 |
| 9 | Conditional drift field E[dρ/dt \| ρ, T] |
| 10 | Probabilistic state estimator P(regime\|ρ,T,C,k) |
| 11 | Main |

---

## State estimator

```python
import cloud_manifold_v3_0 as mod

# After fitting (see main()):
probs = se.query_state(rho=0.92, T=-5.0, C=0.9, k=1)
# → {'clear': 0.10, 'convective': 0.18, 'mixed': 0.63, ...}
```

**Note:** This is a probability estimate from historical statistics — not a weather forecast.

---

## Contributing

This repository is maintained as an open library for cloud thermodynamic analysis
within the TMF framework.

**Contributions are welcome.** If you develop an extension — regional model, pressure-level
version, validation against CERES/MODIS, etc. — please submit a pull request or contact
the maintainer. **Author attribution is strictly maintained for all contributions.**

Guidelines:
- One subdirectory per contribution
- Include author name and ORCID in file headers
- Include a short README in your subdirectory
- Do not modify existing files — extend them or add new ones

---

## Citation

```bibtex
@software{tatai_cloud_manifold_2026,
  author  = {Tatai, László},
  title   = {Cloud Energetic Classification within the TMF Framework},
  year    = {2026},
  version = {v3.0},
  url     = {https://github.com/LesliePi/ClimateManifoldDynamics/tree/main/cloud},
  note    = {Part of ClimateManifoldDynamics. Related theory: doi:10.5281/zenodo.19430594}
}
```
