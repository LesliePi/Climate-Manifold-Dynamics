# Zenodo Upload Metadata
# Cloud Energetic Classification and Edge Dynamics within the TMF Framework

## Title
Cloud Energetic Classification and Edge Dynamics within the  
Thermodynamic Manifold Framework — Supplementary Note S1 with Software (v3.0)

## Authors
- László Tatai (ORCID: 0009-0007-5153-6306)

## Description
This deposit contains the Supplementary Note S1 and accompanying open-source software
(cloud_manifold_v3_0.py) for the cloud energetic classification module of the
Thermodynamic Manifold Framework (TMF).

The note introduces a dimensionless cloud energy ratio R = F_dyn / F_rad, shows that
ERA5 atmospheric observations are constrained to a curved surface R_max(T,C) in
(T, C, R) state space, and characterises the edge dynamics at the physical boundary
ρ = R / R_max = 1.

Key empirical results from ERA5 (1980–2025, n = 16.3M):
  - Parametric envelope fit R² = 0.937
  - Asymmetry index AI = +0.111 (slow build-up, fast collapse), n = 700,269 events
  - Spontaneous organisation threshold at ρ ≈ 0.75
  - Probabilistic state estimator P(regime | ρ, T, C, k)

Related TMF paper: https://doi.org/10.5281/zenodo.19430594

## Keywords
cloud classification, thermodynamic manifold, ERA5, cloud energy ratio,
convective dynamics, probabilistic state estimation, climate physics,
radiative-convective balance, atmospheric state space

## License
Apache License 2.0 WITH Commons Clause v1.0

## Resource type
Software + Supplementary Note

## Related identifiers
- Is supplement to: https://doi.org/10.5281/zenodo.19430594 (TMF v2.0)
- Is part of: https://github.com/LesliePi/ClimateManifoldDynamics

## Files to upload
1. Supplementary_Note_S1.md        ← theoretical and empirical description
2. cloud_manifold_v3_0.py          ← full pipeline software
3. README.md                       ← usage instructions
4. edge_asymmetry.png              ← Figure 1: asymmetry at ρ = 1
5. drift_field.png                 ← Figure 2: conditional drift field
6. state_estimator_examples.png    ← Figure 3: probabilistic state estimator
7. R_3D_regimes.png                ← Figure 4: cloud regimes in R-T-Cloud space

## Version notes
v3.0 — adds edge dynamics (Section 8), drift field (Section 9),
        probabilistic state estimator (Section 10) to the v2.0 envelope analysis.
