# Supplementary Note S1

## Cloud Energetic Classification and Edge Dynamics  
## within the Thermodynamic Manifold Framework

**Author:** László Tatai  
**ORCID:** 0009-0007-5153-6306  
**Date:** April 2026  
**Repository:** https://github.com/LesliePi/ClimateManifoldDynamics/tree/main/cloud  
**Related work:** *Thermodynamic Manifold Dynamics of the Climate System*  
**DOI (TMF v2.0):** https://doi.org/10.5281/zenodo.19430594  

---

## Abstract

This note introduces a cloud-resolved energetic descriptor $R$ that enables a direct,
measurable bridge between observable atmospheric fields and the effective dissipation
structure $\kappa(\mu)$ of the Thermodynamic Manifold Framework (TMF). We define $R$ as
the ratio of dynamic to radiative energy flux, show that ERA5 observations are constrained
to a well-defined curved surface $R_{\max}(T, C)$ in the $(T, C, R)$ state space, and
characterise the dynamics at the physical boundary $\rho = R / R_{\max} = 1$. Key results
include: (i) an analytic parametric fit to the envelope with $R^2 = 0.937$ on 16.3 million
ERA5 observations spanning 1980–2025; (ii) confirmed asymmetry at the edge
(AI $= +0.111$, $n = 700{,}269$ events), consistent with slow thermodynamic build-up and
fast convective collapse; (iii) a conditional drift field revealing a spontaneous
organisation threshold at $\rho \approx 0.75$; and (iv) a probabilistic state estimator
providing $P(\text{regime}_{t+k} \mid \rho_t, T_t, C_t)$ from historical ERA5 statistics.
This framework does not modify the core TMF formulation but substantially extends its
observational testability and practical utility.

---

## S1. Motivation

The Thermodynamic Manifold Framework (TMF) defines climate dynamics through the coupled
evolution of energy and the hydrological state-space measure $\mu(W, t)$. While this
framework captures the global structure of the system, it does not explicitly resolve the
**mesoscale thermodynamic structures** responsible for energy redistribution.

Clouds represent precisely such structures: localised, dynamically evolving regions where
latent heat exchange, radiative transfer, convective transport, and turbulent mixing interact
simultaneously. Their collective effect on $\kappa(\mu)$ — the effective dissipation term of
the TMF — is real but unresolved in the core formulation.

This note introduces a cloud-resolved energetic descriptor that provides:

- a physically interpretable local order parameter derived from observable ERA5 fields
- an analytic description of the admissible thermodynamic state space
- a characterisation of the edge dynamics at the physical boundary
- a probabilistic state estimator linking manifold position to regime probability

---

## S2. Definition of the Cloud Energy Ratio

We define the dimensionless cloud energy ratio:

$$R(x,t) = \frac{F_{\mathrm{dyn}}(x,t)}{F_{\mathrm{rad}}(x,t)}$$

where $F_{\mathrm{dyn}}$ is the dynamic (non-radiative) energy flux and $F_{\mathrm{rad}}$
is the radiative energy flux. The ratio $R$ acts as a **local order parameter** for cloud
thermodynamic behaviour.

### S2.1 Radiative component

$$F_{\mathrm{rad}}(x,t) = \sigma\, T(x,t)^4$$

where $\sigma = 5.67 \times 10^{-8}$ W m$^{-2}$ K$^{-4}$ is the Stefan–Boltzmann constant
and $T$ is the local temperature in Kelvin.

### S2.2 Dynamic component

$$F_{\mathrm{dyn}}(x,t) \propto S(x,t) \cdot C(x,t)$$

where $S \in [0,1]$ is the normalised saturation ratio (relative humidity) and
$C \in [0,1]$ is the cloud fraction. Convective and latent components are physically
coupled — latent heat release drives convection — so combining them into a single effective
flux $F_{\mathrm{dyn}}$ is both physically justified and analytically convenient. Thus:

$$R(T) = \frac{S(T)\cdot C(T)}{\sigma\, T_K^4}$$

---

## S3. The Physical Envelope

### S3.1 Observed constraint

ERA5 observations are not uniformly distributed in $(T, C, R)$ space. They are constrained
to a curved surface with a well-defined upper boundary:

$$R_{\max}(T, C) = \sup\{R \mid (T, R, C) \in \mathcal{M}\}$$

This boundary reflects the physical impossibility of simultaneously achieving high $S$,
high $C$, and low $T$ without triggering convective instability and precipitation feedback.

### S3.2 Parametric model

Motivated by Clausius–Clapeyron scaling and convective instability constraints, we
propose the two-dimensional parametric form:

$$\boxed{R_{\max}(T, C) = \frac{(a \cdot C + b)\cdot e^{-\alpha T}}{\sigma\, T_K^4}}$$

Fitted to the 99th-percentile envelope of 16,346,138 ERA5 observations (1980–2025) on a
$40 \times 20$ temperature–cloud cover grid:

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| $a$ | 0.8744 | cloud modulation coefficient |
| $b$ | 0.0166 | baseline dynamic flux |
| $\alpha$ | 0.0044 | thermodynamic decay rate |

**Fit quality:** RMSE $= 2.66 \times 10^{-4}$, MAE $= 1.91 \times 10^{-4}$, $R^2 = 0.937$

A non-parametric Gaussian Process fit yielded $R^2 = 0.849$, confirming that the analytic
formula captures the dominant physical structure.

### S3.3 Relative position on the manifold

We define the normalised coordinate:

$$\rho(x,t) = \frac{R(x,t)}{R_{\max}(T,C)} \in [0,\, 1+\varepsilon]$$

where $\varepsilon > 0$ represents reanalysis artefacts and measurement uncertainty.
The value $\rho = 1$ corresponds to the physical boundary of the manifold. All regime
boundaries in this framework are expressed as thresholds on $\rho$, making classification
**invariant to the T-dependent radiative scale**.

---

## S4. Regime Classification

Cloud regimes are defined as regions in the $(T, C, \rho)$ state space, not as
categorical labels. Boundaries are continuous transitions, not discrete jumps.
Classification priority (first matching condition):

| Regime | Condition | Physical character |
|--------|-----------|-------------------|
| Clear | $C < 0.20$ | Radiative cooling dominates; no phase activity |
| Convective | $S > 0.9$ and $\rho > 0.75$ | Instability-driven; dominant latent heat release |
| Stratiform | $S > 0.9$ and $\rho \leq 0.75$ | Layered; radiative–diffusive balance |
| Mixed | $S \leq 0.9$ and $C > 0.50$ | Transitional; intermittent cloud formation |
| Transitional | remainder | Boundary states |

**Observed distribution** (ERA5, 1980–2025, $n = 16.3 \times 10^6$):

| Regime | Count | Fraction |
|--------|-------|----------|
| Clear | 5,345,432 | 32.7% |
| Mixed | 6,348,707 | 38.8% |
| Convective | 2,418,221 | 14.8% |
| Transitional | 2,233,778 | 13.7% |
| Stratiform | 0 | 0.0% |

The absence of stratiform events under the current thresholds is noted as a limitation
(see Section S9).

---

## S5. The ρ = 1 Edge: Asymmetry Analysis

**Figure 1** shows the distribution of rise and fall times at the physical boundary.

### S5.1 Method

From time-ordered ERA5 data, we detect complete cycles in the $\rho$ time series of each
spatial grid cell: $\rho < \rho_{\mathrm{base}} \to \rho > \rho_{\mathrm{edge}} \to
\rho < \rho_{\mathrm{base}}$, with $\rho_{\mathrm{edge}} = 0.90$ and
$\rho_{\mathrm{base}} = 0.50$.

### S5.2 Results

From $n = 700{,}269$ complete cycles:

| Metric | Value |
|--------|-------|
| Median rise time $\tau_{\mathrm{rise}}$ | 4.0 ERA5 steps |
| Median fall time $\tau_{\mathrm{fall}}$ | 5.0 ERA5 steps |
| Asymmetry index AI | $+0.111$ |

$$\mathrm{AI} = \frac{\tau_{\mathrm{fall}} - \tau_{\mathrm{rise}}}{\tau_{\mathrm{fall}} + \tau_{\mathrm{rise}}} = +0.111$$

AI $> 0$ confirms the physical expectation: **slow build-up, fast collapse**.

### S5.3 Physical interpretation

The rise distribution is narrow and sharply peaked — the build-up process is constrained.
The fall distribution has a substantially longer tail, extending to $>20$ steps. This
indicates that some convective systems **partially sustain themselves** through latent heat
feedback before finally collapsing — a signature of the non-equilibrium self-organising
character of deep convection.

---

## S6. Conditional Drift Field

**Figure 2** shows the conditional mean drift $E[d\rho/dt \mid \rho, T]$.

### S6.1 Definition

For each $(\rho_{\mathrm{bin}}, T_{\mathrm{bin}})$ cell, we compute the binned mean of
$d\rho/dt$ from consecutive ERA5 time steps within the same spatial cell. This gives an
empirical vector field on the manifold — the analogue of the Fokker–Planck drift term
for the $\rho$ coordinate.

### S6.2 Key result: spontaneous organisation threshold

The drift profile (Figure 2, centre panel) reveals a **zero-crossing at $\rho \approx 0.75$**:

- For $\rho < 0.75$: mean drift is **negative** — the system spontaneously moves away
  from the envelope and returns to lower-energy states
- For $\rho > 0.75$: mean drift is **positive** — the system is attracted toward $\rho = 1$
- For $\rho > 1$: drift turns negative again — the system is repelled back

This structure constitutes a **non-equilibrium potential well** with an unstable fixed
point at $\rho \approx 0.75$. Below this threshold, radiative dissipation dominates.
Above it, the latent heat feedback loop takes over and drives the system toward the
physical boundary.

The noise level (Figure 2, right panel) is maximal near $\rho = 1$ at warm temperatures
($T = 20$–$40°C$), consistent with the known high variability of tropical convection.

---

## S7. Probabilistic State Estimator

**Figure 3** shows example outputs of the state estimator.

### S7.1 Scope and limitations

> **This module provides probability estimates, not weather forecasts.**  
> It reflects historical transition statistics from ERA5 data (1980–2025).  
> It does not predict specific weather events or individual atmospheric outcomes.

### S7.2 Construction

From consecutive ERA5 time steps within each spatial cell, we build:

1. A $\rho$-transition kernel: $P(\rho_{t+1} \mid \rho_t, T_t)$ — the empirical
   distribution of $\rho$ at the next step given the current state
2. A cell regime distribution: $P(\mathrm{regime} \mid \rho_{\mathrm{bin}}, T_{\mathrm{bin}})$

Multi-step estimates are obtained by $k$ repeated applications of the transition kernel.

### S7.3 Selected results

| Initial state | k | Dominant regime | Probability |
|--------------|---|----------------|-------------|
| $\rho=0.30$, $T=10°C$, $C=0.3$ | 1 | clear | 91% |
| $\rho=0.70$, $T=5°C$, $C=0.7$ | 1 | clear | 58% |
| $\rho=0.92$, $T=-5°C$, $C=0.9$ | 1 | mixed | 63% |
| $\rho=0.92$, $T=-5°C$, $C=0.9$ | 6 | clear | 87% |
| $\rho=0.98$, $T=15°C$, $C=0.85$ | 1 | convective | 56% |
| $\rho=0.98$, $T=15°C$, $C=0.85$ | 6 | convective | 91% |

The last two rows demonstrate **convective persistence**: once the system reaches the
physical boundary at warm temperatures, convective conditions are self-reinforcing over
multiple time steps. The cold near-envelope case ($T = -5°C$) shows the opposite — the
system rapidly transitions to clear conditions because latent heat is insufficient to
sustain the convective structure.

---

## S8. Connection to the TMF Dissipation

The effective dissipation of the TMF is extended as:

$$\kappa(\mu) = \kappa_{\mathrm{rad}} + \Phi(R,\, \mu)$$

where $R$ provides a measurable intermediate variable. Clouds act as **co-moving
thermodynamic vortices** in the atmospheric flow — they concentrate and release latent
energy, modify radiative transfer, and generate localised instabilities. Within the TMF,
they form the mesoscopic carriers of geometric structure in $W$-space.

The drift field (Section S6) provides the empirical form of $d\kappa/d\rho$ along the
manifold, connecting the observable $\rho$ coordinate to the theoretical dissipation
structure.

---

## S9. Limitations

**L1. Stratiform regime.** The current classification yields zero stratiform events in the
ERA5 dataset. The threshold conditions require revision — stratiform clouds likely occupy
a narrow band in $(\rho, S)$ space that the present binary threshold does not resolve.
A continuous probabilistic boundary is recommended for future work.

**L2. Temporal resolution.** ERA5 time steps are used as the unit of dynamics. The
asymmetry index AI and the drift field are therefore expressed in ERA5 steps, not physical
time. Converting to hours requires knowledge of the ERA5 output frequency used.

**L3. Spatial averaging.** ERA5 grid cells represent area averages. Sub-grid convective
organisation — the primary driver of deep convection — is not resolved. This may bias
the rise/fall statistics toward longer time scales.

**L4. Single-level approximation.** The analysis uses 2m temperature and total cloud
cover. Vertical structure — cloud-top temperature, cloud thickness, multi-layer systems
— is not captured. The $R$ descriptor should be extended to pressure levels in future work.

**L5. Parametric envelope validity.** The fitted envelope $R_{\max}(T, C)$ achieves
$R^2 = 0.937$ globally but shows systematic residuals at $T < -25°C$ (cold-cloud regime)
and $T > 30°C$ (deep tropical convection). These regions require either a modified
analytic form or the GP-based alternative.

**L6. No causal inference.** The transition kernel and drift field are statistical
summaries of historical co-occurrence patterns. They do not establish causal mechanisms.
The spontaneous organisation threshold at $\rho \approx 0.75$ is an empirical finding
that requires theoretical explanation within the TMF framework.

**L7. Regional heterogeneity not resolved.** The analysis is global. Tropical, mid-latitude,
and polar regimes are pooled. Regional drift fields and transition kernels are expected to
differ substantially. This is an open direction left for future contributions.

**L8. ERA5 reanalysis uncertainty.** Values $\rho > 1$ (observed in approximately 5% of
points) are attributed to reanalysis interpolation artefacts and humidity overshoot errors.
The physical constraint $R \leq R_{\max}(T,C)$ is assumed to hold exactly.

---

## S10. Reproducibility

All results are fully reproducible from publicly available ERA5 data using the accompanying
open-source software.

**Software:** `cloud_manifold_v3_0.py`  
**Repository:** `https://github.com/LesliePi/ClimateManifoldDynamics/tree/main/cloud`  
**Data:** ERA5 reanalysis, variables: `2m_temperature` (t2m), `total_cloud_cover` (tcc),
`2m_dewpoint_temperature` (d2m). Available via the Copernicus Climate Data Store.  
**Period:** 1980–2025, global coverage  
**Sample size:** 16,346,138 observations after spatial subsampling (step = 200)

---

## S11. Position within the TMF Programme

This Supplementary Note provides:

- a physically interpretable observable ($R$ and $\rho$) computable from standard reanalysis
- an analytic description of the thermodynamic state space boundary
- empirical evidence for asymmetric edge dynamics and a spontaneous organisation threshold
- a probabilistic state estimator for regime probability given manifold position
- a pathway toward cloud-resolved TMF extensions without modifying the core formulation

It does not modify the core TMF equations. It extends their interpretability and testability
through a bridge between observable mesoscale structures and the theoretical dissipation term.

> The atmosphere is not characterised by equilibrium states, but by trajectories within a
> constrained, history-dependent manifold — and the $\rho = 1$ edge is where that constraint
> becomes dynamically active.

---

## References

ERA5 reanalysis data: Hersbach, H. et al. (2020). The ERA5 global reanalysis.
*Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049.
https://doi.org/10.1002/qj.3803

TMF framework: Tatai, L. (2026). *Thermodynamic Manifold Dynamics of the Climate System*
(v2.0). Zenodo. https://doi.org/10.5281/zenodo.19430594

Stefan–Boltzmann law and Clausius–Clapeyron relation: standard atmospheric thermodynamics
references (e.g., Wallace & Hobbs, *Atmospheric Science*, 2nd ed., 2006).

---

*End of Supplementary Note S1*
