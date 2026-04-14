# CMD_Analysis_Unified.py – Program Manual

**Version:** 1.0.1 (2026-04-12)  
**Author:** László Tatai  
**ORCID:** 0009-0007-5153-6306  
**Repository:** [ClimateManifoldDynamics](https://github.com/LesliePi/ClimateManifoldDynamics)  
**License:** Apache 2.0 with Commons Clause  

---

## Overview

`CMD_Analysis_Unified.py` is the core analysis script of the Climate Manifold Dynamics (CMD) project. It processes ERA5 reanalysis data (temperature, relative humidity, cloud cover) and computes multiple geometric, topological, and information‑theoretic metrics to characterise the climate system’s state space and its evolution over time.

The script integrates six independent modules:

| Module | Purpose |
|--------|---------|
| `hexbin` | 2D hexbin visualisations (cloud cover, mass, trajectory, anomaly, topology) |
| `diff` | Differential hexbin density change between baseline and comparison windows |
| `volume` | Convex hull volume (annual and 3‑year windows) in normalised (T, RH, cloud) space |
| `kde` | 95% kernel density estimate volume (sliding 3‑year windows, computationally heavy) |
| `shift` | Jensen–Shannon divergence and Wasserstein distances (2D) |
| `3d` | Full 3D metrics: JS divergence, Wasserstein, persistent homology (requires `gudhi`) |

All outputs are saved as PNG figures in a timestamped directory.

---

## Dependencies

### Core
- Python ≥ 3.8
- numpy, pandas, matplotlib
- scipy (spatial, stats, ndimage)
- scikit‑learn (MinMaxScaler)
- argparse, multiprocessing (optional)

### Optional (but recommended)
- `gudhi` – for persistent homology (3D module)
- `psutil` – for accurate physical CPU core count (if parallel KDE is used)

### Installation (conda example)
```bash
conda create -n cmd_env python=3.10
conda activate cmd_env
conda install numpy pandas matplotlib scipy scikit-learn
conda install -c conda-forge gudhi   # optional
pip install psutil                    # optional
Command‑Line Arguments
Argument	Type	Default	Description
--base Y0 Y1	int, int	1980 1982	Baseline window (first and last year)
--comp Y0 Y1	int, int	2022 2024	Comparison window
--modules	list	hexbin diff volume	Which modules to run (choices: hexbin, diff, volume, kde, shift, 3d)
--parquet PATH	str	output/era5_processed.parquet	Input Parquet file (ERA5 preprocessed)
--outdir PATH	str	output	Output directory for figures
--grid N	int	50	Hexbin grid size (number of bins along one axis)
--kde-res N	int	40	3D histogram / KDE resolution per dimension (e.g., 40 → 40³ cells)
--kde-workers N	int	None	Number of parallel workers for KDE (Windows: set 1)
--no-show	flag	False	Disable interactive matplotlib windows (batch mode)
Example runs
bash

# Default: hexbin + diff + volume (fast)
python CMD_Analysis_Unified.py

# Add KDE (slow) and 3D metrics
python CMD_Analysis_Unified.py --modules hexbin diff volume kde 3d

# Change baseline and comparison windows
python CMD_Analysis_Unified.py --base 1995 1997 --comp 2015 2017

# Batch mode (no pop‑up figures)
python CMD_Analysis_Unified.py --no-show

Module Descriptions
1. hexbin

    2a. Cloud cover: 3‑year windows (1980–2024, 15 windows). Hexbin of T vs RH, colour = mean cloud cover.

    2b. Mass (point count): Same windows, colour = log(count) per bin.

    2c. Trajectory: Annual centroid (mean T, RH) with colour = cloud cover.

    2d. Temperature anomaly: Annual centroid T vs T anomaly (baseline 1980).

    2e. Topology (baseline vs comparison): Common T‑RH extent, mincnt=0 to ensure identical grid. Identifies stable, lost, and new bins.

2. diff

    Differential hexbin between baseline and comparison windows.

    Computes normalised density (counts / max count) for each window, then Δ = comp_norm - base_norm.

    Visualised with diverging colormap (RdBu_r) and TwoSlopeNorm centered at 0.

    Also produces a summary figure with differential map + topological mask.

3. volume

    Annual convex hull volume in normalised (T_norm, RH_norm, cloud_norm) space.

    3D convex hull for baseline vs comparison windows (with scatter and hull edges).

    Reports relative change in volume.

4. kde

    Sliding 3‑year windows (1980–1982, 1983–1985, …, 2022–2024).

    For each window: subsample (max 200,000 points), build 3D Gaussian KDE, evaluate on a regular grid (kde_res³), compute the 95% quantile volume.

    Note: On Windows, parallel execution may fail; use --kde-workers 1 or run serially (the script now defaults to serial for stability).

5. shift

    2D histogram (T, RH) with 50×50 bins, common extent.

    Jensen–Shannon divergence (0–1 scale).

    Wasserstein‑1 distances for T, RH, cloud marginals.

    Contour plot of the two smoothed distributions.

6. 3d

    3D histogram (T_norm, RH_norm, cloud_norm) with resolution kde_res³.

    Jensen–Shannon divergence.

    Combined Wasserstein distance (Euclidean norm of the three marginal Wasserstein distances).

    Persistent homology (if gudhi is installed): cubical complex built from the negative smoothed density field; bottleneck distances for H0 (connected components) and H1 (loops); persistence diagrams plotted.

    2D projection plots (T‑RH, T‑cloud, RH‑cloud) with smoothed density contours.

Data Input Format

The script expects a Parquet file with at least the following columns:
Column	Description	Unit / Range
T	2‑metre temperature	°C
RH	Relative humidity	0–1 (clipped)
cloud	Total cloud cover	0–1 (clipped)
year	Year	integer (e.g., 1980)

The script automatically adds:

    T_anom = T – mean(T in 1980)

    T_norm, RH_norm, cloud_norm via MinMaxScaler (fitted on the whole dataset).

Output Files

All figures are saved with a timestamp prefix (YYYYMMDD_HHMMSS_). Typical outputs:
Suffix	Module	Description
hexbin_cloud.png	hexbin	15‑panel cloud cover
hexbin_mass.png	hexbin	15‑panel point count
trajectory.png	hexbin	Annual centroid trajectory
temp_anomaly.png	hexbin	T vs T anomaly
topology.png	hexbin	Topological mask (baseline vs comp)
diff_hexbin_...png	diff	3‑panel differential (base, comp, Δ)
diff_summary_...png	diff	Summary (Δ map + topology)
volume_annual.png	volume	Annual convex hull volume time series
volume_3d.png	volume	3D hulls (baseline vs comp)
kde_volume.png	kde	Sliding‑window KDE 95% volume
distribution_shift.png	shift	Contour plot (2D) + JS/Wasserstein
persistence_diagrams.png	3d	H0 & H1 persistence diagrams
3d_projections.png	3d	2D projections of 3D density
Performance Considerations

    Fast modules (hexbin, diff, volume, shift): finish in seconds.

    KDE module: ~20–30 minutes for 43 windows with kde_res=40 on a modern CPU (serial). Parallelisation on Linux/macOS works; Windows users should run serially or with --kde-workers 1.

    3D persistent homology: adds ~1–2 minutes per window pair (only two windows, not the sliding series). Requires gudhi.

To speed up KDE, reduce --kde-res (e.g., 30 → 27,000 grid points instead of 64,000).
Troubleshooting
Problem	Likely cause	Solution
AttributeError: module 'multiprocessing.spawn' has no attribute 'get_preparation_data'	Windows + parallel KDE	Use --kde-workers 1 or the serial KDE version (script now includes serial fallback)
ValueError: operands could not be broadcast together	Different hexbin grid sizes	Fixed by using common extent and mincnt=0 in topology/diff modules
ImportError: No module named 'gudhi'	Missing persistent homology library	Install with conda install -c conda-forge gudhi or skip the 3d module
Memory error during KDE	Grid too large	Reduce --kde-res to 30 or 35
Code Architecture
text

CMD_Analysis_Unified.py
├── parse_args()
├── helpers: ts(), save(), dark_ax(), load_data(), window_df(), window_label()
├── run_hexbin()
├── run_diff()
├── run_volume()
├── run_kde()          # serial by default (Windows safe)
├── run_shift()
├── run_3d_metrics()   # requires gudhi for persistence
└── main()

All modules are independent; the main() function calls them in the order specified by --modules.
Extending the Script

To add a new module:

    Define a function run_newmodule(df, args).

    Add its name to the choices list in parse_args().

    Add a conditional call in main().

The function receives the full DataFrame and the parsed arguments. It may produce figures (use save()) and print console output.
References (Programmatic)

    Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. Computing in Science & Engineering, 9(3), 90–95.

    Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing. Nature Methods, 17, 261–272.

    Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python. Journal of Machine Learning Research, 12, 2825–2830.

    The GUDHI Project (2023). GUDHI user and reference manual. https://gudhi.inria.fr/

License and Citation

This script is part of the Climate Manifold Dynamics project. If you use it in published work, please cite:

Tatai, L. (2026). Climate Manifold Dynamics v3.0: Unified analysis framework (Version 1.0.1) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.19430594

(DOI will be updated for v3.0.)
text


---
## 2. Matematikai és elméleti dokumentáció (CMD_Mathematical_Foundations.md)

```markdown
# Mathematical Foundations of the Climate Manifold Dynamics (CMD) Framework

**Author:** László Tatai  
**ORCID:** 0009-0007-5153-6306  
**Date:** April 2026  
**Related theory DOI:** https://doi.org/10.5281/zenodo.19430594 (v2.0)  

---

## Abstract

This document provides a self‑contained mathematical description of the methods implemented in `CMD_Analysis_Unified.py`. It covers the theoretical underpinnings of hexbin density estimation, convex hull volumes, kernel density estimation (KDE), Jensen–Shannon divergence, Wasserstein distance, and persistent homology – all applied to the climate state space spanned by temperature (T), relative humidity (RH), and cloud cover (CF). Connections to the Thermodynamic Manifold Framework (TMF) are highlighted, and key references are given in APA format.

---

## 1. The Climate State Space

Let the climate system be described by a state vector  

\[
\mathbf{W}(t) = \big(T(t),\, p(t),\, e(t),\, q(t),\, \text{CF}(t),\, \text{ICE}(t),\, \text{OHC}(t)\big)
\]

as introduced in Tatai (2026a). In the current implementation (v3.0), the accessible variables are temperature \(T\) (K or °C), relative humidity \(\text{RH} \in [0,1]\), and total cloud cover \(\text{CF} \in [0,1]\). These are normalised to the unit cube \([0,1]^3\) via MinMax scaling to ensure comparability across variables:

\[
T_{\text{norm}} = \frac{T - T_{\min}}{T_{\max} - T_{\min}}, \quad
\text{RH}_{\text{norm}} = \text{RH}, \quad
\text{CF}_{\text{norm}} = \text{CF}.
\]

The empirical distribution of \(\mathbf{W}\) is represented by a set of \(N\) samples \(\{\mathbf{w}_i\}_{i=1}^N\) obtained from ERA5 reanalysis (Hersbach et al., 2020).

---

## 2. Hexbin Density Estimation

Hexbin binning (Carr et al., 1987) partitions the 2D plane (T, RH) into hexagonal cells. For a given grid size \(G\), the number of cells is approximately \(G^2\). For each cell \(j\) we compute:

- **Mass:** \(m_j = \#\{\mathbf{w}_i \text{ in cell } j\}\)
- **Mean cloud:** \(\bar{c}_j = \frac{1}{m_j}\sum_{\mathbf{w}_i \in \text{cell } j} \text{CF}_i\)

The hexagonal grid is more isotropic than square grids and reduces directional bias.

**Topological comparison** between two windows uses a common extent and `mincnt=0` to guarantee identical cell sets. The binary mask \(M_j = \mathbb{1}(m_j > 0)\) indicates occupied cells. The overlap between baseline (B) and comparison (C) is:

\[
\text{Stable} = M_B \land M_C,\quad
\text{Lost} = M_B \land \lnot M_C,\quad
\text{New} = M_C \land \lnot M_B.
\]

**Differential density** normalises each window’s counts by its maximum count:

\[
\rho_B = \frac{m_j^{(B)}}{\max m^{(B)}},\quad
\rho_C = \frac{m_j^{(C)}}{\max m^{(C)}},\quad
\Delta_j = \rho_C - \rho_B \in [-1,1].
\]

This visualises where the comparison window has relatively higher or lower occupancy.

---

## 3. Convex Hull Volume

Given a set of points in \(\mathbb{R}^3\), the convex hull is the smallest convex set containing them. The volume of the convex hull is computed using the Quickhull algorithm (Barber et al., 1996) as implemented in `scipy.spatial.ConvexHull` (Virtanen et al., 2020).

For a set of \(n\) points, the volume is the sum of the volumes of tetrahedra formed by the hull’s facets and an interior reference point. The volume in normalised coordinates quantifies the **extent** of the climate state space occupied in a given year or window. A larger hull indicates greater variability.

---

## 4. Kernel Density Estimation (KDE) and 95% Volume

The Gaussian KDE estimates the probability density function from a sample \(\{\mathbf{w}_i\}_{i=1}^n\) (Rosenblatt, 1956; Parzen, 1962):

\[
\hat{f}(\mathbf{w}) = \frac{1}{n h^3} \sum_{i=1}^n K\!\left(\frac{\mathbf{w} - \mathbf{w}_i}{h}\right),
\]

where \(K\) is the standard Gaussian kernel and \(h\) is the bandwidth. The optimal bandwidth is estimated by Scott’s rule (Scott, 1992) for multivariate data.

The 95% **volume** is defined as the Lebesgue measure of the smallest region \(R\) such that

\[
\int_R \hat{f}(\mathbf{w})\,d\mathbf{w} = 0.95.
\]

In practice, we evaluate \(\hat{f}\) on a regular grid of resolution \(R^3\) (e.g., \(40^3\)), sort the density values, and find the threshold \(t\) that captures 95% of the total probability mass. The volume is then

\[
V_{95} = \frac{\#\{\hat{f}(\mathbf{w}_k) > t\}}{R^3}.
\]

This metric is less sensitive to outliers than the convex hull volume and better reflects the “typical” occupied region.

---

## 5. Information‑Theoretic Distances

### 5.1 Jensen–Shannon Divergence

Given two probability distributions \(P\) and \(Q\) over a discrete sample space \(\mathcal{X}\), the Jensen–Shannon divergence (Lin, 1991) is

\[
\text{JSD}(P \| Q) = \frac{1}{2} D_{\text{KL}}(P \| M) + \frac{1}{2} D_{\text{KL}}(Q \| M),
\]

where \(M = (P+Q)/2\) is the mixture distribution and \(D_{\text{KL}}\) is the Kullback–Leibler divergence:

\[
D_{\text{KL}}(P \| M) = \sum_{x\in\mathcal{X}} P(x) \log\frac{P(x)}{M(x)}.
\]

JSD is symmetric, bounded between 0 and 1 (when using the natural logarithm, or \(\log_2\) gives bits), and does not require absolute continuity. In the script, we compute JSD from 3D histograms of the normalised state space using `scipy.spatial.distance.jensenshannon` (which returns the square root of JSD, i.e., the Jensen–Shannon distance). The value \(0\) indicates identical distributions; \(1\) indicates completely disjoint support.

### 5.2 Wasserstein‑1 Distance (Earth Mover’s Distance)

For two one‑dimensional empirical distributions, the Wasserstein‑1 distance (also called the Earth Mover’s Distance) is the integral of the absolute difference between their cumulative distribution functions (Villani, 2009):

\[
W_1(P,Q) = \int_{-\infty}^{\infty} |F_P(x) - F_Q(x)|\,dx,
\]

where \(F_P\) and \(F_Q\) are the CDFs. For samples, it equals the average absolute difference between sorted samples. In the script, we compute \(W_1\) separately for \(T\), RH, and CF using `scipy.stats.wasserstein_distance`. The combined distance is defined as the Euclidean norm:

\[
W_{\text{comb}} = \sqrt{W_1(T)^2 + W_1(\text{RH})^2 + W_1(\text{CF})^2}.
\]

This gives a single scalar representing the overall shift in the distribution’s location.

---

## 6. Persistent Homology (Topological Data Analysis)

Persistent homology (Edelsbrunner et al., 2002; Zomorodian & Carlsson, 2005) extracts topological features – connected components (H₀), loops (H₁), voids (H₂) – from a filtered space.

**Cubical complex approach:** We discretise the density field \(\hat{f}\) on a 3D grid (size \(R^3\)). The negative smoothed density \(- \tilde{f}\) (where \(\tilde{f}\) is a Gaussian‑smoothed version) serves as a height function. The sublevel sets

\[
L_t = \{\mathbf{x} : -\tilde{f}(\mathbf{x}) \le t\}
\]

are examined as \(t\) increases. Each topological feature appears at a birth value \(b\) and disappears at a death value \(d\). The persistence \(d-b\) measures how “robust” the feature is.

**Persistence diagrams** plot points \((b, d)\) for each feature. The bottleneck distance between two diagrams (Cohen‑Steiner et al., 2007) is:

\[
d_B(D_1, D_2) = \inf_{\gamma} \sup_{p \in D_1} \|p - \gamma(p)\|_\infty,
\]

where the infimum is over all bijections \(\gamma\) between the points of \(D_1\) and \(D_2\) (allowing matches to the diagonal). In the script, we compute bottleneck distances separately for H₀ and H₁ and then combine them quadratically:

\[
d_{\text{topo}} = \sqrt{ d_B(H_0)^2 + d_B(H_1)^2 }.
\]

A larger value indicates a greater topological change between the two windows.

This method is implemented using the `gudhi` library (The GUDHI Project, 2023).

---

## 7. Connection to the Thermodynamic Manifold Framework (TMF)

The empirical state‑space density \(\mu(\mathbf{W})\) (Tatai, 2026b) can be interpreted as a stability landscape via the potential

\[
V(\mathbf{W}) = -\log \mu(\mathbf{W}).
\]

Minima of \(V\) correspond to frequently visited (stable) states; saddle points indicate regime boundaries. The geometric measures computed above – fractal dimension of the λ–κ manifold, curvature, and persistent homology – are all signatures of this underlying landscape.

The effective dissipation \(\kappa(\mu)\) in the TMF energy equation

\[
\frac{dE}{dt} = F_{\text{in}} - \kappa(\mu) E
\]

is hypothesised to depend on the entire distribution \(\mu\), not only on mean temperature. The residual mutual information \(\text{MI}(\lambda, \kappa \mid T) = 0.185\) nats (Tatai, 2026a) empirically supports this claim.

---

## 8. Summary of Computed Metrics and Their Meaning

| Metric | Symbol | Range | Interpretation |
|--------|--------|-------|----------------|
| Hexbin mass | \(m_j\) | ℕ | Local data density |
| Hexbin Δ density | \(\Delta_j\) | [-1,1] | Relative occupancy change |
| Convex hull volume | \(V_{\text{hull}}\) | [0,1] | Extent of state space |
| KDE 95% volume | \(V_{95}\) | [0,1] | Typical occupied volume |
| Jensen–Shannon distance | JSD | [0,1] | Distributional difference |
| Wasserstein combined | \(W_{\text{comb}}\) | [0,√3] | Average shift in state space |
| Bottleneck distance | \(d_B\) | [0,∞) | Topological change |

---

## References (APA format)

Barber, C. B., Dobkin, D. P., & Huhdanpaa, H. (1996). The quickhull algorithm for convex hulls. *ACM Transactions on Mathematical Software*, 22(4), 469–483. https://doi.org/10.1145/235815.235821

Carr, D. B., Littlefield, R. J., Nicholson, W. L., & Littlefield, J. S. (1987). Scatterplot matrix techniques for large N. *Journal of the American Statistical Association*, 82(398), 424–436. https://doi.org/10.1080/01621459.1987.10478445

Cohen‑Steiner, D., Edelsbrunner, H., & Harer, J. (2007). Stability of persistence diagrams. *Discrete & Computational Geometry*, 37(1), 103–120. https://doi.org/10.1007/s00454-006-1276-5

Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2002). Topological persistence and simplification. *Discrete & Computational Geometry*, 28(4), 511–533. https://doi.org/10.1007/s00454-002-2885-2

Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‑Sabater, J., … Thépaut, J.‑N. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049. https://doi.org/10.1002/qj.3803

Lin, J. (1991). Divergence measures based on the Shannon entropy. *IEEE Transactions on Information Theory*, 37(1), 145–151. https://doi.org/10.1109/18.61115

Parzen, E. (1962). On estimation of a probability density function and mode. *The Annals of Mathematical Statistics*, 33(3), 1065–1076. https://doi.org/10.1214/aoms/1177704472

Rosenblatt, M. (1956). Remarks on some nonparametric estimates of a density function. *The Annals of Mathematical Statistics*, 27(3), 832–837. https://doi.org/10.1214/aoms/1177728190

Scott, D. W. (1992). *Multivariate density estimation: Theory, practice, and visualization*. John Wiley & Sons. https://doi.org/10.1002/9780470316849

Tatai, L. (2026a). *Climate feedback manifold: A geometric approach to time‑varying λ–κ dynamics (1880–2025)* (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.19421325

Tatai, L. (2026b). *State‑dependent climate feedback dynamics: Forward simulation and κ calibration* (v2.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.19430594

The GUDHI Project. (2023). *GUDHI user and reference manual*. https://gudhi.inria.fr/

Villani, C. (2009). *Optimal transport: Old and new* (Vol. 338). Springer. https://doi.org/10.1007/978-3-540-71050-9

Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., … SciPy 1.0 Contributors. (2020). SciPy 1.0: fundamental algorithms for scientific computing. *Nature Methods*, 17, 261–272. https://doi.org/10.1038/s41592-019-0686-2

Zomorodian, A., & Carlsson, G. (2005). Computing persistent homology. *Discrete & Computational Geometry*, 33(2), 249–274. https://doi.org/10.1007/s00454-004-1146-y

---

*This document is part of the Climate Manifold Dynamics project. Version 1.0 – April 2026.*

