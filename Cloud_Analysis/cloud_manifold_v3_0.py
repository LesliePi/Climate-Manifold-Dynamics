# Climate Manifold Dynamic - Cloud Manifold Pipeline
# cloud_manifold_v3_0.py
# Variables: 2m temperature (t2m) + total cloud cover (tcc) + relative humidity (d2m / q)
# Version: v3.0  (2026-04-14)
# Author: László Tatai
# ORCID:  0009-0007-5153-6306
# Part of: Climate Manifold Dynamics
# GitHub:  https://github.com/LesliePi/ClimateManifoldDynamics
# Theory:  https://doi.org/10.5281/zenodo.19430594  (v2.0)
# License: Apache License 2.0 WITH Commons Clause v1.0
#
# Licensed under the Apache License, Version 2.0 (the "License")
# with the addition of the Commons Clause License Condition v1.0.
# You may not use this file except in compliance with the License.
# You may obtain a copy of the Apache License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# The Commons Clause condition (full text below) applies:
#
# ──────────────────────────────────────────────────────────────
# "Commons Clause" License Condition v1.0
#
# Without limiting other conditions in the License, the grant of
# rights under the License will not include, and the License does not
# grant to you, the right to Sell the Software.
#
# Software: Climate Manifold Dynamic -- cloud_manifold_v2_0.py
# License:  Apache License 2.0
# Licensor: László Tatai
# ──────────────────────────────────────────────────────────────
# -*- coding: utf-8 -*-

# ============================================================
# CLIMATE MANIFOLD DYNAMICS – CLOUD MANIFOLD PIPELINE  v3.0
#
# v2.0 additions (retained):
#   Section 5  – Surface fit: parametric R_max(T,C) + GP comparison
#   Section 6  – ρ-based regime classification
#   Section 7  – P(regime | T,C) probabilistic module + entropy map
#
# v3.0 additions — Edge Dynamics:
#   Section 8  – Asymmetry analysis at ρ = 1
#     8.1  Build-up / collapse event detection in time-ordered data
#     8.2  Rise time vs fall time distributions per (T,C) cell
#     8.3  Asymmetry index  AI = (τ_fall − τ_rise) / (τ_fall + τ_rise)
#          AI > 0 → slow build, fast collapse (expected physical result)
#
#   Section 9  – Conditional drift field  dρ/dt | ρ, T, C
#     9.1  Compute ρ time-derivative from consecutive ERA5 time steps
#     9.2  Binned mean drift: E[dρ/dt | ρ_bin, T_bin] on 2D grid
#     9.3  Drift sign map: where does the manifold push the system?
#          Positive drift → system moves toward envelope (ρ→1)
#          Negative drift → system relaxes away from envelope
#
#   Section 10 – Probabilistic state estimator
#     10.1  Transition kernel: P(ρ_{t+k} | ρ_t, T_t, C_t)
#           Built from observed consecutive-step ρ transitions
#     10.2  Regime transition matrix: P(regime_{t+k} | regime_t, T, C)
#     10.3  query_state(rho, T, C, k) → regime probability vector
#           This is the weather probability estimator:
#           "Given current manifold position, what is the likely
#            regime distribution k steps ahead?"
#           NOTE: This is a probability estimate, NOT a forecast.
#                 It does not predict specific weather events.
# ============================================================

import os
import glob
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import xarray as xr

from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde

def _r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """R² coefficient of determination."""
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    return 1.0 - ss_res / (ss_tot + 1e-30)


# ── Minimal pure-NumPy GP (RBF kernel, closed-form posterior mean) ──────────
#
# This avoids the sklearn dependency for GP fitting, making the code portable
# to any NumPy environment.  The implementation is deliberately minimal:
# only the posterior mean is computed (no uncertainty estimate).
#
# Kernel: k(x,x') = σ²·exp(−0.5·||x−x'||²/l²) + noise·δ
# Optimisation: marginal log-likelihood maximised via grid search over (l, σ, noise).

class _NumpyGP:
    """
    Pure-NumPy Gaussian Process regressor with RBF kernel.

    Parameters
    ----------
    l_init     : initial length-scale
    sigma_init : initial signal variance
    noise      : observation noise variance (fixed)
    """

    def __init__(
        self,
        l_init:     float = 1.0,
        sigma_init: float = 1.0,
        noise:      float = 1e-3,
    ):
        self.l     = l_init
        self.sigma = sigma_init
        self.noise = noise
        self._alpha  = None   # (K + noise·I)⁻¹ y
        self._X_train = None
        self._y_mean  = 0.0
        self._X_scale = None

    @staticmethod
    def _rbf(X1: np.ndarray, X2: np.ndarray, l: float, sigma: float) -> np.ndarray:
        """Compute RBF kernel matrix K(X1, X2)."""
        # X1: (n,d), X2: (m,d) → (n,m)
        diff  = X1[:, np.newaxis, :] - X2[np.newaxis, :, :]   # (n,m,d)
        sq    = np.sum(diff**2, axis=-1)                        # (n,m)
        return sigma**2 * np.exp(-0.5 * sq / (l**2 + 1e-30))

    def _log_marginal(self, X: np.ndarray, y: np.ndarray, l: float, sigma: float) -> float:
        n = len(y)
        K = self._rbf(X, X, l, sigma)
        K_noise = K + (self.noise + 1e-8) * np.eye(n)
        try:
            L    = np.linalg.cholesky(K_noise)
            alph = np.linalg.solve(L.T, np.linalg.solve(L, y))
            lml  = -0.5 * y @ alph - np.sum(np.log(np.diag(L))) - 0.5 * n * np.log(2 * np.pi)
            return float(lml)
        except np.linalg.LinAlgError:
            return -np.inf

    def fit(self, X: np.ndarray, y: np.ndarray) -> None:
        """
        Fit GP by optimising (l, sigma) via coarse grid search on log-marginal likelihood.
        X must be standardised (zero mean, unit variance) before calling.
        """
        self._y_mean  = float(np.mean(y))
        y_c = y - self._y_mean

        # Grid search over length-scale and signal std
        best_lml, best_l, best_s = -np.inf, self.l, self.sigma
        for l in np.logspace(-1, 1, 8):
            for s in np.logspace(-2, 1, 6):
                lml = self._log_marginal(X, y_c, l, s)
                if lml > best_lml:
                    best_lml, best_l, best_s = lml, l, s

        self.l, self.sigma = best_l, best_s

        K = self._rbf(X, X, self.l, self.sigma)
        K_noise = K + (self.noise + 1e-8) * np.eye(len(y_c))
        L        = np.linalg.cholesky(K_noise + 1e-10 * np.eye(len(y_c)))
        self._alpha   = np.linalg.solve(L.T, np.linalg.solve(L, y_c))
        self._X_train = X.copy()

    def predict(self, X_new: np.ndarray) -> np.ndarray:
        """Return posterior mean at X_new."""
        K_star = self._rbf(X_new, self._X_train, self.l, self.sigma)  # (m, n)
        return K_star @ self._alpha + self._y_mean

warnings.filterwarnings("ignore")


# ============================================================
# 0. LOG BUFFER
# ============================================================

LOG_BUFFER: list[str] = []


def log(msg: str) -> None:
    """Print to console AND buffer for Markdown export."""
    print(msg)
    LOG_BUFFER.append(str(msg))


# ============================================================
# 1. CONFIG
# ============================================================

DATA_PATH       = "data/raw/era5/*.nc"
OUTPUT_DIR      = "output"
PROCESSED_FILE  = os.path.join(OUTPUT_DIR, "era5_cloud_processed.parquet")
FORCE_REPROCESS = False
COMPRESSION     = "zstd"
SAMPLE_STEP     = 200        # spatial subsampling of flattened ERA5 arrays
P0              = 101325.0   # reference surface pressure [Pa]

# Stefan–Boltzmann constant
SIGMA = 5.67e-8  # W m⁻² K⁻⁴

# Grid resolution for surface fitting and probability maps
T_BINS = 40   # temperature bins  (typically −45 … +50 °C)
C_BINS = 20   # cloud cover bins  (0 … 1)

# ρ-based regime thresholds (fraction of local R_max)
RHO_CONVECTIVE  = 0.75   # ρ > this AND S > S_SAT  → convective
S_SAT_THRESHOLD = 0.90   # saturation threshold for S
CLOUD_MIXED_THR = 0.50   # cloud > this (when S ≤ S_SAT) → mixed
CLOUD_CLEAR_THR = 0.20   # cloud < this → clear

# GP fitting: max number of envelope grid points used for GP training
GP_MAX_TRAIN    = 400

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# 2. UTILITIES
# ============================================================

def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def save_figure(fig: plt.Figure, name: str) -> str:
    """Save figure to OUTPUT_DIR and return the relative filename."""
    fname = f"{_ts()}_{name}.png"
    path  = os.path.join(OUTPUT_DIR, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    log(f"  Saved: {path}")
    return fname


# ============================================================
# 3. PHYSICS
# ============================================================

def saturation_pressure(T_K: np.ndarray) -> np.ndarray:
    """Magnus formula: saturation vapour pressure [Pa] from temperature [K]."""
    T_C = np.asarray(T_K, dtype=float) - 273.15
    return 6.112 * np.exp((17.67 * T_C) / (T_C + 243.5)) * 100.0


def rh_from_specific_humidity(
    q: np.ndarray,
    T_K: np.ndarray,
    p: np.ndarray | float = P0,
) -> np.ndarray:
    q   = np.asarray(q, dtype=float)
    T_K = np.asarray(T_K, dtype=float)
    e   = (q * p) / (0.622 + 0.378 * q)
    return np.clip(e / saturation_pressure(T_K), 0.0, 1.5)


def rh_from_dewpoint(Td_K: np.ndarray, T_K: np.ndarray) -> np.ndarray:
    return np.clip(
        saturation_pressure(Td_K) / saturation_pressure(T_K),
        0.0, 1.5,
    )


def compute_R(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add F_rad, S, F_dyn, R columns.

    R = F_dyn / F_rad  =  (S · C) / (σ · T_K⁴)

    where S = RH clipped to [0, 1.5].
    """
    df    = df.copy()
    T_K   = df["T"].values + 273.15
    df["F_rad"] = SIGMA * T_K**4
    df["S"]     = df["RH"].clip(0.0, 1.5)
    df["F_dyn"] = df["S"] * df["cloud"]
    df["R"]     = df["F_dyn"] / (df["F_rad"] + 1e-30)
    return df


# ============================================================
# 4. ERA5 LOAD
# ============================================================

def _load_era5_raw(data_path: str, sample_step: int) -> pd.DataFrame:
    files = sorted(glob.glob(data_path))
    log(f"\nFiles found: {len(files)}")
    if not files:
        raise FileNotFoundError(f"No .nc files at: {data_path}")

    VAR_RENAME = {
        "2m_temperature":          "t2m",
        "total_cloud_cover":       "tcc",
        "surface_pressure":        "sp",
        "2m_dewpoint_temperature": "d2m",
        "specific_humidity":       "q",
    }

    chunks: list[pd.DataFrame] = []

    for i, fpath in enumerate(files):
        fname = os.path.basename(fpath)
        log(f"  [{i+1:3d}/{len(files)}] {fname}")

        try:
            ds = xr.open_dataset(fpath)
        except Exception as exc:
            log(f"    -> SKIP ({exc})")
            continue

        rename = {old: new for old, new in VAR_RENAME.items() if old in ds}
        if rename:
            ds = ds.rename(rename)

        if "t2m" not in ds or "tcc" not in ds:
            log("    -> SKIP (missing t2m or tcc)")
            ds.close()
            continue

        time_coord = "valid_time" if "valid_time" in ds.coords else "time"
        times = pd.to_datetime(ds[time_coord].values)
        lats  = ds.latitude.values
        lons  = ds.longitude.values
        ntim, nlat, nlon = len(times), len(lats), len(lons)
        N = ntim * nlat * nlon

        t2m_f = ds["t2m"].values.reshape(-1)
        tcc_f = ds["tcc"].values.reshape(-1)
        lat_f = np.broadcast_to(
            lats[np.newaxis, :, np.newaxis], (ntim, nlat, nlon)).reshape(-1)
        lon_f = np.broadcast_to(
            lons[np.newaxis, np.newaxis, :], (ntim, nlat, nlon)).reshape(-1)
        tim_f = np.broadcast_to(
            np.array(times)[:, np.newaxis, np.newaxis], (ntim, nlat, nlon)).reshape(-1)

        idx = np.arange(0, N, sample_step)

        t2m_s = t2m_f[idx]
        tcc_s = tcc_f[idx]
        lat_s = lat_f[idx]
        lon_s = lon_f[idx]
        tim_s = tim_f[idx]

        if "d2m" in ds:
            d2m_s  = ds["d2m"].values.reshape(-1)[idx]
            rh_s   = rh_from_dewpoint(d2m_s, t2m_s)
            rh_src = "dewpoint (d2m)"
        elif "q" in ds:
            q_s    = ds["q"].values.reshape(-1)[idx]
            p_s    = (ds["sp"].values.reshape(-1)[idx]
                      if "sp" in ds else np.full(len(idx), P0))
            rh_s   = rh_from_specific_humidity(q_s, t2m_s, p_s)
            rh_src = "specific humidity (q)"
        else:
            rh_s   = np.full(len(idx), np.nan)
            rh_src = "none"

        log(f"    -> {len(idx):,} pts | RH source: {rh_src}")
        ds.close()

        chunks.append(pd.DataFrame({
            "T":     t2m_s - 273.15,
            "cloud": tcc_s,
            "RH":    rh_s,
            "lat":   lat_s,
            "lon":   lon_s,
            "time":  pd.to_datetime(tim_s),
        }))

    df = pd.concat(chunks, ignore_index=True)
    df = df.dropna(subset=["T", "cloud"])
    df["month"]  = df["time"].dt.month
    df["year"]   = df["time"].dt.year

    n_high = (df["RH"] > 1.2).sum()
    n_low  = (df["RH"] < 0.0).sum()
    log(f"\nDataFrame: {len(df):,} rows | years {df['year'].min()}–{df['year'].max()}"
        f" | RH valid: {df['RH'].notna().sum():,}")
    log(f"RH validation: RH > 1.2 = {n_high:,}  |  RH < 0 = {n_low:,}")
    return df


def load_era5(
    data_path: str   = DATA_PATH,
    sample_step: int = SAMPLE_STEP,
    force: bool      = FORCE_REPROCESS,
) -> pd.DataFrame:
    """Load ERA5 with Parquet cache. Set force=True to reprocess raw .nc files."""
    if not force and os.path.exists(PROCESSED_FILE):
        log(f"Loading cached Parquet: {PROCESSED_FILE}")
        df = pd.read_parquet(PROCESSED_FILE)
        log(f"  Loaded {len(df):,} rows | years {df['year'].min()}–{df['year'].max()}")
        return df

    log("Processing raw .nc files ...")
    df = _load_era5_raw(data_path, sample_step)
    df.to_parquet(PROCESSED_FILE, compression=COMPRESSION, index=False, engine="pyarrow")
    mb = os.path.getsize(PROCESSED_FILE) / 1024 / 1024
    log(f"Parquet cache saved: {PROCESSED_FILE} ({mb:.1f} MB)")
    return df


# ============================================================
# 5. SURFACE FIT MODULE
# ============================================================
#
# The 3D scatter plot shows that the cloud regime data lies on a
# curved surface in (T, C, R) space.  The upper boundary of this
# surface — R_max(T, C) — is the physical envelope.
#
# We estimate R_max two ways and compare:
#   A) Parametric:  R_max(T,C) = (a·C + b)·exp(−α·T) / σ·T_K⁴
#   B) GP:          GaussianProcessRegressor with RBF kernel
#
# The relative position  ρ = R / R_max(T,C)  normalises each
# observation to [0,1] and makes regime thresholds scale-invariant.
# ============================================================

class SurfaceFit:
    """
    Container for the 2D physical envelope R_max(T, C).

    Attributes
    ----------
    grid_T, grid_C, grid_Rmax : 1D arrays — observed envelope grid
    params_parametric          : (a, b, alpha) from parametric fit
    gp_model, gp_scaler        : trained GP regressor + input scaler
    metrics                    : dict with RMSE / MAE / R² for both models
    """

    def __init__(self):
        self.grid_T    = None
        self.grid_C    = None
        self.grid_Rmax = None
        self.params_parametric = None
        self.gp_model  = None
        self._gp_X_mean = None
        self._gp_X_std  = None
        self.metrics   = {}

    # ----------------------------------------------------------
    # 5.1  Observed envelope grid
    # ----------------------------------------------------------

    def build_envelope_grid(self, df: pd.DataFrame) -> None:
        """
        Compute 99th-percentile R in each (T_bin, C_bin) cell.

        Populates self.grid_T, self.grid_C, self.grid_Rmax.
        """
        T_edges = np.linspace(df["T"].min(),     df["T"].max(),     T_BINS + 1)
        C_edges = np.linspace(df["cloud"].min(), df["cloud"].max(), C_BINS + 1)

        T_vals, C_vals, R_vals = [], [], []

        for i in range(len(T_edges) - 1):
            for j in range(len(C_edges) - 1):
                mask = (
                    (df["T"]     >= T_edges[i]) & (df["T"]     < T_edges[i + 1]) &
                    (df["cloud"] >= C_edges[j]) & (df["cloud"] < C_edges[j + 1])
                )
                if mask.sum() < 30:
                    continue
                R_cell = df.loc[mask, "R"]
                T_vals.append(0.5 * (T_edges[i] + T_edges[i + 1]))
                C_vals.append(0.5 * (C_edges[j] + C_edges[j + 1]))
                R_vals.append(float(np.percentile(R_cell, 99)))

        self.grid_T    = np.array(T_vals)
        self.grid_C    = np.array(C_vals)
        self.grid_Rmax = np.array(R_vals)

        log(f"\nEnvelope grid: {len(T_vals)} cells with ≥ 30 observations")

    # ----------------------------------------------------------
    # 5.2  Parametric fit
    # ----------------------------------------------------------

    @staticmethod
    def _parametric_model(X: np.ndarray, a: float, b: float, alpha: float) -> np.ndarray:
        """
        R_max(T, C) = (a·C + b) · exp(−α·T) / (σ · T_K⁴)

        X : (N, 2) array with columns [T_celsius, C]
        """
        T   = X[:, 0]
        C   = X[:, 1]
        T_K = T + 273.15
        return (a * C + b) * np.exp(-alpha * T) / (SIGMA * T_K**4)

    def fit_parametric(self) -> None:
        """Fit (a, b, alpha) via scipy curve_fit."""
        X    = np.column_stack([self.grid_T, self.grid_C])
        popt, _ = curve_fit(
            self._parametric_model, X, self.grid_Rmax,
            p0=[0.5, 0.4, 0.03],
            bounds=([0, 0, 0], [5, 5, 1]),
            maxfev=10_000,
        )
        self.params_parametric = tuple(popt)
        pred = self._parametric_model(X, *popt)
        self.metrics["parametric"] = _regression_metrics(self.grid_Rmax, pred)

        a, b, alpha = popt
        log(f"\nParametric fit:  a = {a:.5f}   b = {b:.5f}   alpha = {alpha:.5f}")
        log(f"  R_max(0°C, C=0.5)  = {self.predict_parametric(np.array([[0.0, 0.5]]))[0]:.6f}")
        log(f"  R_max(20°C, C=0.5) = {self.predict_parametric(np.array([[20.0, 0.5]]))[0]:.6f}")
        m = self.metrics["parametric"]
        log(f"  RMSE={m['rmse']:.2e}   MAE={m['mae']:.2e}   R²={m['r2']:.4f}")

    def predict_parametric(self, X: np.ndarray) -> np.ndarray:
        """Predict R_max from (N,2) array [T_celsius, C]."""
        return self._parametric_model(X, *self.params_parametric)

    # ----------------------------------------------------------
    # 5.3  GP fit
    # ----------------------------------------------------------

    def fit_gp(self) -> None:
        """
        Fit a Gaussian Process regressor on the envelope grid (pure-NumPy RBF GP).

        Inputs are standardised before fitting.
        A random subsample of GP_MAX_TRAIN points is used if the grid is large
        (GP inference is O(n³) in the number of training points).
        """
        X = np.column_stack([self.grid_T, self.grid_C])
        y = self.grid_Rmax.copy()

        # Standardise inputs
        self._gp_X_mean  = X.mean(axis=0)
        self._gp_X_std   = X.std(axis=0) + 1e-10
        X_sc = (X - self._gp_X_mean) / self._gp_X_std

        # Subsample if needed
        if len(X_sc) > GP_MAX_TRAIN:
            rng   = np.random.default_rng(42)
            idx   = rng.choice(len(X_sc), GP_MAX_TRAIN, replace=False)
            X_tr, y_tr = X_sc[idx], y[idx]
        else:
            X_tr, y_tr = X_sc, y

        gp = _NumpyGP(noise=1e-4)
        gp.fit(X_tr, y_tr)
        self.gp_model = gp

        # Evaluate on full grid
        pred = gp.predict(X_sc)
        self.metrics["gp"] = _regression_metrics(y, pred)

        m = self.metrics["gp"]
        log(f"\nGP fit:  l={gp.l:.3f}   sigma={gp.sigma:.3f}   noise={gp.noise:.1e}")
        log(f"  RMSE={m['rmse']:.2e}   MAE={m['mae']:.2e}   R²={m['r2']:.4f}")

    def predict_gp(self, X: np.ndarray) -> np.ndarray:
        """Predict R_max from (N,2) array [T_celsius, C]."""
        X_sc = (X - self._gp_X_mean) / self._gp_X_std
        return self.gp_model.predict(X_sc)

    # ----------------------------------------------------------
    # 5.4  Fit comparison summary
    # ----------------------------------------------------------

    def print_comparison(self) -> None:
        log("\n── Envelope fit comparison ──────────────────────────────")
        log(f"{'Model':<14s}  {'RMSE':>10s}  {'MAE':>10s}  {'R²':>8s}")
        log(f"{'─'*14}  {'─'*10}  {'─'*10}  {'─'*8}")
        for name, m in self.metrics.items():
            log(f"{name:<14s}  {m['rmse']:>10.3e}  {m['mae']:>10.3e}  {m['r2']:>8.4f}")
        log("─" * 50)

    # ----------------------------------------------------------
    # 5.5  Relative position ρ
    # ----------------------------------------------------------

    def add_rho(self, df: pd.DataFrame, model: str = "parametric") -> pd.DataFrame:
        """
        Add column rho = R / R_max(T, C).

        model : "parametric" or "gp"

        rho ∈ [0, 1]  under the physical envelope.
        Values > 1 indicate reanalysis artefacts or outliers.
        """
        df   = df.copy()
        X    = df[["T", "cloud"]].values
        if model == "gp":
            Rmax = self.predict_gp(X)
        else:
            Rmax = self.predict_parametric(X)
        Rmax        = np.maximum(Rmax, 1e-30)
        df["R_max"] = Rmax
        df["rho"]   = df["R"].values / Rmax
        return df


def _regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    residuals = y_true - y_pred
    return {
        "rmse": float(np.sqrt(np.mean(residuals**2))),
        "mae":  float(np.mean(np.abs(residuals))),
        "r2":   _r2_score(y_true, y_pred),
    }


# ============================================================
# 6. REGIME CLASSIFICATION (ρ-based)
# ============================================================
#
# Boundaries are fractions of the local physical envelope R_max(T,C),
# not fixed absolute R thresholds.  This makes the classification
# invariant to the T-dependent radiative scale.
#
# Priority order (first matching condition wins):
#   1. clear        : cloud < CLOUD_CLEAR_THR
#   2. convective   : S > S_SAT  AND  rho > RHO_CONVECTIVE
#   3. stratiform   : S > S_SAT  AND  rho ≤ RHO_CONVECTIVE
#   4. mixed        : S ≤ S_SAT  AND  cloud > CLOUD_MIXED_THR
#   5. transitional : remainder
# ============================================================

REGIMES = ["clear", "convective", "stratiform", "mixed", "transitional"]

REGIME_COLORS = {
    "clear":        "orange",
    "convective":   "green",
    "stratiform":   "red",
    "mixed":        "steelblue",
    "transitional": "purple",
}


def classify_cloud_regime(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign cloud regime using ρ-based thresholds.

    Requires columns: cloud, S, rho.
    """
    df = df.copy()
    conditions = [
        df["cloud"] < CLOUD_CLEAR_THR,
        (df["S"] > S_SAT_THRESHOLD) & (df["rho"] > RHO_CONVECTIVE),
        (df["S"] > S_SAT_THRESHOLD) & (df["rho"] <= RHO_CONVECTIVE),
        (df["S"] <= S_SAT_THRESHOLD) & (df["cloud"] > CLOUD_MIXED_THR),
    ]
    labels = ["clear", "convective", "stratiform", "mixed"]
    df["cloud_regime"] = np.select(conditions, labels, default="transitional")
    return df


def print_regime_summary(df: pd.DataFrame) -> None:
    counts = df["cloud_regime"].value_counts()
    total  = len(df)
    log("\nCloud regime distribution:")
    for regime in REGIMES:
        n = counts.get(regime, 0)
        log(f"  {regime:<14s}: {n:>10,}  ({100*n/total:.1f} %)")


# ============================================================
# 7. PROBABILISTIC MODULE
# ============================================================
#
# 7.1  Conditional KDE: P(R | T_bin, C_bin)
#      For each (T,C) grid cell, fit a 1D KDE over the R values.
#      Stored as a 2D dict of scipy gaussian_kde objects.
#
# 7.2  Regime probability: P(regime | T, C)
#      Fraction of each regime in each (T,C) cell.
#      Output: dict of 2D arrays keyed by regime name.
#
# 7.3  Entropy map: H(T, C) = −Σ p·log₂(p)
#      Shannon entropy of regime distribution per cell.
#      High H → uncertain/transitional; Low H → dominated by one regime.
# ============================================================

class ProbabilisticModule:
    """
    Estimate P(regime | T, C) and H(T, C) on a 2D grid.

    Usage
    -----
    pm = ProbabilisticModule()
    pm.fit(df)          # build all probability maps
    pm.print_summary()
    p = pm.query(T=15.0, C=0.7)   # → dict of regime probs at (T=15, C=0.7)
    """

    def __init__(self, t_bins: int = T_BINS, c_bins: int = C_BINS):
        self.t_bins  = t_bins
        self.c_bins  = c_bins
        self.T_edges = None
        self.C_edges = None
        self.T_mid   = None
        self.C_mid   = None
        # P(regime | T, C): 2D arrays  shape (t_bins, c_bins)
        self.regime_prob: dict[str, np.ndarray] = {}
        # H(T, C): 2D array
        self.entropy: np.ndarray | None = None
        # count of samples per cell
        self.cell_count: np.ndarray | None = None

    # ----------------------------------------------------------
    # 7.2 / 7.3  Regime probability and entropy maps
    # ----------------------------------------------------------

    def fit(self, df: pd.DataFrame, min_count: int = 30) -> None:
        """
        Build P(regime | T, C) and entropy H(T, C) maps.

        Cells with fewer than min_count observations are set to NaN.
        """
        self.T_edges = np.linspace(df["T"].min(),     df["T"].max(),     self.t_bins + 1)
        self.C_edges = np.linspace(df["cloud"].min(), df["cloud"].max(), self.c_bins + 1)
        self.T_mid   = 0.5 * (self.T_edges[:-1] + self.T_edges[1:])
        self.C_mid   = 0.5 * (self.C_edges[:-1] + self.C_edges[1:])

        nt, nc = self.t_bins, self.c_bins
        self.cell_count = np.zeros((nt, nc), dtype=int)

        # Accumulate regime counts
        regime_counts: dict[str, np.ndarray] = {
            r: np.zeros((nt, nc), dtype=float) for r in REGIMES
        }

        # Vectorised bin assignment
        T_idx = np.searchsorted(self.T_edges[1:-1], df["T"].values)
        C_idx = np.searchsorted(self.C_edges[1:-1], df["cloud"].values)

        for r in REGIMES:
            mask = (df["cloud_regime"] == r).values
            np.add.at(regime_counts[r], (T_idx[mask], C_idx[mask]), 1)

        # Total count per cell
        total = np.zeros((nt, nc), dtype=float)
        for r in REGIMES:
            total += regime_counts[r]
        self.cell_count = total.astype(int)

        # Compute probabilities; NaN for sparse cells
        for r in REGIMES:
            p = np.full((nt, nc), np.nan)
            valid = total >= min_count
            p[valid] = regime_counts[r][valid] / total[valid]
            self.regime_prob[r] = p

        # Shannon entropy
        H = np.zeros((nt, nc))
        for r in REGIMES:
            p = self.regime_prob[r]
            valid = ~np.isnan(p) & (p > 0)
            H[valid] -= p[valid] * np.log2(p[valid])

        # NaN where total < min_count
        sparse = total < min_count
        H[sparse] = np.nan
        self.entropy = H

        n_valid = int(np.sum(~sparse))
        log(f"\nProbabilistic module:")
        log(f"  Grid: {nt} × {nc} cells | valid cells: {n_valid}")
        log(f"  Mean entropy (valid cells): {np.nanmean(H):.3f} bits")
        log(f"  Max  entropy (valid cells): {np.nanmax(H):.3f} bits")

    # ----------------------------------------------------------
    # Query: P(regime | T, C)
    # ----------------------------------------------------------

    def query(self, T: float, C: float) -> dict[str, float]:
        """
        Return P(regime | T, C) for a single point.

        If the point falls in a sparse or out-of-range cell, returns NaN
        for all regimes.
        """
        ti = int(np.clip(
            np.searchsorted(self.T_edges[1:-1], T),
            0, self.t_bins - 1,
        ))
        ci = int(np.clip(
            np.searchsorted(self.C_edges[1:-1], C),
            0, self.c_bins - 1,
        ))
        return {r: float(self.regime_prob[r][ti, ci]) for r in REGIMES}

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------

    def print_summary(self) -> None:
        log("\n── Probabilistic module summary ─────────────────────────")
        dominant = {}
        for r in REGIMES:
            p     = self.regime_prob[r]
            valid = ~np.isnan(p)
            dominant[r] = (
                f"mean={np.nanmean(p):.3f}  "
                f"max={np.nanmax(p):.3f}  "
                f"cells_dominant={int(np.sum(p[valid] > 0.5))}"
            )
            log(f"  {r:<14s}: {dominant[r]}")
        log("─" * 54)


# ============================================================
# 8. PLOTS
# ============================================================

def plot_surface_comparison(sf: SurfaceFit) -> None:
    """
    3-D surface plot: observed envelope grid vs parametric fit vs GP fit.
    """
    T_grid = sf.grid_T
    C_grid = sf.grid_C
    R_obs  = sf.grid_Rmax

    X = np.column_stack([T_grid, C_grid])
    R_par = sf.predict_parametric(X)
    R_gp  = sf.predict_gp(X)

    fig = plt.figure(figsize=(16, 6))
    fig.suptitle("Physical Envelope R_max(T, C) — Model Comparison", fontsize=13)

    for k, (title, R_vals, color) in enumerate([
        ("Observed (99th pct)",   R_obs,  "black"),
        ("Parametric fit",        R_par,  "crimson"),
        ("GP fit",                R_gp,   "steelblue"),
    ], 1):
        ax = fig.add_subplot(1, 3, k, projection="3d")
        sc = ax.scatter(T_grid, C_grid, R_vals,
                        c=R_vals, cmap="plasma", s=20, alpha=0.8)
        ax.set_xlabel("T (°C)", fontsize=9)
        ax.set_ylabel("Cloud Cover", fontsize=9)
        ax.set_zlabel("R_max", fontsize=9)
        ax.set_title(title, fontsize=10)
        fig.colorbar(sc, ax=ax, shrink=0.5, label="R_max")

    plt.tight_layout()
    save_figure(fig, "surface_comparison")
    plt.close(fig)


def plot_surface_residuals(sf: SurfaceFit) -> None:
    """
    2-D scatter: residuals (observed − predicted) for parametric and GP.
    """
    X     = np.column_stack([sf.grid_T, sf.grid_C])
    R_par = sf.predict_parametric(X)
    R_gp  = sf.predict_gp(X)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Envelope Fit Residuals  (observed − predicted)", fontsize=12)

    for ax, (label, R_pred) in zip(axes, [
        ("Parametric", R_par),
        ("GP",         R_gp),
    ]):
        resid = sf.grid_Rmax - R_pred
        vmax  = np.percentile(np.abs(resid), 95)
        sc = ax.scatter(sf.grid_T, sf.grid_C,
                        c=resid, cmap="RdBu", vmin=-vmax, vmax=vmax, s=30, alpha=0.9)
        plt.colorbar(sc, ax=ax, label="Residual")
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Cloud Cover")
        ax.set_title(f"{label}  |  R²={sf.metrics[label.lower()]['r2']:.4f}")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_figure(fig, "surface_residuals")
    plt.close(fig)


def plot_rho_scatter(df: pd.DataFrame) -> None:
    """
    2-D scatter in T–C space coloured by ρ = R / R_max(T,C).
    ρ = 1 means the state is exactly on the physical envelope.
    """
    sample = df.sample(min(150_000, len(df)), random_state=42)
    rho    = sample["rho"].clip(0, 1.05)

    fig, ax = plt.subplots(figsize=(10, 6))
    sc = ax.scatter(sample["T"], sample["cloud"],
                    c=rho, cmap="YlOrRd", s=2, alpha=0.3, vmin=0, vmax=1)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("ρ = R / R_max(T, C)")

    # Mark ρ = 1 contour via 2D histogram
    H, T_e, C_e = np.histogram2d(
        df["T"].values, df["cloud"].values,
        bins=[60, 30],
        weights=(df["rho"].values > 0.9).astype(float),
    )
    H_tot, _, _ = np.histogram2d(df["T"].values, df["cloud"].values, bins=[60, 30])
    with np.errstate(invalid="ignore"):
        frac = np.where(H_tot > 10, H / H_tot, np.nan)
    ax.contour(
        0.5*(T_e[:-1]+T_e[1:]),
        0.5*(C_e[:-1]+C_e[1:]),
        frac.T,
        levels=[0.05],
        colors=["red"], linewidths=1.5, linestyles="--",
    )
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Cloud Cover")
    ax.set_title("Relative Position on the Manifold  ρ = R / R_max(T, C)\n"
                 "dashed red: ρ = 0.9 contour (near-envelope states)")
    save_figure(fig, "rho_scatter")
    plt.close(fig)


def plot_R_3D_regimes(df: pd.DataFrame, n_sample: int = 80_000) -> None:
    """3-D regime scatter in R–T–Cloud space."""
    sample = df.sample(min(n_sample, len(df)), random_state=42)

    fig = plt.figure(figsize=(10, 8))
    ax  = fig.add_subplot(111, projection="3d")

    for regime, color in REGIME_COLORS.items():
        sub = sample[sample["cloud_regime"] == regime]
        if len(sub) == 0:
            continue
        ax.scatter(sub["T"], sub["cloud"], sub["R"],
                   s=2, alpha=0.35, color=color, label=regime)

    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Cloud Cover")
    ax.set_zlabel("R (Dynamics / Radiation)")
    ax.set_title("Cloud Regimes in R–T–Cloud Space  (ρ-based boundaries)")
    ax.legend(loc="upper left", fontsize=9, markerscale=4)
    save_figure(fig, "R_3D_regimes")
    plt.close(fig)


def plot_regime_probability_maps(pm: ProbabilisticModule) -> None:
    """
    One heatmap per regime showing P(regime | T, C).
    """
    TT, CC = np.meshgrid(pm.T_mid, pm.C_mid, indexing="ij")

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle("P(regime | T, C) — Conditional Regime Probability Maps", fontsize=13)

    regime_ax = {r: axes.flat[i] for i, r in enumerate(REGIMES)}

    for regime, ax in regime_ax.items():
        P   = pm.regime_prob[regime]  # (t_bins, c_bins)
        im  = ax.pcolormesh(pm.T_mid, pm.C_mid, P.T,
                            cmap="YlOrRd", vmin=0, vmax=1)
        plt.colorbar(im, ax=ax, label="P")
        ax.set_xlabel("Temperature (°C)", fontsize=9)
        ax.set_ylabel("Cloud Cover",      fontsize=9)
        ax.set_title(f"P({regime} | T, C)", fontsize=10,
                     color=REGIME_COLORS[regime])
        ax.grid(True, alpha=0.2)

    # Last panel: dominant regime map
    ax_dom = axes.flat[5]
    regime_idx = np.full(pm.entropy.shape, np.nan)
    for i, r in enumerate(REGIMES):
        p = pm.regime_prob[r]
        valid = ~np.isnan(p)
        # Mark dominant regime per cell
        # We'll build this by argmax over stacked arrays
    P_stack = np.stack(
        [np.nan_to_num(pm.regime_prob[r], nan=0.0) for r in REGIMES],
        axis=-1,
    )
    dom_idx = np.argmax(P_stack, axis=-1).astype(float)
    # NaN sparse cells
    sparse = pm.cell_count < 30
    dom_idx[sparse] = np.nan

    cmap_dom = matplotlib.colors.ListedColormap(
        [REGIME_COLORS[r] for r in REGIMES]
    )
    im = ax_dom.pcolormesh(pm.T_mid, pm.C_mid, dom_idx.T,
                           cmap=cmap_dom, vmin=-0.5, vmax=len(REGIMES) - 0.5)
    cbar = plt.colorbar(im, ax=ax_dom, ticks=range(len(REGIMES)))
    cbar.ax.set_yticklabels(REGIMES, fontsize=8)
    ax_dom.set_xlabel("Temperature (°C)", fontsize=9)
    ax_dom.set_ylabel("Cloud Cover",      fontsize=9)
    ax_dom.set_title("Dominant Regime", fontsize=10)
    ax_dom.grid(True, alpha=0.2)

    plt.tight_layout()
    save_figure(fig, "regime_probability_maps")
    plt.close(fig)


def plot_entropy_map(pm: ProbabilisticModule) -> None:
    """
    Shannon entropy H(T, C) = −Σ P(regime)·log₂P(regime).

    High H (→ max log₂5 ≈ 2.32 bits): all regimes equally probable.
    Low H (→ 0): one regime dominates.
    """
    fig, ax = plt.subplots(figsize=(9, 6))
    im = ax.pcolormesh(pm.T_mid, pm.C_mid, pm.entropy.T,
                       cmap="inferno_r", vmin=0, vmax=np.log2(len(REGIMES)))
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Shannon entropy H(T, C)  [bits]")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Cloud Cover")
    ax.set_title("Regime Uncertainty Map  H(T, C)\n"
                 "Dark = one regime dominates | Bright = high uncertainty")
    ax.grid(True, alpha=0.2)
    save_figure(fig, "entropy_map")
    plt.close(fig)


def plot_rho_distribution_by_regime(df: pd.DataFrame) -> None:
    """
    Distribution of ρ = R / R_max per cloud regime.
    Shows how each regime occupies the manifold.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    rho_max_plot = 1.1

    for regime, color in REGIME_COLORS.items():
        sub = df[df["cloud_regime"] == regime]["rho"].dropna()
        if len(sub) < 100:
            continue
        sub  = sub.clip(0, rho_max_plot)
        bins = np.linspace(0, rho_max_plot, 60)
        ax.hist(sub, bins=bins, alpha=0.5, density=True,
                color=color, label=f"{regime}  (n={len(sub):,})")

    ax.axvline(1.0, color="black", lw=1.5, ls="--", label="ρ = 1 (envelope)")
    ax.set_xlabel("ρ = R / R_max(T, C)")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of Relative Position ρ per Regime\n"
                 "ρ = 1 means the state lies exactly on the physical envelope")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    save_figure(fig, "rho_distribution_by_regime")
    plt.close(fig)


# ============================================================
# 8. ASYMMETRY ANALYSIS AT ρ = 1
# ============================================================
#
# Physical hypothesis: the ρ = 1 edge is an asymmetric attractor.
#   - Build-up (ρ → 1⁻): slow — moisture accumulation, cloud organisation
#   - Collapse  (ρ → 1⁺ then ρ ↓): fast — convective burst, precipitation
#
# Method:
#   Work on time-ordered ERA5 data grouped by grid cell (lat, lon).
#   For each cell, detect "edge events": consecutive steps where ρ
#   crosses above RHO_EDGE_THRESHOLD, then falls back below it.
#
#   Rise time  τ_rise : steps to go from RHO_BASE → RHO_EDGE_THRESHOLD
#   Fall time  τ_fall : steps to return from RHO_EDGE_THRESHOLD → RHO_BASE
#
#   Asymmetry index:
#     AI = (τ_fall − τ_rise) / (τ_fall + τ_rise)  ∈ [−1, +1]
#     AI > 0 → slow build, fast collapse   (physical expectation)
#     AI < 0 → fast build, slow collapse   (would be surprising)
#     AI = 0 → symmetric
# ============================================================

RHO_EDGE_THRESHOLD = 0.90   # "near-envelope" threshold
RHO_BASE           = 0.50   # "recovered" baseline
MIN_EVENT_STEPS    = 2      # minimum steps to count as a real event


def compute_rho_timeseries(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort DataFrame by (lat, lon, time) and compute ρ first-difference.

    Adds column:  drho_dt  (ρ change per ERA5 time step)
    Only consecutive rows from the same (lat, lon) cell are connected.
    """
    df = df.copy()

    # Round lat/lon to grid resolution to define spatial cells
    df["lat_r"] = df["lat"].round(1)
    df["lon_r"] = df["lon"].round(1)

    df = df.sort_values(["lat_r", "lon_r", "time"]).reset_index(drop=True)

    # Compute ρ difference only within the same spatial cell
    same_cell = (
        (df["lat_r"] == df["lat_r"].shift(1)) &
        (df["lon_r"] == df["lon_r"].shift(1))
    )
    df["drho_dt"] = np.where(same_cell, df["rho"].diff(), np.nan)

    return df


def detect_edge_events(
    df: pd.DataFrame,
    min_events: int = 100,
) -> dict:
    """
    Detect ρ build-up / collapse events in time-ordered data.

    For each spatial cell that has enough data, scan the ρ time series
    for complete cycles:  ρ < RHO_BASE  →  ρ > RHO_EDGE_THRESHOLD  →  ρ < RHO_BASE

    Returns a dict with:
      rise_times   : array of build-up durations [steps]
      fall_times   : array of collapse durations [steps]
      asymmetry_index : scalar AI
      n_events     : number of complete cycles found
    """
    rise_times: list[int] = []
    fall_times: list[int] = []

    cells = df.groupby(["lat_r", "lon_r"], sort=False)

    for _, cell_df in cells:
        rho = cell_df["rho"].values
        n   = len(rho)
        if n < 20:
            continue

        i = 0
        while i < n:
            # Find start: ρ below base
            if rho[i] >= RHO_BASE:
                i += 1
                continue

            # Scan upward: how many steps to reach threshold?
            j = i + 1
            while j < n and rho[j] < RHO_EDGE_THRESHOLD:
                j += 1
            if j >= n:
                break
            rise = j - i
            if rise < MIN_EVENT_STEPS:
                i = j
                continue

            # Scan downward: how many steps to return to base?
            k = j + 1
            while k < n and rho[k] >= RHO_BASE:
                k += 1
            fall = k - j
            if fall < MIN_EVENT_STEPS or k >= n:
                i = k
                continue

            rise_times.append(rise)
            fall_times.append(fall)
            i = k

    rise_arr = np.array(rise_times, dtype=float)
    fall_arr = np.array(fall_times, dtype=float)
    n_events = len(rise_arr)

    if n_events < min_events:
        log(f"\nEdge events: only {n_events} complete cycles found "
            f"(need {min_events}). Asymmetry index not computed.")
        return {
            "rise_times": rise_arr,
            "fall_times": fall_arr,
            "asymmetry_index": np.nan,
            "n_events": n_events,
        }

    tau_rise = float(np.median(rise_arr))
    tau_fall = float(np.median(fall_arr))
    ai       = (tau_fall - tau_rise) / (tau_fall + tau_rise + 1e-10)

    log(f"\nEdge event analysis:")
    log(f"  Complete cycles detected : {n_events:,}")
    log(f"  Median rise time  τ_rise : {tau_rise:.2f} steps")
    log(f"  Median fall time  τ_fall : {tau_fall:.2f} steps")
    log(f"  Asymmetry index   AI     : {ai:+.4f}  "
        f"({'slow build / fast collapse' if ai > 0.05 else 'fast build / slow collapse' if ai < -0.05 else 'approximately symmetric'})")

    return {
        "rise_times": rise_arr,
        "fall_times": fall_arr,
        "asymmetry_index": ai,
        "n_events": n_events,
    }


def plot_asymmetry(events: dict) -> None:
    """
    Side-by-side histogram of rise vs fall times + asymmetry summary.
    """
    if events["n_events"] < 10:
        log("  Asymmetry plot: insufficient events, skipping")
        return

    rise = events["rise_times"]
    fall = events["fall_times"]
    ai   = events["asymmetry_index"]

    clip = int(np.percentile(np.concatenate([rise, fall]), 97))
    bins = np.arange(0, clip + 2, 1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        f"Asymmetry at the ρ = 1 Edge  (AI = {ai:+.4f})\n"
        f"AI > 0 → slow build-up, fast collapse  |  "
        f"n = {events['n_events']:,} events",
        fontsize=12,
    )

    # Rise time distribution
    axes[0].hist(rise.clip(0, clip), bins=bins,
                 color="steelblue", alpha=0.8, density=True)
    axes[0].axvline(np.median(rise), color="navy", lw=2, ls="--",
                    label=f"median = {np.median(rise):.1f}")
    axes[0].set_xlabel("Rise time [ERA5 steps]")
    axes[0].set_ylabel("Density")
    axes[0].set_title("Build-up duration  (ρ → 1)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Fall time distribution
    axes[1].hist(fall.clip(0, clip), bins=bins,
                 color="crimson", alpha=0.8, density=True)
    axes[1].axvline(np.median(fall), color="darkred", lw=2, ls="--",
                    label=f"median = {np.median(fall):.1f}")
    axes[1].set_xlabel("Fall time [ERA5 steps]")
    axes[1].set_ylabel("Density")
    axes[1].set_title("Collapse duration  (ρ → base)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    # Overlay comparison
    axes[2].hist(rise.clip(0, clip), bins=bins,
                 color="steelblue", alpha=0.5, density=True, label="Rise")
    axes[2].hist(fall.clip(0, clip), bins=bins,
                 color="crimson",   alpha=0.5, density=True, label="Fall")
    axes[2].axvline(np.median(rise), color="steelblue", lw=2, ls="--")
    axes[2].axvline(np.median(fall), color="crimson",   lw=2, ls="--")
    axes[2].set_xlabel("Duration [ERA5 steps]")
    axes[2].set_ylabel("Density")
    axes[2].set_title("Rise vs Fall — Overlap")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    save_figure(fig, "edge_asymmetry")
    plt.close(fig)


# ============================================================
# 9. CONDITIONAL DRIFT FIELD  dρ/dt | ρ, T, C
# ============================================================
#
# For each point in (ρ, T) space, compute the mean observed ρ-derivative.
# This gives a vector field on the manifold:
#   - Where drift > 0: system is pushed toward the envelope (ρ → 1)
#   - Where drift < 0: system relaxes away from the envelope
#   - The ρ = 1 line should act as an attractor from below and
#     a repeller from above (consistent with asymmetry hypothesis)
#
# The drift field is the empirical analogue of the Fokker–Planck
# drift term for the ρ coordinate on the manifold.
# ============================================================

RHO_DRIFT_BINS = 20    # bins along ρ axis for drift computation
T_DRIFT_BINS   = 30    # bins along T axis for drift computation


def compute_drift_field(df: pd.DataFrame) -> dict:
    """
    Compute binned mean drift  E[dρ/dt | ρ_bin, T_bin].

    Returns a dict with:
      rho_mid    : (RHO_DRIFT_BINS,)  bin centres along ρ
      T_mid      : (T_DRIFT_BINS,)    bin centres along T
      drift_mean : (RHO_DRIFT_BINS, T_DRIFT_BINS)  mean dρ/dt
      drift_std  : (RHO_DRIFT_BINS, T_DRIFT_BINS)  std of dρ/dt
      n_count    : (RHO_DRIFT_BINS, T_DRIFT_BINS)  sample count
    """
    valid = df.dropna(subset=["drho_dt", "rho", "T"])
    # Focus on physically meaningful ρ range
    valid = valid[(valid["rho"] >= 0) & (valid["rho"] <= 1.2)]

    rho_edges = np.linspace(0.0, 1.2, RHO_DRIFT_BINS + 1)
    T_edges   = np.linspace(valid["T"].min(), valid["T"].max(), T_DRIFT_BINS + 1)
    rho_mid   = 0.5 * (rho_edges[:-1] + rho_edges[1:])
    T_mid     = 0.5 * (T_edges[:-1]   + T_edges[1:])

    shape = (RHO_DRIFT_BINS, T_DRIFT_BINS)
    drift_sum = np.zeros(shape)
    drift_sq  = np.zeros(shape)
    n_count   = np.zeros(shape, dtype=int)

    rho_idx = np.searchsorted(rho_edges[1:-1], valid["rho"].values)
    T_idx   = np.searchsorted(T_edges[1:-1],   valid["T"].values)
    drho    = valid["drho_dt"].values

    for ri, ti, d in zip(rho_idx, T_idx, drho):
        if 0 <= ri < RHO_DRIFT_BINS and 0 <= ti < T_DRIFT_BINS:
            drift_sum[ri, ti] += d
            drift_sq[ri, ti]  += d * d
            n_count[ri, ti]   += 1

    with np.errstate(invalid="ignore"):
        drift_mean = np.where(n_count > 5,
                              drift_sum / n_count, np.nan)
        variance   = np.where(n_count > 5,
                              drift_sq / n_count - drift_mean**2, np.nan)
        drift_std  = np.sqrt(np.maximum(variance, 0))

    n_valid = int(np.sum(n_count > 5))
    log(f"\nDrift field:")
    log(f"  Grid: {RHO_DRIFT_BINS} × {T_DRIFT_BINS} | valid cells: {n_valid}")

    # Find zero-crossing line (where drift changes sign along ρ)
    for ti in range(T_DRIFT_BINS):
        col = drift_mean[:, ti]
        sign_changes = np.where(np.diff(np.sign(col[~np.isnan(col)])))[0]
        if len(sign_changes) > 0:
            rho_zero = rho_mid[~np.isnan(col)][sign_changes[0]]
        else:
            rho_zero = np.nan

    # Summary
    pos_frac = float(np.nanmean(drift_mean > 0))
    log(f"  Fraction of cells with positive drift (push toward ρ=1): {pos_frac:.3f}")

    return {
        "rho_mid":    rho_mid,
        "T_mid":      T_mid,
        "drift_mean": drift_mean,
        "drift_std":  drift_std,
        "n_count":    n_count,
    }


def plot_drift_field(drift: dict) -> None:
    """
    Three-panel plot of the conditional drift field.

    Panel 1: 2D heatmap of mean dρ/dt in (ρ, T) space
    Panel 2: Mean drift as function of ρ only (T-averaged)
    Panel 3: Drift std (noise level) in (ρ, T) space
    """
    dm   = drift["drift_mean"]
    ds   = drift["drift_std"]
    rmid = drift["rho_mid"]
    tmid = drift["T_mid"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    fig.suptitle(
        "Conditional Drift Field  E[dρ/dt | ρ, T]\n"
        "Positive drift → system moves toward ρ = 1  |  "
        "Negative drift → system relaxes",
        fontsize=12,
    )

    # Panel 1: 2D drift map
    vmax = np.nanpercentile(np.abs(dm), 95)
    im1  = axes[0].pcolormesh(tmid, rmid, dm,
                               cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axes[0].axhline(RHO_EDGE_THRESHOLD, color="black", lw=1.5, ls="--",
                    label=f"ρ = {RHO_EDGE_THRESHOLD}")
    axes[0].axhline(1.0, color="red", lw=2, ls="-", label="ρ = 1 (envelope)")
    plt.colorbar(im1, ax=axes[0], label="Mean dρ/dt")
    axes[0].set_xlabel("Temperature (°C)")
    axes[0].set_ylabel("ρ = R / R_max")
    axes[0].set_title("Mean Drift  E[dρ/dt | ρ, T]")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.2)

    # Panel 2: T-averaged drift vs ρ
    mean_over_T = np.nanmean(dm, axis=1)
    std_over_T  = np.nanstd(dm,  axis=1)
    axes[1].fill_between(rmid,
                         mean_over_T - std_over_T,
                         mean_over_T + std_over_T,
                         alpha=0.25, color="steelblue")
    axes[1].plot(rmid, mean_over_T, color="steelblue", lw=2.5)
    axes[1].axhline(0, color="black", lw=1, ls="--")
    axes[1].axvline(RHO_EDGE_THRESHOLD, color="gray",  lw=1.5, ls="--",
                    label=f"ρ = {RHO_EDGE_THRESHOLD}")
    axes[1].axvline(1.0, color="red", lw=2, label="ρ = 1")
    axes[1].set_xlabel("ρ = R / R_max")
    axes[1].set_ylabel("Mean dρ/dt  (T-averaged)")
    axes[1].set_title("Drift Profile along ρ")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)

    # Mark zero-crossing
    sign_change = np.where(np.diff(np.sign(mean_over_T[~np.isnan(mean_over_T)])))[0]
    if len(sign_change):
        rho_eq = rmid[~np.isnan(mean_over_T)][sign_change[0]]
        axes[1].axvline(rho_eq, color="green", lw=2, ls=":",
                        label=f"drift = 0 at ρ ≈ {rho_eq:.2f}")
        axes[1].legend(fontsize=8)
        log(f"  Drift zero-crossing (equilibrium): ρ ≈ {rho_eq:.3f}")

    # Panel 3: drift std map
    im3 = axes[2].pcolormesh(tmid, rmid, ds, cmap="YlOrRd")
    axes[2].axhline(1.0, color="red", lw=2, ls="-")
    plt.colorbar(im3, ax=axes[2], label="Std of dρ/dt")
    axes[2].set_xlabel("Temperature (°C)")
    axes[2].set_ylabel("ρ = R / R_max")
    axes[2].set_title("Drift Noise Level  std[dρ/dt | ρ, T]")
    axes[2].grid(True, alpha=0.2)

    plt.tight_layout()
    save_figure(fig, "drift_field")
    plt.close(fig)


# ============================================================
# 10. PROBABILISTIC STATE ESTIMATOR
# ============================================================
#
# This is the weather probability estimator module.
#
# IMPORTANT SCOPE NOTE:
#   This module provides probability estimates, NOT forecasts.
#   It does not predict specific weather events or outcomes.
#   It answers: "Given the current manifold position (ρ, T, C),
#   what is the probability distribution over regimes k steps ahead,
#   based on observed historical transition statistics?"
#
# Method:
#   10.1  ρ-transition kernel:
#         For each (ρ_bin, T_bin) cell, compute the empirical distribution
#         of ρ at the next time step.  This gives P(ρ_{t+1} | ρ_t, T_t).
#
#   10.2  Regime transition matrix:
#         P(regime_{t+1} | regime_t) — marginalised over (T, C).
#         Shows which regimes are "sticky" and which are transient.
#
#   10.3  query_state(rho, T, C, k=1):
#         Multi-step probability estimate via k applications of the
#         transition kernel.  Returns a regime probability vector.
# ============================================================

RHO_TRANS_BINS = 15    # ρ bins for transition kernel
T_TRANS_BINS   = 20    # T bins for transition kernel


class StateEstimator:
    """
    Probabilistic state estimator for the cloud manifold.

    Builds empirical transition kernels from consecutive ERA5 time steps
    and provides multi-step regime probability estimates.

    Usage
    -----
    se = StateEstimator()
    se.fit(df)
    probs = se.query_state(rho=0.85, T=15.0, C=0.7, k=1)
    se.print_summary()
    """

    def __init__(
        self,
        rho_bins: int = RHO_TRANS_BINS,
        t_bins:   int = T_TRANS_BINS,
    ):
        self.rho_bins = rho_bins
        self.t_bins   = t_bins

        self.rho_edges: np.ndarray | None = None
        self.T_edges:   np.ndarray | None = None
        self.rho_mid:   np.ndarray | None = None
        self.T_mid:     np.ndarray | None = None

        # P(ρ_{t+1} | ρ_t, T_t): shape (rho_bins, t_bins, rho_bins)
        # transition_kernel[i, j, :] = probability distribution over next ρ bin
        self.transition_kernel: np.ndarray | None = None

        # P(regime_{t+1} | regime_t): shape (n_regimes, n_regimes)
        self.regime_transition: np.ndarray | None = None

        # Marginal regime counts per (ρ_bin, T_bin) cell
        # shape: (rho_bins, t_bins, n_regimes)
        self.cell_regime_prob: np.ndarray | None = None

    # ----------------------------------------------------------
    # 10.1  ρ-transition kernel
    # ----------------------------------------------------------

    def fit(self, df: pd.DataFrame, min_count: int = 20) -> None:
        """
        Build transition kernel and regime transition matrix from data.

        Requires columns: rho, drho_dt, T, cloud_regime, time.
        Consecutive rows must be from the same spatial cell (use after
        compute_rho_timeseries).
        """
        valid = df.dropna(subset=["rho", "drho_dt", "T", "cloud_regime"])
        valid = valid[(valid["rho"] >= 0) & (valid["rho"] <= 1.2)]

        self.rho_edges = np.linspace(0.0, 1.2, self.rho_bins + 1)
        self.T_edges   = np.linspace(valid["T"].min(),
                                     valid["T"].max(), self.t_bins + 1)
        self.rho_mid   = 0.5 * (self.rho_edges[:-1] + self.rho_edges[1:])
        self.T_mid     = 0.5 * (self.T_edges[:-1]   + self.T_edges[1:])

        nr, nt = self.rho_bins, self.t_bins
        n_reg  = len(REGIMES)

        # Transition kernel: counts[i, j, k] = transitions from (ρ_bin=i, T_bin=j) → ρ_bin=k
        counts = np.zeros((nr, nt, nr), dtype=float)

        # Regime transition counts: reg_counts[i, j] = transitions from regime i → regime j
        reg_counts = np.zeros((n_reg, n_reg), dtype=float)

        # Cell regime distribution: shape (nr, nt, n_reg)
        cell_reg = np.zeros((nr, nt, n_reg), dtype=float)

        rho_t  = valid["rho"].values
        T_t    = valid["T"].values
        drho   = valid["drho_dt"].values
        rho_t1 = rho_t + drho   # ρ at next step

        reg_t  = valid["cloud_regime"].values
        reg_map = {r: i for i, r in enumerate(REGIMES)}

        ri_t  = np.searchsorted(self.rho_edges[1:-1], rho_t)
        ri_t1 = np.searchsorted(self.rho_edges[1:-1],
                                 np.clip(rho_t1, 0, 1.2))
        ti    = np.searchsorted(self.T_edges[1:-1], T_t)

        for i in range(len(rho_t)):
            r0, r1, tj = ri_t[i], ri_t1[i], ti[i]
            if not (0 <= r0 < nr and 0 <= r1 < nr and 0 <= tj < nt):
                continue
            counts[r0, tj, r1] += 1

            # Regime transition
            reg_i = reg_map.get(reg_t[i], -1)
            if i + 1 < len(reg_t) and reg_i >= 0:
                reg_j = reg_map.get(reg_t[i + 1] if i + 1 < len(reg_t)
                                    else reg_t[i], -1)
                # Use next-step regime from the data directly
                # (approximate: use drho sign to estimate)
                if reg_j >= 0:
                    reg_counts[reg_i, reg_j] += 1

            cell_reg[r0, tj, reg_i] += 1

        # Normalise to probabilities
        row_sums = counts.sum(axis=2, keepdims=True)
        self.transition_kernel = np.where(
            row_sums > min_count,
            counts / (row_sums + 1e-30),
            np.nan,
        )

        reg_row_sums = reg_counts.sum(axis=1, keepdims=True)
        self.regime_transition = np.where(
            reg_row_sums > 0,
            reg_counts / (reg_row_sums + 1e-30),
            1.0 / n_reg,
        )

        cell_reg_sums = cell_reg.sum(axis=2, keepdims=True)
        self.cell_regime_prob = np.where(
            cell_reg_sums > min_count,
            cell_reg / (cell_reg_sums + 1e-30),
            np.nan,
        )

        n_valid_cells = int(np.sum(
            ~np.isnan(self.transition_kernel[:, :, 0])
        ))
        log(f"\nState estimator:")
        log(f"  Transition kernel: {nr} × {nt} × {nr} | "
            f"valid (ρ, T) cells: {n_valid_cells}")
        log(f"  Regime transition matrix built from "
            f"{int(reg_counts.sum()):,} observed transitions")

    # ----------------------------------------------------------
    # 10.2  Regime transition matrix print
    # ----------------------------------------------------------

    def print_regime_transition(self) -> None:
        if self.regime_transition is None:
            return
        log("\nRegime transition matrix  P(regime_{t+1} | regime_t):")
        header = f"  {'':14s}" + "".join(f"{r[:6]:>8s}" for r in REGIMES)
        log(header)
        log("  " + "─" * (14 + 8 * len(REGIMES)))
        for i, r_from in enumerate(REGIMES):
            row = "  " + f"{r_from:<14s}"
            for j in range(len(REGIMES)):
                p = self.regime_transition[i, j]
                row += f"{p:8.3f}"
            log(row)

    # ----------------------------------------------------------
    # 10.3  Multi-step query
    # ----------------------------------------------------------

    def query_state(
        self,
        rho: float,
        T:   float,
        C:   float,
        k:   int = 1,
    ) -> dict[str, float]:
        """
        Estimate P(regime_{t+k} | ρ_t=rho, T_t=T, C_t=C).

        NOTE: This is a probability estimate, NOT a weather forecast.
        It reflects historical transition statistics from ERA5 data.
        k=1 means one ERA5 time step ahead.

        Parameters
        ----------
        rho : current ρ value
        T   : current temperature [°C]
        C   : current cloud cover [0,1]
        k   : number of steps ahead

        Returns
        -------
        dict mapping regime name → probability
        """
        if self.transition_kernel is None:
            return {r: np.nan for r in REGIMES}

        ri = int(np.clip(
            np.searchsorted(self.rho_edges[1:-1], rho),
            0, self.rho_bins - 1,
        ))
        ti = int(np.clip(
            np.searchsorted(self.T_edges[1:-1], T),
            0, self.t_bins - 1,
        ))

        # Get ρ distribution after k steps by repeated application
        # of the transition kernel (marginalised over T, treating T as
        # approximately stationary over short horizons)
        rho_prob = np.zeros(self.rho_bins)
        rho_prob[ri] = 1.0

        for _ in range(k):
            new_prob = np.zeros(self.rho_bins)
            for ri2 in range(self.rho_bins):
                if rho_prob[ri2] == 0:
                    continue
                kernel_row = self.transition_kernel[ri2, ti, :]
                if np.any(~np.isnan(kernel_row)):
                    kernel_row = np.nan_to_num(kernel_row, nan=0.0)
                    s = kernel_row.sum()
                    if s > 0:
                        new_prob += rho_prob[ri2] * kernel_row / s
                    else:
                        new_prob[ri2] += rho_prob[ri2]
                else:
                    new_prob[ri2] += rho_prob[ri2]
            rho_prob = new_prob

        # Convert ρ distribution to regime probabilities
        # using cell_regime_prob[ρ_bin, T_bin, regime]
        regime_probs = np.zeros(len(REGIMES))
        for ri2 in range(self.rho_bins):
            if rho_prob[ri2] == 0:
                continue
            cell_p = self.cell_regime_prob[ri2, ti, :]
            if np.any(~np.isnan(cell_p)):
                cell_p = np.nan_to_num(cell_p, nan=0.0)
                regime_probs += rho_prob[ri2] * cell_p

        # Normalise
        s = regime_probs.sum()
        if s > 0:
            regime_probs /= s
        else:
            regime_probs[:] = 1.0 / len(REGIMES)

        return {r: float(p) for r, p in zip(REGIMES, regime_probs)}

    def print_summary(self) -> None:
        self.print_regime_transition()


def plot_transition_kernel(se: StateEstimator) -> None:
    """
    Visualise the ρ-transition kernel as a 2D heatmap for selected T bins.

    Each row of the kernel  P(ρ_{t+1} | ρ_t, T)  is shown as a column
    in the heatmap.  The diagonal would mean "no change"; off-diagonal
    structure reveals systematic drift.
    """
    if se.transition_kernel is None:
        return

    # Pick three representative T slices: cold, mid, warm
    nt    = se.t_bins
    t_idx = [nt // 6, nt // 2, 5 * nt // 6]
    t_labels = [f"T ≈ {se.T_mid[i]:.0f}°C" for i in t_idx]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        "ρ-Transition Kernel  P(ρ_{t+1} | ρ_t, T)\n"
        "Columns = initial ρ bin  |  Rows = next-step ρ bin\n"
        "Diagonal = no change; below diagonal = collapse; above = build-up",
        fontsize=11,
    )

    for ax, ti, label in zip(axes, t_idx, t_labels):
        K = se.transition_kernel[:, ti, :].T   # (next_ρ, current_ρ)
        K_plot = np.nan_to_num(K, nan=0.0)
        im = ax.pcolormesh(se.rho_mid, se.rho_mid, K_plot,
                           cmap="YlOrRd", vmin=0, vmax=1)
        ax.plot([0, 1.2], [0, 1.2], color="white", lw=1, ls="--", alpha=0.6)
        ax.axhline(1.0, color="red", lw=1.5, ls="--", alpha=0.8)
        ax.axvline(1.0, color="red", lw=1.5, ls="--", alpha=0.8)
        plt.colorbar(im, ax=ax, label="P")
        ax.set_xlabel("Current ρ (at t)")
        ax.set_ylabel("Next ρ (at t+1)")
        ax.set_title(label)

    plt.tight_layout()
    save_figure(fig, "transition_kernel")
    plt.close(fig)


def plot_regime_transition_matrix(se: StateEstimator) -> None:
    """Heatmap of the regime-to-regime transition matrix."""
    if se.regime_transition is None:
        return

    M = se.regime_transition
    n = len(REGIMES)

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(M, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
    plt.colorbar(im, ax=ax, label="Transition probability")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(REGIMES, rotation=30, ha="right")
    ax.set_yticklabels(REGIMES)
    ax.set_xlabel("Regime at t+1")
    ax.set_ylabel("Regime at t")
    ax.set_title("Regime Transition Matrix  P(regime_{t+1} | regime_t)\n"
                 "Diagonal = self-persistence  |  Off-diagonal = transitions")

    # Annotate cells
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{M[i,j]:.2f}",
                    ha="center", va="center", fontsize=9,
                    color="white" if M[i, j] > 0.5 else "black")

    plt.tight_layout()
    save_figure(fig, "regime_transition_matrix")
    plt.close(fig)


def plot_state_estimator_example(
    se:     StateEstimator,
    pm:     "ProbabilisticModule",
    states: list[tuple[float, float, float]] | None = None,
) -> None:
    """
    Show P(regime | ρ, T, C) for several example states at k = 1, 3, 6 steps.

    states : list of (rho, T, C) tuples to query.
             If None, uses four representative states.
    """
    if states is None:
        states = [
            (0.3,  10.0, 0.3),   # low ρ, mild, low cloud
            (0.7,   5.0, 0.7),   # mid ρ, cool, high cloud
            (0.92, -5.0, 0.9),   # near-envelope, cold
            (0.98, 15.0, 0.85),  # at edge, warm
        ]

    steps    = [1, 3, 6]
    n_states = len(states)
    n_steps  = len(steps)

    fig, axes = plt.subplots(
        n_states, n_steps,
        figsize=(5 * n_steps, 3.5 * n_states),
        sharey="row",
    )
    fig.suptitle(
        "Probabilistic State Estimator  P(regime_{t+k} | ρ_t, T_t, C_t)\n"
        "NOTE: probability estimate from historical ERA5 statistics — not a forecast",
        fontsize=12,
    )

    colors = [REGIME_COLORS[r] for r in REGIMES]

    for si, (rho, T, C) in enumerate(states):
        for ki, k in enumerate(steps):
            ax    = axes[si][ki] if n_states > 1 else axes[ki]
            probs = se.query_state(rho=rho, T=T, C=C, k=k)
            vals  = [probs[r] for r in REGIMES]

            bars = ax.bar(REGIMES, vals, color=colors, alpha=0.85, edgecolor="white")
            ax.set_ylim(0, 1)
            ax.set_ylabel("P" if ki == 0 else "")
            ax.set_title(
                f"ρ={rho:.2f}  T={T:.0f}°C  C={C:.1f}\n"
                f"k = {k} step{'s' if k > 1 else ''}",
                fontsize=9,
            )
            ax.tick_params(axis="x", rotation=35, labelsize=8)
            ax.grid(True, axis="y", alpha=0.3)

            # Annotate dominant regime
            dom_idx = int(np.argmax(vals))
            ax.text(
                dom_idx, vals[dom_idx] + 0.02,
                f"{vals[dom_idx]:.2f}",
                ha="center", va="bottom", fontsize=8, fontweight="bold",
            )

    plt.tight_layout()
    save_figure(fig, "state_estimator_examples")
    plt.close(fig)


# ============================================================
# 9. MARKDOWN EXPORT
# ============================================================

def export_markdown(
    sf: SurfaceFit,
    pm: "ProbabilisticModule",
    events: dict | None = None,
    se: "StateEstimator | None" = None,
    version: str = "v3.0",
) -> None:
    """Write full log + fit summary + asymmetry + estimator to Markdown."""
    path = os.path.join(OUTPUT_DIR, f"cloud_manifold_{version}_report.md")

    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# Cloud Manifold Pipeline — Run Report\n\n")
        f.write(f"**Version:** {version}  \n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  \n\n")

        f.write("---\n\n## Surface Fit\n\n")
        f.write("| Model | RMSE | MAE | R² |\n")
        f.write("|-------|------|-----|----|\n")
        for name, m in sf.metrics.items():
            f.write(f"| {name} | {m['rmse']:.3e} | {m['mae']:.3e} | {m['r2']:.4f} |\n")
        if sf.params_parametric:
            a, b, alpha = sf.params_parametric
            f.write(f"\n**Parametric model:**  "
                    f"R_max(T,C) = ({a:.4f}·C + {b:.4f})·exp(−{alpha:.4f}·T) / σ·T_K⁴\n\n")

        f.write("---\n\n## Probabilistic Module\n\n")
        f.write(f"Grid: {pm.t_bins} × {pm.c_bins}  |  "
                f"Valid cells: {int(np.sum(pm.cell_count >= 30))}  |  "
                f"Mean H: {np.nanmean(pm.entropy):.3f} bits\n\n")

        if events is not None and not np.isnan(events.get("asymmetry_index", np.nan)):
            ai = events["asymmetry_index"]
            f.write("---\n\n## Edge Asymmetry\n\n")
            f.write(f"| Metric | Value |\n|--------|-------|\n")
            f.write(f"| Events detected | {events['n_events']:,} |\n")
            f.write(f"| Median rise time | {np.median(events['rise_times']):.2f} steps |\n")
            f.write(f"| Median fall time | {np.median(events['fall_times']):.2f} steps |\n")
            f.write(f"| Asymmetry index AI | {ai:+.4f} |\n")
            direction = ("slow build / fast collapse" if ai > 0.05
                         else "fast build / slow collapse" if ai < -0.05
                         else "approximately symmetric")
            f.write(f"| Interpretation | {direction} |\n\n")

        f.write("---\n\n## Full Log\n\n```\n")
        f.write("\n".join(LOG_BUFFER))
        f.write("\n```\n")

    log(f"\nMarkdown report: {path}")


# ============================================================
# 11. MAIN
# ============================================================

def main() -> None:
    log("=" * 65)
    log("CLOUD MANIFOLD PIPELINE  v3.0")
    log(f"Run started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 65)

    # ── Load ────────────────────────────────────────────────────
    df = load_era5(DATA_PATH, SAMPLE_STEP, force=FORCE_REPROCESS)

    # ── Physics ─────────────────────────────────────────────────
    df = compute_R(df)

    # ── Surface fit ─────────────────────────────────────────────
    log("\n[Section 5] Surface fit module")
    sf = SurfaceFit()
    sf.build_envelope_grid(df)
    sf.fit_parametric()
    sf.fit_gp()
    sf.print_comparison()
    df = sf.add_rho(df, model="parametric")

    # ── Regime classification ────────────────────────────────────
    log("\n[Section 6] Regime classification (ρ-based)")
    df = classify_cloud_regime(df)
    print_regime_summary(df)

    # ── Probabilistic module ─────────────────────────────────────
    log("\n[Section 7] Probabilistic module")
    pm = ProbabilisticModule(t_bins=T_BINS, c_bins=C_BINS)
    pm.fit(df)
    pm.print_summary()

    T_ex, C_ex = 15.0, 0.7
    q = pm.query(T=T_ex, C=C_ex)
    log(f"\nExample query  P(regime | T={T_ex}°C, C={C_ex}):")
    for r, p in q.items():
        bar = "█" * int(round(p * 20)) if not np.isnan(p) else "—"
        log(f"  {r:<14s}: {p:.3f}  {bar}")

    # ── ρ time series (needed for Sections 8–10) ─────────────────
    log("\n[Section 8–10 prep] Building ρ time series ...")
    df = compute_rho_timeseries(df)

    # ── Section 8: Asymmetry analysis ────────────────────────────
    log("\n[Section 8] Asymmetry analysis at ρ = 1")
    events = detect_edge_events(df)
    plot_asymmetry(events)

    # ── Section 9: Drift field ────────────────────────────────────
    log("\n[Section 9] Conditional drift field")
    drift = compute_drift_field(df)
    plot_drift_field(drift)

    # ── Section 10: State estimator ───────────────────────────────
    log("\n[Section 10] Probabilistic state estimator")
    se = StateEstimator(rho_bins=RHO_TRANS_BINS, t_bins=T_TRANS_BINS)
    se.fit(df)
    se.print_summary()

    # Example multi-step queries
    log(f"\nState estimator queries:")
    example_states = [
        (0.3,  10.0, 0.3,  "low ρ, mild"),
        (0.7,   5.0, 0.7,  "mid ρ, cool"),
        (0.92, -5.0, 0.9,  "near-envelope, cold"),
        (0.98, 15.0, 0.85, "at edge, warm"),
    ]
    for rho, T, C, label in example_states:
        log(f"\n  State: {label}  (ρ={rho}, T={T}°C, C={C})")
        for k in [1, 3, 6]:
            probs = se.query_state(rho=rho, T=T, C=C, k=k)
            dom   = max(probs, key=lambda x: probs[x])
            log(f"    k={k}: dominant={dom:<14s} "
                + "  ".join(f"{r[:4]}={p:.2f}" for r, p in probs.items()))

    # ── Plots ────────────────────────────────────────────────────
    log("\n[Section 11] Generating all plots ...")
    plot_surface_comparison(sf)
    plot_surface_residuals(sf)
    plot_rho_scatter(df)
    plot_R_3D_regimes(df)
    plot_regime_probability_maps(pm)
    plot_entropy_map(pm)
    plot_rho_distribution_by_regime(df)
    plot_transition_kernel(se)
    plot_regime_transition_matrix(se)
    plot_state_estimator_example(se, pm)

    # ── Export ───────────────────────────────────────────────────
    export_markdown(sf, pm, events=events, se=se, version="v3.0")

    log("\n" + "=" * 65)
    log(f"Done. All outputs in: {os.path.abspath(OUTPUT_DIR)}/")
    log("=" * 65)


if __name__ == "__main__":
    main()
