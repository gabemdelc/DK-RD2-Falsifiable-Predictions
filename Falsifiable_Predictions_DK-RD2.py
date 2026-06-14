# coding=utf-8
"""
##########################################################################################
#    Program:       Falsifiable_ Predictions_DK-RD2
#    Main script:   Falsifiable_ Predictions_DK-RD2.py
#    Author:        Gabriel Martín del Campo Flores & Jazmin Olvera Zacarias
#    Contact:       gabemdelc@gmail.com                 Jaolveraz@ipn.mx. Mx
#    Created:       6-7/june/2026
#    License:       MIT License
##########################################################################################
##########################################################################################
#
# Falsifiable Predictions of the DK-RD2 Framework
# Observational Tests Beyond ΛCDM
#
# ----------------------------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------------------------
#
# This script implements the numerical pipeline used in the paper:
#
# "Falsifiable Predictions of the DK-RD2 Framework:
# Observational Tests Beyond ΛCDM"
#
# The objective is to generate observational predictions capable of confirming
# or falsifying the DK-RD2 framework through direct comparison with current and
# future cosmological observations.
#
# The effective gravitational coupling is:
#
# Gab(T,v) ≈ G0 · [1 + αDK · (v²/c²)(T0/T)]
#
# where:
#
# G0      = Newtonian gravitational constant
# αDK     = geometric projection factor
# v       = characteristic relativistic velocity scale
# T       = thermodynamic state of the system
# T0      = present-day CMB temperature
#
# For cosmological applications:
#
# T(z) ≈ T0 · (1 + z)
#
# ----------------------------------------------------------------------------------------
# Scientific Focus
# ----------------------------------------------------------------------------------------
#
# • Statistical calibration of αDK and r_d
# • DESI BAO + Cosmic Chronometer constraints
# • Expansion-history reconstruction
# • Emergent dark-energy-like behavior
# • Geometric gravitational-lensing signatures
# • Large-scale structure growth
# • Compact-halo amplification effects
# • Hubble-tension projection test
# • Model discrimination against alternative DESI explanations
# • Reproducible falsification roadmap
#
# ----------------------------------------------------------------------------------------
# Code Structure
# ----------------------------------------------------------------------------------------
#
# DK_RD2_Core.py
#   • Gab(T,v)
#   • T(z)
#   • H(z)
#   • Distance relations
#   • BAO observables
#   • Growth observables
#   • Lensing diagnostics
#   • Hubble-tension diagnostics
#   • Statistical likelihood calculations
#
# DK_RD2_CLASS.py
#   • Optional CMB consistency tools
#   • Auxiliary geometric calculations
#
# DK_RD2_Framework.py
#   • Dataset loading
#   • DK-RD2 predictions
#   • Figure generation
#   • Statistical diagnostics
#   • Reproducible outputs
#
# ----------------------------------------------------------------------------------------
# Core Outputs
# ----------------------------------------------------------------------------------------
#
# ✓ Calibrated αDK
# ✓ Calibrated sound horizon r_d
# ✓ BAO likelihood diagnostics
# ✓ χ² profiles
# ✓ H(z) reconstruction
# ✓ Effective equation of state w_eff(z)
# ✓ Lensing predictions
# ✓ fσ8(z) growth predictions
# ✓ Compact-halo amplification μeff(r)
# ✓ Effective H0 reconstructions
# ✓ Multi-model comparison diagnostics
#
# ----------------------------------------------------------------------------------------
# Figures Generated
# ----------------------------------------------------------------------------------------
#
# Figure 01 — Statistical Calibration of αDK
# Figure 02 — Likelihood Constraint on αDK
# Figure 03 — Effective Expansion History and Dynamical Dark-Energy Reconstruction
# Figure 04 — Geometric Gravitational-Lensing Signatures
# Figure 05 — Growth of Cosmic Structure
# Figure 06 — Compact-Halo Prediction and Effective Gravitational Amplification
# Figure 07 — Hubble-Tension Projection Test
# Figure 08 — Model Discrimination: DK-RD2 versus Alternative Explanations of Recent DESI Observations
#
# ----------------------------------------------------------------------------------------
# Outputs
# ----------------------------------------------------------------------------------------
#
# • Publication figures
# • Statistical tables
# • χ² diagnostics
# • Calibration summaries
# • Reproducibility metadata
#
# ----------------------------------------------------------------------------------------
# Falsifiability
# ----------------------------------------------------------------------------------------
#
# DK-RD2 predicts:
# • A statistically preferred non-zero αDK
# • Dynamical dark-energy-like evolution
# • Observable lensing deviations
# • Distinct structure-growth histories
# • Compact-halo amplification signatures
# • Thermodynamic projection effects in Hubble reconstructions
# • Measurable differences relative to competing DESI interpretations
#
# DK-RD2 is challenged if:
# • αDK → 0 within uncertainties
# • Expansion remains fully ΛCDM-compatible
# • No lensing deviations are observed
# • fσ8(z) follows ΛCDM without DK-RD2 signatures
# • Compact-halo amplification is absent
# • No Hubble-tension projection effects exist
#
# ----------------------------------------------------------------------------------------
# Scientific Interpretation
# ----------------------------------------------------------------------------------------
#
# Gravity is treated as an emergent thermodynamic-relativistic phenomenon whose
# observable behavior depends on the thermal, energetic, and dynamical state of
# the system.
#
# The framework is evaluated through a multi-scale observational falsification
# program spanning cosmological expansion, gravitational lensing, structure
# formation, compact systems, the Hubble-tension problem, and competing
# interpretations of recent DESI observations.
#
# ----------------------------------------------------------------------------------------
# Motto
# ----------------------------------------------------------------------------------------
#
# Gravity is not assumed — it emerges.
# GabE = mc²  —  Luludns = ∞Ψ
#
##########################################################################################
"""
from DK_RD2_Core import * # DK-RD2 Core Utilities – Constants, Functions, and Relativistic Dynamic Gravitational Engine
##########################################################################################


import matplotlib.pyplot as plt

Core_git_gabe = "https://github.com/gabemdelc/DK-RD2-Falsifiable-Predictions"
Core_zenodo = "Zenodo DOI: https://doi.org/10.5281/zenodo.20637136" # DK-RD2 Falsifiable Predictions

Core_autor_text = (
    f"{Core_author}. "
    "Reproducible from public DK-RD2 code.\n"
    f"GitHub: {Core_git_gabe} | Zenodo: {Core_zenodo}"
)

# Global linestyle & color for use in all Figures to models comparison
DK_RD2_color      = "blue"     # DK-RD2 main curves
LCDM_color        = "orange"

DK_LIGHT_color    = "#7fe7ff"
LCDM_LIGHT_color  = "#ffd280"

DATA_color        = "white"
ERROR_color       = "0.75"

SECONDARY_color   = "#b57cff"

GRID_color        = "0.22"
TEXT_color        = "0.92"

DK_RD2_linestyle="dashdot"
LCDM_linestyle="--"

# Global DK-RD2 calibrated projection factor cache
alpha_dk_best_global = None
rd_dk_best_global = None
alpha_dk_source_global = "not computed"


def calibrate_alpha_DK_global(
    desi_bao_mean_path: str,
    desi_bao_cov_path: str,
    *,
    force: bool = False,
    verbose: bool = True,
):
    """
    Compute alpha_DK and r_d once from DESI BAO + Cosmic Chronometer calibration.

    This value is the official paper calibration and must be reused by all figures.
    """

    global alpha_dk_best_global
    global rd_dk_best_global
    global alpha_dk_source_global

    if alpha_dk_best_global is not None and rd_dk_best_global is not None and not force:
        return float(alpha_dk_best_global), float(rd_dk_best_global)

    alpha_dk_best_global, rd_dk_best_global = compute_alpha_DK_from_DESI(
        desi_bao_mean_path=desi_bao_mean_path,
        desi_bao_cov_path=desi_bao_cov_path,
    )

    alpha_dk_source_global = "DESI BAO + Cosmic Chronometer global calibration"

    if verbose:
        print("============================================================")
        print("Official DK-RD2 global calibration")
        print(f"alpha_DK = {alpha_dk_best_global:.6f}")
        print(f"r_d      = {rd_dk_best_global:.6f} Mpc")
        print("Source   = DESI BAO + CC")
        print("============================================================")

    return float(alpha_dk_best_global), float(rd_dk_best_global)

def get_alpha_DK(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
):
    """
    Return the globally calibrated alpha_DK and r_d.

    If the values were not computed yet, compute them once.
    """

    if alpha_dk_best_global is None or rd_dk_best_global is None:
        if desi_bao_mean_path is None or desi_bao_cov_path is None:
            raise RuntimeError(
                "alpha_DK has not been calibrated yet. "
                "Call calibrate_alpha_DK_global(...) in main first."
            )

        return calibrate_alpha_DK_global(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
            force=False,
            verbose=True,
        )

    return float(alpha_dk_best_global), float(rd_dk_best_global)

def generate_figure01(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    Core_H0: float | None = None,
    Core_Omega_m_LCDM: float | None = None,
    Omega_L_LCDM: float | None = None,
    rd_planck: float | None = None,
    rd_bounds: tuple[float, float] | None = None,
    v_model=None,
    T_model=None,
):
    """
    Figure 01 — DESI Expansion and BAO Calibration of alpha_DK.
    Uses the official global DK-RD2 calibration computed once in main().
    """

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure01 requires desi_bao_mean_path and desi_bao_cov_path.")

    if alpha_dk_best_global is None or rd_dk_best_global is None:
        raise RuntimeError(
            "Global DK-RD2 calibration has not been computed. "
            "Run calibrate_alpha_DK_global() in main first."
        )

    if Core_H0 is None:
        Core_H0 = float(Core_Hubble_H0)

    if Core_Omega_m_LCDM is None:
        Core_Omega_m_LCDM = float(Core_OMEGA_M_LCDM)

    if Omega_L_LCDM is None:
        Omega_L_LCDM = float(Core_OMEGA_L_LCDM)

    if rd_planck is None:
        rd_planck = float(Core_rd_planck_mpc)

    if rd_bounds is None:
        rd_bounds = tuple(Core_rd_fit_bounds_mpc)

    alpha_dk_best = float(alpha_dk_best_global)
    rd_dk_best = float(rd_dk_best_global)

    file_fig = generate_evidence("image", 1)
    file_table = generate_evidence("table", 1)
    file_stats = file_table.replace(".csv", "_stats.csv")

    z_bao, bao_obs, bao_type, bao_cov, bao_cov_inv = load_desi_gaussian_bao(
        desi_bao_mean_path,
        desi_bao_cov_path,
    )

    sigma_bao = np.sqrt(np.diag(bao_cov))

    def E_DK_alpha(z_in, alpha_DK):
        z_arr = np.asarray(z_in, dtype=float)
        E_raw = E_Relativistic(
            z_arr,
            Core_Omega_m=None,
            Omega_L_value=None,
            v_model=v_model,
            T_model=T_model,
        )
        E2_alpha = 1.0 + float(alpha_DK) * (E_raw**2 - 1.0)
        return np.sqrt(np.clip(E2_alpha, 1e-300, None))

    def H_DK_alpha(z_in, alpha_DK):
        return float(Core_H0) * E_DK_alpha(z_in, alpha_DK)

    lcdm_fit = fit_rd_for_bao(
        z_bao,
        bao_obs,
        bao_type,
        bao_cov_inv,
        E_LCDM,
        Core_H0,
        rd_bounds=rd_bounds,
        Core_Omega_m=Core_Omega_m_LCDM,
        Core_Omega_L=Omega_L_LCDM,
    )

    bao_lcdm_planck = bao_distance_vector_over_rs(
        z_bao,
        bao_type,
        E_LCDM,
        Core_H0,
        rd_planck,
        Core_Omega_m=Core_Omega_m_LCDM,
        Core_Omega_L=Omega_L_LCDM,
    )

    chi2_lcdm_planck = chi2_gaussian_bao(
        bao_lcdm_planck,
        bao_obs,
        bao_cov_inv,
    )

    z_cc = float(Core_DESI_CC_z)
    H_cc = float(Core_DESI_CC_H)
    H_cc_err = float(Core_DESI_CC_H_err)

    H_lcdm_cc = float(H_LCDM(
        z_cc,
        Core_H0,
        Core_Omega_m_LCDM,
        Omega_L_LCDM,
    ))

    chi2_cc_lcdm = chi2_cosmic_chronometer_point(
        H_lcdm_cc,
        H_obs=H_cc,
        sigma_H=H_cc_err,
    )

    def dk_model_vector(rd_mpc, alpha_DK):
        return bao_distance_vector_over_rs(
            z_bao,
            bao_type,
            lambda zz: E_DK_alpha(zz, alpha_DK),
            Core_H0,
            rd_mpc,
        )

    bao_dk_best = dk_model_vector(rd_dk_best, alpha_dk_best)
    dk_residuals = bao_dk_best - bao_obs
    dk_pulls = dk_residuals / sigma_bao

    chi2_dk_bao = chi2_gaussian_bao(
        bao_dk_best,
        bao_obs,
        bao_cov_inv,
    )

    H_dk_cc = float(H_DK_alpha(z_cc, alpha_dk_best))

    chi2_cc_dk = chi2_cosmic_chronometer_point(
        H_dk_cc,
        H_obs=H_cc,
        sigma_H=H_cc_err,
    )

    z_plot = np.linspace(0.001, 2.5, 900)

    H_lcdm = H_LCDM(
        z_plot,
        Core_H0,
        Core_Omega_m_LCDM,
        Omega_L_LCDM,
    )

    H_dk = H_DK_alpha(z_plot, alpha_dk_best)

    z_bao_curve = np.linspace(0.05, 2.5, 600)

    DM_lcdm_curve = np.array([
        comoving_distance(
            zi,
            E_LCDM,
            Core_H0,
            Core_c_km_s=Core_c_km_s,
            Core_Omega_m=Core_Omega_m_LCDM,
            Core_Omega_L=Omega_L_LCDM,
        ) / lcdm_fit["best_rd"]
        for zi in z_bao_curve
    ])

    DH_lcdm_curve = np.array([
        (
            Core_c_km_s / H_LCDM(
                zi,
                Core_H0,
                Core_Omega_m_LCDM,
                Omega_L_LCDM,
            )
        ) / lcdm_fit["best_rd"]
        for zi in z_bao_curve
    ])

    DM_dk_curve = np.array([
        comoving_distance(
            zi,
            lambda zz: E_DK_alpha(zz, alpha_dk_best),
            Core_H0,
            Core_c_km_s=Core_c_km_s,
        ) / rd_dk_best
        for zi in z_bao_curve
    ])

    DH_dk_curve = np.array([
        (
            Core_c_km_s / H_DK_alpha(zi, alpha_dk_best)
        ) / rd_dk_best
        for zi in z_bao_curve
    ])

    table = pd.DataFrame({
        "z": z_bao,
        "observable": bao_type,
        "DESI_DR2_value": bao_obs,
        "sigma_diag": sigma_bao,
        "LCDM_reference_rd": bao_lcdm_planck,
        "LCDM_fit_rd": lcdm_fit["model"],
        "DKRD2_global_alpha_rd": bao_dk_best,
        "LCDM_fit_residual": lcdm_fit["residuals"],
        "DKRD2_global_alpha_residual": dk_residuals,
        "LCDM_fit_pull": lcdm_fit["pulls"],
        "DKRD2_global_alpha_pull": dk_pulls,
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="DESI_BAO_DR2_ALPHA_DK",
        figure_id=1,
        strict=False,
        index=False,
        meta={
            "rd_reference_mpc": rd_planck,
            "rd_lcdm_best_mpc": lcdm_fit["best_rd"],
            "rd_dkrd2_global_mpc": rd_dk_best,
            "alpha_dkrd2_global": alpha_dk_best,
            "chi2_lcdm_reference_bao": chi2_lcdm_planck,
            "chi2_lcdm_fit_bao": lcdm_fit["best_chi2"],
            "chi2_dkrd2_bao": chi2_dk_bao,
            "chi2_dkrd2_cc": chi2_cc_dk,
            "chi2_dkrd2_total": chi2_dk_bao + chi2_cc_dk,
            "H_cc_z": z_cc,
            "H_cc": H_cc,
            "H_cc_err": H_cc_err,
            "alpha_source": "alpha_DK and r_d are obtained from the global BAO + CC calibration.",
        },
    )

    fig, (ax_h, ax_bao) = plt.subplots(
        1,
        2,
        figsize=(15, 6.5),
        dpi=130,
    )

    ax_h.plot(
        z_plot,
        H_lcdm,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=2.1,
        label=r"$\Lambda$CDM reference",
    )

    ax_h.plot(
        z_plot,
        H_dk,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.3,
        label=rf"DK-RD2 calibrated coupling: $\alpha_{{\rm DK}}={alpha_dk_best:.4f}$",
    )

    ax_h.errorbar(
        [z_cc],
        [H_cc],
        yerr=[H_cc_err],
        fmt="o",
        color="black",
        capsize=4,
        label=rf"DESI DR1 CC: $H({z_cc:.2f})={H_cc:.2f}\pm{H_cc_err:.2f}$",
    )

    ax_h.text(
        0.03,
        0.96,
        r"DK-RD2 global calibration:" "\n"
        r"$E_{\rm DK,\alpha}^2(z)=1+\alpha_{\rm DK}[E_{\rm DK,raw}^2(z)-1]$" "\n"
        r"$\alpha_{\rm DK}$ and $r_d$ are obtained from the global BAO + CC calibration.",
        transform=ax_h.transAxes,
        fontsize=8,
        va="top",
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="white",
            edgecolor="0.5",
            alpha=0.85,
        ),
    )

    ax_h.set_xlabel("Redshift z")
    ax_h.set_ylabel(r"$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]")
    ax_h.set_title(r"Panel A — Expansion Rate")
    ax_h.grid(alpha=0.3)
    ax_h.legend(fontsize=8, loc="best")

    ax_bao.plot(
        z_bao_curve,
        DM_lcdm_curve,
        color=LCDM_color,
        linestyle="-",
        linewidth=1.8,
        label=rf"$\Lambda$CDM $D_M/r_d$, best $r_d={lcdm_fit['best_rd']:.2f}$ Mpc",
    )

    ax_bao.plot(
        z_bao_curve,
        DH_lcdm_curve,
        color=LCDM_color,
        linestyle="--",
        linewidth=1.8,
        label=r"$\Lambda$CDM $D_H/r_d$",
    )

    ax_bao.plot(
        z_bao_curve,
        DM_dk_curve,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=2.0,
        label=rf"DK-RD2 $D_M/r_d$, $r_d={rd_dk_best:.2f}$ Mpc",
    )

    ax_bao.plot(
        z_bao_curve,
        DH_dk_curve,
        color=DK_RD2_color,
        linestyle="--",
        linewidth=2.0,
        label=r"DK-RD2 $D_H/r_d$",
    )

    mask_dm = bao_type == "DM_over_rs"
    mask_dh = bao_type == "DH_over_rs"
    mask_dv = bao_type == "DV_over_rs"

    ax_bao.errorbar(
        z_bao[mask_dm],
        bao_obs[mask_dm],
        yerr=sigma_bao[mask_dm],
        fmt="s",
        color="black",
        capsize=3,
        label=r"DESI DR2 $D_M/r_d$",
    )

    ax_bao.errorbar(
        z_bao[mask_dh],
        bao_obs[mask_dh],
        yerr=sigma_bao[mask_dh],
        fmt="^",
        color="0.35",
        capsize=3,
        label=r"DESI DR2 $D_H/r_d$",
    )

    ax_bao.errorbar(
        z_bao[mask_dv],
        bao_obs[mask_dv],
        yerr=sigma_bao[mask_dv],
        fmt="o",
        color="0.15",
        capsize=3,
        label=r"DESI DR2 $D_V/r_d$",
    )

    delta_chi2 = (chi2_dk_bao + chi2_cc_dk) - (lcdm_fit["best_chi2"] + chi2_cc_lcdm)

    fit_text = (
        r"DK-RD2 global fit:" "\n"
        rf"$\chi^2_{{\rm BAO,DK}} = {chi2_dk_bao:.2f}$" "\n"
        rf"$\chi^2_{{\rm CC,DK}} = {chi2_cc_dk:.3f}$" "\n"
        rf"$\chi^2_{{\rm total,DK}} = {(chi2_dk_bao + chi2_cc_dk):.2f}$" "\n"
        "\n"
        rf"$\chi^2_{{\Lambda CDM}} = {(lcdm_fit['best_chi2'] + chi2_cc_lcdm):.2f}$" "\n"
        rf"$\Delta \chi^2_{{\rm DK-\Lambda CDM}} = {delta_chi2:.2f}$"
    )

    ax_bao.text(
        0.70,
        0.77,
        fit_text,
        transform=ax_bao.transAxes,
        fontsize=9.0,
        va="top",
        ha="left",
        linespacing=1.25,
        bbox=dict(
            boxstyle="round,pad=0.55",
            facecolor="white",
            edgecolor="black",
            linewidth=1.1,
            alpha=0.92,
        ),
    )

    ax_bao.set_xlabel("Redshift z")
    ax_bao.set_ylabel(r"BAO compressed observable")
    ax_bao.set_title(r"Panel B — DESI DR2 Gaussian BAO")
    ax_bao.grid(alpha=0.3)
    ax_bao.legend(fontsize=7, loc="best")

    fig.suptitle(
        rf"Figure 01 — DESI Expansion and BAO Calibration of $\alpha_{{\rm DK}}$"
        "\n"
        rf"$\alpha_{{\rm DK}} = {alpha_dk_best:.4f},\quad r_d = {rd_dk_best:.2f}\,\mathrm{{Mpc}}$",
        fontsize=14,
    )

    fig.text(
        0.5,
        0.043,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=6,
        color=DK_RD2_color,
    )

    plt.tight_layout(rect=(0, 0.06, 1, 1.0))
    plt.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    stats = pd.DataFrame([
        {
            "model": "ΛCDM reference-calibrated",
            "rd_mpc": float(rd_planck),
            "alpha_DK": np.nan,
            "chi2_BAO": float(chi2_lcdm_planck),
            "chi2_CC": float(chi2_cc_lcdm),
            "chi2_total": float(chi2_lcdm_planck + chi2_cc_lcdm),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "reference_rd_fixed",
        },
        {
            "model": "ΛCDM BAO-calibrated",
            "rd_mpc": float(lcdm_fit["best_rd"]),
            "alpha_DK": np.nan,
            "chi2_BAO": float(lcdm_fit["best_chi2"]),
            "chi2_CC": float(chi2_cc_lcdm),
            "chi2_total": float(lcdm_fit["best_chi2"] + chi2_cc_lcdm),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "rd_free",
        },
        {
            "model": "DK-RD2 global calibrated coupling",
            "rd_mpc": float(rd_dk_best),
            "alpha_DK": float(alpha_dk_best),
            "chi2_BAO": float(chi2_dk_bao),
            "chi2_CC": float(chi2_cc_dk),
            "chi2_total": float(chi2_dk_bao + chi2_cc_dk),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "global_alpha_rd_from_main",
        },
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=1,
        fit_mode="DESI_BAO_PLUS_CC_GLOBAL_ALPHA_DK",
        index=False,
    )

    return file_fig, file_table, file_stats

def generate_figure02(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    n_alpha: int = 220,
    alpha_min: float | None = None,
    alpha_max: float | None = None,
    rd_bounds: tuple[float, float] | None = None,
    Core_H0: float | None = None,
    Omega_L_DK: float | None = None,
    v_model=None,
    T_model=None,
):
    """
    Figure 02 — Statistical Constraint on the DK Projection Factor.

    Purpose
    -------
    Demonstrate that alpha_DK is not arbitrary by profiling
    chi2_BAO+CC(alpha_DK, r_d) over r_d for each alpha_DK.

    Notes
    -----
    The official alpha_DK and r_d values are computed once in main()
    through calibrate_alpha_DK_global(...) and reused here.
    This figure scans alpha_DK only to show the likelihood structure.
    """

    from scipy.optimize import minimize_scalar

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure02 requires desi_bao_mean_path and desi_bao_cov_path.")

    if alpha_dk_best_global is None or rd_dk_best_global is None:
        raise RuntimeError(
            "Global DK-RD2 calibration has not been computed. "
            "Run calibrate_alpha_DK_global() in main first."
        )

    if Core_H0 is None:
        Core_H0 = float(Core_Hubble_H0)

    if rd_bounds is None:
        rd_bounds = tuple(Core_rd_fit_bounds_mpc)

    if alpha_min is None:
        alpha_min = float(Core_alpha_DK_bounds[0])

    if alpha_max is None:
        alpha_max = float(Core_alpha_DK_bounds[1])

    Omega_L_DK = None

    alpha_best_profile = float(alpha_dk_best_global)
    rd_best_profile = float(rd_dk_best_global)

    file_fig = generate_evidence("image", 2)
    file_table = generate_evidence("table", 2)
    file_stats = file_table.replace(".csv", "_stats.csv")

    z_bao, bao_obs, bao_type, bao_cov, bao_cov_inv = load_desi_gaussian_bao(
        desi_bao_mean_path,
        desi_bao_cov_path,
    )

    z_cc = float(Core_DESI_CC_z)
    H_cc = float(Core_DESI_CC_H)
    H_cc_err = float(Core_DESI_CC_H_err)

    def E_DK_alpha(z_in, alpha_DK):
        z_arr = np.asarray(z_in, dtype=float)

        E_raw = E_Relativistic(
            z_arr,
            Core_Omega_m=None,
            Omega_L_value=Omega_L_DK,
            v_model=v_model,
            T_model=T_model,
        )

        E2_alpha = 1.0 + float(alpha_DK) * (E_raw**2 - 1.0)
        return np.sqrt(np.clip(E2_alpha, 1e-300, None))

    def H_DK_alpha(z_in, alpha_DK):
        return float(Core_H0) * E_DK_alpha(z_in, alpha_DK)

    def model_vector_for_alpha_rd(alpha_DK, rd_mpc):
        return bao_distance_vector_over_rs(
            z_bao,
            bao_type,
            lambda zz: E_DK_alpha(zz, alpha_DK),
            Core_H0,
            rd_mpc,
        )

    def chi2_total_for_alpha_rd(alpha_DK, rd_mpc):
        bao_model = model_vector_for_alpha_rd(alpha_DK, rd_mpc)

        chi2_bao = chi2_gaussian_bao(
            bao_model,
            bao_obs,
            bao_cov_inv,
        )

        H_model_cc = float(H_DK_alpha(z_cc, alpha_DK))

        chi2_cc = chi2_cosmic_chronometer_point(
            H_model_cc,
            H_obs=H_cc,
            sigma_H=H_cc_err,
        )

        return float(chi2_bao + chi2_cc), float(chi2_bao), float(chi2_cc)

    alpha_grid = np.linspace(float(alpha_min), float(alpha_max), int(n_alpha))
    rows = []

    for alpha_val in alpha_grid:

        def objective_rd(rd_mpc):
            chi2_total, _, _ = chi2_total_for_alpha_rd(alpha_val, rd_mpc)
            return chi2_total

        res_rd = minimize_scalar(
            objective_rd,
            bounds=rd_bounds,
            method="bounded",
        )

        rd_best = float(res_rd.x)
        chi2_total, chi2_bao, chi2_cc = chi2_total_for_alpha_rd(
            alpha_val,
            rd_best,
        )

        rows.append({
            "alpha_DK": float(alpha_val),
            "rd_profiled_mpc": rd_best,
            "chi2_BAO": chi2_bao,
            "chi2_CC": chi2_cc,
            "chi2_total": chi2_total,
            "profile_success": bool(res_rd.success),
        })

    profile_df = pd.DataFrame(rows)

    idx_min = int(profile_df["chi2_total"].idxmin())
    alpha_best_grid = float(profile_df.loc[idx_min, "alpha_DK"])
    rd_best_grid = float(profile_df.loc[idx_min, "rd_profiled_mpc"])
    chi2_min_grid = float(profile_df.loc[idx_min, "chi2_total"])

    chi2_min, chi2_bao_min, chi2_cc_min = chi2_total_for_alpha_rd(
        alpha_best_profile,
        rd_best_profile,
    )

    profile_df["delta_chi2"] = profile_df["chi2_total"] - chi2_min

    delta_floor = 1.0e-2
    profile_df["delta_chi2_plot"] = np.clip(
        profile_df["delta_chi2"].to_numpy(dtype=float),
        delta_floor,
        None,
    )

    delta_1sigma = 1.0
    delta_2sigma = 4.0
    delta_3sigma = 9.0

    def interval_from_delta(delta_value):
        mask = profile_df["delta_chi2"].values <= float(delta_value)
        if not np.any(mask):
            return np.nan, np.nan

        alpha_vals = profile_df["alpha_DK"].values[mask]
        return float(np.min(alpha_vals)), float(np.max(alpha_vals))

    alpha_1sig_low, alpha_1sig_high = interval_from_delta(delta_1sigma)
    alpha_2sig_low, alpha_2sig_high = interval_from_delta(delta_2sigma)
    alpha_3sig_low, alpha_3sig_high = interval_from_delta(delta_3sigma)

    def diagnostic_fixed_alpha(alpha_fixed):
        def objective_rd(rd_mpc):
            chi2_total, _, _ = chi2_total_for_alpha_rd(alpha_fixed, rd_mpc)
            return chi2_total

        res = minimize_scalar(
            objective_rd,
            bounds=rd_bounds,
            method="bounded",
        )

        chi2_total, chi2_bao, chi2_cc = chi2_total_for_alpha_rd(
            alpha_fixed,
            float(res.x),
        )

        return {
            "alpha_DK": float(alpha_fixed),
            "rd_mpc": float(res.x),
            "chi2_BAO": chi2_bao,
            "chi2_CC": chi2_cc,
            "chi2_total": chi2_total,
            "success": bool(res.success),
        }

    diag_alpha0 = diagnostic_fixed_alpha(0.0)
    diag_alpha1 = diagnostic_fixed_alpha(1.0)

    dkrd2_to_csv(
        profile_df,
        file_table,
        table_kind="FIG02_ALPHA_DK_CHI2_PROFILE",
        figure_id=2,
        strict=False,
        index=False,
        meta={
            "alpha_DK_official": alpha_best_profile,
            "rd_DK_official_mpc": rd_best_profile,
            "chi2_official": chi2_min,
            "chi2_official_BAO": chi2_bao_min,
            "chi2_official_CC": chi2_cc_min,
            "alpha_grid_min": alpha_min,
            "alpha_grid_max": alpha_max,
            "rd_bounds_mpc": rd_bounds,
            "profiled_parameter": "rd_mpc",
            "profile_curve_purpose": (
                "Visual chi2(alpha_DK) diagnostic; official alpha_DK and r_d "
                "are obtained from the global BAO + CC calibration."
            ),
            "alpha_source": "Official global DESI BAO + CC calibration",
            "cc_anchor_note": "Single low-redshift cosmic chronometer anchor.",
            "dk_background": "Omega_L_value=None; no explicit Lambda term.",
            "diagnostic_alpha0_chi2_total": diag_alpha0["chi2_total"],
            "diagnostic_alpha0_rd_mpc": diag_alpha0["rd_mpc"],
            "diagnostic_alpha1_chi2_total": diag_alpha1["chi2_total"],
            "diagnostic_alpha1_rd_mpc": diag_alpha1["rd_mpc"],
            "grid_min_alpha_DK": alpha_best_grid,
            "grid_min_rd_mpc": rd_best_grid,
            "grid_min_chi2_total": chi2_min_grid,
        },
    )

    stats = pd.DataFrame([
        {
            "model": "DK-RD2 official global calibration",
            "alpha_DK": alpha_best_profile,
            "rd_mpc": rd_best_profile,
            "chi2_total": chi2_min,
            "chi2_BAO": chi2_bao_min,
            "chi2_CC": chi2_cc_min,
            "delta_chi2": 0.0,
            "fit_mode": "official_global_alpha_rd_from_main",
        },
        {
            "model": "DK-RD2 profiled grid minimum",
            "alpha_DK": alpha_best_grid,
            "rd_mpc": rd_best_grid,
            "chi2_total": chi2_min_grid,
            "chi2_BAO": np.nan,
            "chi2_CC": np.nan,
            "delta_chi2": chi2_min_grid - chi2_min,
            "fit_mode": "grid_profile_minimum",
        },
        {
            "model": "DK-RD2 diagnostic alpha=0",
            "alpha_DK": diag_alpha0["alpha_DK"],
            "rd_mpc": diag_alpha0["rd_mpc"],
            "chi2_total": diag_alpha0["chi2_total"],
            "chi2_BAO": diag_alpha0["chi2_BAO"],
            "chi2_CC": diag_alpha0["chi2_CC"],
            "delta_chi2": diag_alpha0["chi2_total"] - chi2_min,
            "fit_mode": "alpha_fixed_0_rd_profiled",
        },
        {
            "model": "DK-RD2 diagnostic alpha=1",
            "alpha_DK": diag_alpha1["alpha_DK"],
            "rd_mpc": diag_alpha1["rd_mpc"],
            "chi2_total": diag_alpha1["chi2_total"],
            "chi2_BAO": diag_alpha1["chi2_BAO"],
            "chi2_CC": diag_alpha1["chi2_CC"],
            "delta_chi2": diag_alpha1["chi2_total"] - chi2_min,
            "fit_mode": "alpha_fixed_1_rd_profiled",
        },
        {
            "model": "DK-RD2 alpha_DK confidence profile",
            "alpha_DK": alpha_best_profile,
            "rd_mpc": rd_best_profile,
            "chi2_total": chi2_min,
            "chi2_BAO": chi2_bao_min,
            "chi2_CC": chi2_cc_min,
            "alpha_1sigma_low_delta_chi2_1": alpha_1sig_low,
            "alpha_1sigma_high_delta_chi2_1": alpha_1sig_high,
            "alpha_2sigma_low_delta_chi2_4": alpha_2sig_low,
            "alpha_2sigma_high_delta_chi2_4": alpha_2sig_high,
            "alpha_3sigma_low_delta_chi2_9": alpha_3sig_low,
            "alpha_3sigma_high_delta_chi2_9": alpha_3sig_high,
            "n_alpha_grid": int(n_alpha),
            "rd_bounds_low_mpc": float(rd_bounds[0]),
            "rd_bounds_high_mpc": float(rd_bounds[1]),
            "cc_anchor_z": z_cc,
            "cc_anchor_H": H_cc,
            "cc_anchor_sigma_H": H_cc_err,
            "fit_mode": "profiled_rd_for_each_alpha",
        },
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=2,
        fit_mode="ALPHA_DK_CHI2_PROFILE_RD_PROFILED",
        index=False,
    )

    fig, (ax_chi, ax_rd) = plt.subplots(
        2,
        1,
        num="Figure 02 — DK Projection Factor Constraint",
        figsize=(10.8, 8.4),
        dpi=130,
        sharex=True,
        gridspec_kw={"height_ratios": [2.55, 1.25], "hspace": 0.10},
    )

    if np.isfinite(alpha_3sig_low) and np.isfinite(alpha_3sig_high):
        ax_chi.axvspan(
            alpha_3sig_low,
            alpha_3sig_high,
            color=DK_RD2_color,
            alpha=0.045,
        )

    if np.isfinite(alpha_2sig_low) and np.isfinite(alpha_2sig_high):
        ax_chi.axvspan(
            alpha_2sig_low,
            alpha_2sig_high,
            color=DK_RD2_color,
            alpha=0.075,
        )

    if np.isfinite(alpha_1sig_low) and np.isfinite(alpha_1sig_high):
        ax_chi.axvspan(
            alpha_1sig_low,
            alpha_1sig_high,
            color=DK_RD2_color,
            alpha=0.12,
        )

    ax_chi.plot(
        profile_df["alpha_DK"],
        profile_df["delta_chi2_plot"],
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.4,
        label=r"Profiled $\Delta\chi^2(\alpha_{\rm DK})$",
    )

    ax_chi.axhline(
        delta_1sigma,
        color="black",
        linestyle=":",
        linewidth=1.1,
        alpha=0.85,
        label=r"$1\sigma$: $\Delta\chi^2=1$",
    )

    ax_chi.axhline(
        delta_2sigma,
        color="0.35",
        linestyle="--",
        linewidth=1.0,
        alpha=0.75,
        label=r"$2\sigma$: $\Delta\chi^2=4$",
    )

    ax_chi.axhline(
        delta_3sigma,
        color="0.55",
        linestyle="-.",
        linewidth=1.0,
        alpha=0.70,
        label=r"$3\sigma$: $\Delta\chi^2=9$",
    )

    ax_chi.axvline(
        alpha_best_profile,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=1.6,
        alpha=0.95,
    )

    ax_chi.scatter(
        [alpha_best_profile],
        [delta_floor],
        color=DK_RD2_color,
        s=45,
        zorder=5,
    )

    ax_chi.set_yscale("log")
    ax_chi.set_ylim(delta_floor * 0.8, profile_df["delta_chi2_plot"].max() * 1.35)
    ax_chi.set_ylabel(r"$\Delta\chi^2(\alpha_{\rm DK})$")
    ax_chi.set_title(
        rf"Figure 02 — Statistical Constraint on the DK Projection Factor"
        "\n"
        rf"Official DESI BAO + CC Calibration: "
        rf"$\alpha_{{\rm DK}}={alpha_best_profile:.4f}$, "
        rf"$r_d={rd_best_profile:.2f}\,\mathrm{{Mpc}}$, "
        rf"$\chi^2_{{min}}={chi2_min:.2f}$",
        fontsize=14,
    )

    ax_chi.grid(alpha=0.30, which="both")
    ax_chi.legend(fontsize=8, loc="upper right")

    result_text = (
        r"Official global calibration:" "\n"
        rf"$\alpha_{{\rm DK}}={alpha_best_profile:.4f}$" "\n"
        rf"$r_d={rd_best_profile:.2f}\,\mathrm{{Mpc}}$" "\n"
        "\n"
        rf"$\chi^2_{{\rm BAO}}={chi2_bao_min:.2f}$" "\n"
        rf"$\chi^2_{{\rm CC}}={chi2_cc_min:.3f}$" "\n"
        rf"$\chi^2_{{\rm tot}}={chi2_min:.2f}$"
    )

    ax_chi.text(
        0.03,
        0.30,
        result_text,
        transform=ax_chi.transAxes,
        ha="left",
        va="top",
        fontsize=5,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.20",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.90,
        ),
    )

    method_text = (
        r"For each $\alpha_{\rm DK}$, $r_d$ is profiled by minimizing"
        "\n"
        r"$\chi^2_{\rm BAO+CC}(\alpha_{\rm DK},r_d)$."
        "\n"
        r"The minimum is not imposed; it emerges from the likelihood."
        "\n"
        r"Diagnostics include $\alpha_{\rm DK}=0$ and $\alpha_{\rm DK}=1$."
    )

    ax_chi.text(
        0.50,
        0.055,
        method_text,
        transform=ax_chi.transAxes,
        ha="center",
        va="bottom",
        fontsize=5,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.18",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.78,
        ),
    )

    ax_rd.plot(
        profile_df["alpha_DK"],
        profile_df["rd_profiled_mpc"],
        color="black",
        linestyle="-",
        linewidth=1.8,
        label=r"Profiled $r_d(\alpha_{\rm DK})$",
    )

    ax_rd.axvline(
        alpha_best_profile,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=1.4,
        alpha=0.90,
    )

    ax_rd.axhline(
        rd_best_profile,
        color="0.35",
        linestyle="--",
        linewidth=1.2,
        alpha=0.80,
    )

    ax_rd.set_ylim(
        float(rd_bounds[0]) - 2.0,
        float(rd_bounds[1]) + 2.0,
    )

    ax_rd.set_xlabel(r"Projection factor $\alpha_{\rm DK}$", labelpad=8)
    ax_rd.set_ylabel(r"Profiled $r_d$ [Mpc]")
    ax_rd.grid(alpha=0.30)
    ax_rd.legend(fontsize=8, loc="best")

    fig.text(
        0.5,
        0.030,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=5.4,
        color=DK_RD2_color,
    )

    fig.subplots_adjust(
        top=0.88,
        bottom=0.16,
        hspace=0.14,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats

def generate_figure03(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
):
    """
    Figure 03 — DK-RD2 Coupling Structure and Cosmological Projection.

    This figure uses the official global alpha_DK and r_d calibration
    computed once in main(). No calibration is performed here.
    """

    if alpha_dk_best_global is None or rd_dk_best_global is None:
        raise RuntimeError(
            "Global DK-RD2 calibration has not been computed. "
            "Run calibrate_alpha_DK_global() in main first."
        )

    alpha_DK = float(alpha_dk_best_global)
    rd_DK = float(rd_dk_best_global)

    file_fig = generate_evidence("image", 3)
    file_table = generate_evidence("table", 3)
    file_stats = file_table.replace(".csv", "_stats.csv")

    # ------------------------------------------------------------
    # Panel A data — temperature coupling structure.
    # ------------------------------------------------------------
    T_vals = np.logspace(
        np.log10(Core_TCMB_K),
        np.log10(Core_T_plot_max),
        700,
    )

    beta = 1.0
    v_rel = beta * Core_c_light

    Gab_max = Gab(T_vals, v_rel)
    Gab_max_over_G0 = Gab_max / Core_G0

    thermo_rel_excess = Gab_max_over_G0 - 1.0

    Gab_eff_over_G0 = 1.0 + alpha_DK * thermo_rel_excess
    Gab_eff = Core_G0 * Gab_eff_over_G0

    projected_excess_over_G0 = Gab_eff_over_G0 - 1.0
    unprojected_excess_over_G0 = Gab_max_over_G0 - Gab_eff_over_G0

    projection_ratio_excess = np.where(
        thermo_rel_excess > 0.0,
        projected_excess_over_G0 / thermo_rel_excess,
        np.nan,
    )

    # ------------------------------------------------------------
    # Panel B data — cosmological redshift projection.
    # ------------------------------------------------------------
    z_cosmo = np.linspace(0.0, 3.0, 700)

    T_cosmo = Core_TCMB_K * (1.0 + z_cosmo)
    thermal_projection = Core_TCMB_K / T_cosmo

    delta_G_max_cosmo_over_G0 = beta**2 * thermal_projection
    delta_G_eff_cosmo_over_G0 = alpha_DK * beta**2 * thermal_projection

    Gab_max_cosmo_over_G0 = 1.0 + delta_G_max_cosmo_over_G0
    Gab_eff_cosmo_over_G0 = 1.0 + delta_G_eff_cosmo_over_G0

    # ------------------------------------------------------------
    # Save evidence table.
    # ------------------------------------------------------------
    n_rows = max(len(T_vals), len(z_cosmo))

    def pad_array(values, n=n_rows):
        arr = np.asarray(values, dtype=float)
        out = np.full(n, np.nan, dtype=float)
        out[:len(arr)] = arr
        return out

    df = pd.DataFrame({
        "Temperature_K": pad_array(T_vals),
        "beta_temperature_panel": pad_array(np.full_like(T_vals, beta)),
        "alpha_DK_temperature_panel": pad_array(np.full_like(T_vals, alpha_DK)),
        "rd_DK_Mpc_temperature_panel": pad_array(np.full_like(T_vals, rd_DK)),

        "Gab_max_over_G0_temperature": pad_array(Gab_max_over_G0),
        "Gab_eff_over_G0_temperature": pad_array(Gab_eff_over_G0),
        "thermo_rel_excess_over_G0_temperature": pad_array(thermo_rel_excess),
        "projected_excess_over_G0_temperature": pad_array(projected_excess_over_G0),
        "unprojected_excess_over_G0_temperature": pad_array(unprojected_excess_over_G0),
        "projection_ratio_excess_temperature": pad_array(projection_ratio_excess),
        "Gab_max_SI_temperature": pad_array(Gab_max),
        "Gab_eff_SI_temperature": pad_array(Gab_eff),

        "z_cosmological_projection": pad_array(z_cosmo),
        "T_cosmological_K": pad_array(T_cosmo),
        "T0_over_Tz": pad_array(thermal_projection),
        "delta_G_max_cosmo_over_G0": pad_array(delta_G_max_cosmo_over_G0),
        "delta_G_eff_cosmo_over_G0": pad_array(delta_G_eff_cosmo_over_G0),
        "Gab_max_cosmo_over_G0": pad_array(Gab_max_cosmo_over_G0),
        "Gab_eff_cosmo_over_G0": pad_array(Gab_eff_cosmo_over_G0),
    })

    dkrd2_to_csv(
        df,
        file_table,
        table_kind="FIG03_COUPLING_STRUCTURE_AND_COSMOLOGICAL_PROJECTION",
        figure_id=3,
        strict=False,
        index=False,
        meta={
            "alpha_DK": alpha_DK,
            "rd_DK_Mpc": rd_DK,
            "velocity_reference": "beta = v/c = 1",
            "temperature_range_K": f"{Core_TCMB_K} to {Core_T_plot_max}",
            "redshift_projection_range": f"{float(z_cosmo.min())} to {float(z_cosmo.max())}",
            "T0_source": "Core_TCMB_K from DK_RD2_Core.py",
            "G0_source": "Core_G0 from DK_RD2_Core.py",
            "c_source": "Core_c_light from DK_RD2_Core.py",
            "alpha_source": "Official global DESI BAO + Cosmic Chronometer calibration computed once in main().",
            "cosmological_projection": "T(z)=T0(1+z), T0/T(z)=1/(1+z)",
            "interpretation": (
                "Panel A shows the projected excess coupling structure. "
                "Panel B shows the cosmological redshift suppression of the DK correction. "
                "No parameter is fitted in this figure."
            ),
        },
    )

    # ------------------------------------------------------------
    # Plot configuration.
    # ------------------------------------------------------------
    fig, (ax_T, ax_z) = plt.subplots(
        2,
        1,
        figsize=(12.2, 10.0),
        dpi=130,
        gridspec_kw={
            "height_ratios": [2.35, 1.15],
            "hspace": 0.36,
        },
    )

    try:
        fig.canvas.manager.set_window_title("DK-RD2 Figures")
    except Exception:
        pass

    title_fs = 15
    panel_title_fs = 11.5
    label_fs = 10
    tick_fs = 8.5
    legend_fs = 7.0
    box_fs = 7.2
    footer_fs = 5.4

    # ============================================================
    # Panel A — Temperature-dependent coupling structure.
    # ============================================================
    ax_T.plot(
        T_vals,
        Gab_max_over_G0,
        color="red",
        linestyle=LCDM_linestyle,
        linewidth=2.1,
        label=r"$G_{ab}^{max}/G_0$  $(\beta=1,\ \alpha_{\rm DK}=1)$",
    )

    ax_T.plot(
        T_vals,
        Gab_eff_over_G0,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.4,
        label=rf"$G_{{ab}}^{{eff}}/G_0$  $(\alpha_{{\rm DK}}={alpha_DK:.4f})$",
    )

    ax_T.fill_between(
        T_vals,
        Gab_eff_over_G0,
        Gab_max_over_G0,
        alpha=0.22,
        label=r"Unprojected excess coupling",
    )

    ax_T.axhline(
        1.0,
        color="black",
        linewidth=1.0,
        alpha=0.80,
        label=r"Newtonian limit $G_0$",
    )

    ax_T.axhline(
        2.0,
        color="red",
        linestyle=":",
        linewidth=1.0,
        alpha=0.60,
        label=r"$G_{ab}^{max}\to 2G_0$  $(\beta=1,\ T\to T_0)$",
    )

    ax_T.axvline(
        Core_TCMB_K,
        color="black",
        linestyle=":",
        linewidth=1.2,
        alpha=0.85,
        label=rf"$T_0={Core_TCMB_K:.4f}\,\mathrm{{K}}$",
    )

    ax_T.set_ylim(0.995, 2.05)
    ax_T.set_xscale("log")
    ax_T.set_xlabel("Temperature T [K]", fontsize=label_fs)
    ax_T.set_ylabel(r"Normalized coupling $G_{ab}/G_0$", fontsize=label_fs)
    ax_T.set_title(
        r"Panel A — Thermodynamic–relativistic coupling structure",
        fontsize=panel_title_fs,
        pad=6,
    )

    ax_T.tick_params(axis="both", labelsize=tick_fs)
    ax_T.grid(True, which="both", linestyle="--", alpha=0.26)

    # ------------------------------------------------------------
    # Secondary axis — projected excess fraction.
    # ------------------------------------------------------------
    ax_frac = ax_T.twinx()

    ax_frac.plot(
        T_vals,
        projection_ratio_excess,
        color="black",
        linestyle=":",
        alpha=0.68,
        linewidth=1.3,
        label=r"Projected excess fraction",
    )

    ax_frac.set_ylim(0.0, 1.05)
    ax_frac.set_ylabel(r"Projected excess fraction", fontsize=label_fs)
    ax_frac.tick_params(axis="y", labelsize=tick_fs)

    lines_1, labels_1 = ax_T.get_legend_handles_labels()
    lines_2, labels_2 = ax_frac.get_legend_handles_labels()

    ax_T.legend(
        lines_1 + lines_2,
        labels_1 + labels_2,
        fontsize=legend_fs,
        loc="upper right",
        framealpha=0.90,
    )

    # ------------------------------------------------------------
    # Equation box.
    # ------------------------------------------------------------
    eq_text = (
        r"$G_{ab}^{max}=G_0\left[1+\frac{v^2}{c^2}\frac{T_0}{T}\right]$" "\n"
        r"$G_{ab}^{eff}=G_0\left[1+\alpha_{\rm DK}\frac{v^2}{c^2}\frac{T_0}{T}\right]$" "\n"
        rf"$\beta=v/c=1,\quad \alpha_{{\rm DK}}={alpha_DK:.4f}$"
    )

    ax_T.text(
        0.065,
        0.945,
        eq_text,
        transform=ax_T.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs,
        linespacing=1.08,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="black",
            linewidth=0.8,
            alpha=0.88,
        ),
    )

    insight_text = (
        "Projection interpretation:\n"
        f"{alpha_DK * 100:.2f}% of the coupling excess is projected into observable gravity.\n"
        "In the beta = 1, T -> T0 limit: "
        r"$G_{ab}^{max}\to 2G_0$."
    )

    ax_T.text(
        0.31,
        0.43,
        insight_text,
        transform=ax_T.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs,
        linespacing=1.08,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.8,
            alpha=0.88,
        ),
    )

    # ============================================================
    # Panel B — Cosmological thermal projection.
    # ============================================================
    ax_z.plot(
        z_cosmo,
        thermal_projection,
        color="black",
        linestyle="--",
        linewidth=1.7,
        label=r"$T_0/T(z)=1/(1+z)$",
    )

    ax_z.plot(
        z_cosmo,
        delta_G_eff_cosmo_over_G0,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=rf"$\Delta G_{{eff}}/G_0=\alpha_{{\rm DK}}/(1+z)$  $(\beta=1)$",
    )

    ax_z.plot(
        z_cosmo,
        Gab_eff_cosmo_over_G0,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=1.8,
        alpha=0.72,
        label=r"$G_{ab}^{eff}(z)/G_0=1+\Delta G_{eff}/G_0$",
    )

    ax_z.axhline(
        1.0,
        color="black",
        linewidth=0.85,
        alpha=0.45,
    )

    ax_z.set_xlabel("Redshift z", fontsize=label_fs)
    ax_z.set_ylabel("Projected contribution", fontsize=label_fs)
    ax_z.set_title(
        r"Panel B — Cosmological redshift projection of the DK correction",
        fontsize=panel_title_fs,
        pad=6,
    )

    ax_z.tick_params(axis="both", labelsize=tick_fs)
    ax_z.grid(alpha=0.28)

    ax_z.legend(
        fontsize=legend_fs,
        loc="center right",
        framealpha=0.88,
    )

    projection_note = (
        r"$T(z)=T_0(1+z)$  $\Rightarrow$  $T_0/T(z)=1/(1+z)$" "\n"
        "The DK correction is suppressed at high redshift\n"
        "and becomes maximal toward the late Universe."
    )

    ax_z.text(
        0.04,
        0.94,
        projection_note,
        transform=ax_z.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.26",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.8,
            alpha=0.86,
        ),
    )

    fig.suptitle(
        r"Figure 03 — DK-RD2 Coupling Structure and Cosmological Projection"
        "\n"
        rf"$\alpha_{{\rm DK}}={alpha_DK:.4f}$ from official global DESI BAO + CC calibration",
        fontsize=title_fs-2,
        y=0.975,
    )

    fig.text(
        0.5,
        0.0080,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=footer_fs,
        color=DK_RD2_color,
    )

    fig.subplots_adjust(
        top=0.875,
        bottom=0.11,
        left=0.080,
        right=0.930,
        hspace=0.38,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    # ------------------------------------------------------------
    # Stats CSV.
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "DK-RD2 coupling structure and cosmological projection",
            "alpha_DK": float(alpha_DK),
            "rd_mpc": float(rd_DK),
            "T_min_K": float(T_vals.min()),
            "T_max_K": float(T_vals.max()),
            "z_min_projection": float(z_cosmo.min()),
            "z_max_projection": float(z_cosmo.max()),
            "beta_v_over_c": float(beta),
            "Gab_max_over_G0_at_T0_beta1": float(Gab_max_over_G0[0]),
            "Gab_eff_over_G0_at_T0_beta1": float(Gab_eff_over_G0[0]),
            "delta_G_eff_over_G0_at_z0_beta1": float(delta_G_eff_cosmo_over_G0[0]),
            "delta_G_eff_over_G0_at_z3_beta1": float(delta_G_eff_cosmo_over_G0[-1]),
            "projected_excess_fraction": float(alpha_DK),
            "thermal_projection_z0": float(thermal_projection[0]),
            "thermal_projection_z3": float(thermal_projection[-1]),
            "fit_source": "alpha_DK and r_d are obtained from the global BAO + CC calibration.",
            "fit_mode": "COUPLING_STRUCTURE_PLUS_COSMOLOGICAL_PROJECTION",
        }
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=3,
        fit_mode="COUPLING_STRUCTURE_PLUS_COSMOLOGICAL_PROJECTION",
        index=False,
    )

    return file_fig, file_table, file_stats

def generate_figure04(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    z_min: float = 0.001,
    # Limit the diagnostic reconstruction to z ≤ 2.0.
    # Beyond this range Ω_DE_eff(z) approaches zero and the
    # reconstructed FRW equation-of-state becomes ill-defined.
    # This is a limitation of the diagnostic mapping, not of DK-RD2.
    z_max = 2.0,
    n_z: int = 900,
    Core_H0: float | None = None,
    Core_Omega_m_LCDM: float | None = None,
    Omega_L_LCDM: float | None = None,
    Omega_L_DK: float | None = None,
    v_model=None,
    T_model=None,
):
    """
    # ============================================================
    # Figure 04 — Explicit Diagnostic Reconstruction of w(z)
    # ============================================================
    #
    # Purpose
    # -------
    # Reconstruct the effective equation of state w(z) that would be
    # inferred if the DK-RD2 calibrated expansion history were interpreted
    # through a standard FRW / ΛCDM matter reference.
    #
    # Important
    # ---------
    # This is a diagnostic reconstruction only.
    # It does not imply that DK-RD2 contains a fundamental dark-energy fluid.
    #
    # Definitions
    # -----------
    #     Omega_DE_eff(z) = E^2(z) - Omega_m,ref (1+z)^3
    #
    #     w(z) = -1 + (1+z)/3 * d ln[Omega_DE_eff(z)] / dz
    #
    # with:
    #     Omega_m,ref = Omega_m,LCDM = 0.315
    #
    # alpha_DK is not hardcoded.
    # It is obtained from the official global
    # DESI BAO + Cosmic Chronometer calibration
    # ============================================================
    """

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure04 requires desi_bao_mean_path and desi_bao_cov_path.")

    if Core_H0 is None:
        Core_H0 = float(Core_Hubble_H0)

    if Core_Omega_m_LCDM is None:
        Core_Omega_m_LCDM = float(Core_OMEGA_M_LCDM)

    if Omega_L_LCDM is None:
        Omega_L_LCDM = float(Core_OMEGA_L_LCDM)

    # ------------------------------------------------------------
    # Global calibrated alpha_DK access from main()
    alpha_dk_best = alpha_dk_best_global
    rd_dk_best = rd_dk_best_global
    # ------------------------------------------------------------

    file_fig = generate_evidence("image", 4)
    file_table = generate_evidence("table", 4)
    file_stats = file_table.replace(".csv", "_stats.csv")

    # ------------------------------------------------------------
    # Redshift grid.
    # ------------------------------------------------------------
    z = np.linspace(float(z_min), float(z_max), int(n_z))

    # ------------------------------------------------------------
    # ΛCDM reference.
    # ------------------------------------------------------------
    E_lcdm = E_LCDM(
        z,
        Core_Omega_m=Core_Omega_m_LCDM,
        Core_Omega_L=Omega_L_LCDM,
    )

    Omega_DE_lcdm = Omega_DE_eff_from_E(
        z,
        E_lcdm,
        Omega_m_ref=Core_Omega_m_LCDM,
    )

    w_lcdm = -1.0 * np.ones_like(z)

    # ------------------------------------------------------------
    # DK-RD2 calibrated projected expansion history.
    #
    # The core DK-RD2 background gives the raw thermodynamic envelope.
    # The observable large-scale projection uses the globally calibrated
    # alpha_DK:
    #
    #     E_DK,eff^2 = 1 + alpha_DK * (E_DK,raw^2 - 1)
    #
    # ------------------------------------------------------------
    E_dk_raw = E_Relativistic(
        z,
        Core_Omega_m=None,
        Omega_L_value=Omega_L_DK,
        v_model=v_model,
        T_model=T_model,
    )

    E2_dk_eff = 1.0 + alpha_dk_best * (E_dk_raw ** 2 - 1.0)
    E_dk_eff = np.sqrt(np.clip(E2_dk_eff, 1e-300, None))

    # ------------------------------------------------------------
    # Explicit diagnostic reconstruction.
    # ------------------------------------------------------------
    Omega_DE_dk = Omega_DE_eff_from_E(
        z,
        E_dk_eff,
        Omega_m_ref=Core_Omega_m_LCDM,
    )

    w_dk = w_eff_from_E(
        z,
        E_dk_eff,
        Omega_m_ref=Core_Omega_m_LCDM,
    )

    delta_w = w_dk - w_lcdm

    finite_w = np.isfinite(w_dk)
    finite_omega = np.isfinite(Omega_DE_dk)

    # ------------------------------------------------------------
    # Helper for representative values.
    # ------------------------------------------------------------
    def interp_safe(x_grid, y_grid, x_value):
        mask = np.isfinite(x_grid) & np.isfinite(y_grid)
        if np.count_nonzero(mask) < 2:
            return np.nan
        return float(np.interp(x_value, x_grid[mask], y_grid[mask]))

    w0 = interp_safe(z, w_dk, 0.01)
    w05 = interp_safe(z, w_dk, 0.5)
    w1 = interp_safe(z, w_dk, 1.0)
    w2 = interp_safe(z, w_dk, 2.0)

    omega0 = interp_safe(z, Omega_DE_dk, 0.01)
    omega05 = interp_safe(z, Omega_DE_dk, 0.5)
    omega1 = interp_safe(z, Omega_DE_dk, 1.0)
    omega2 = interp_safe(z, Omega_DE_dk, 2.0)

    # ------------------------------------------------------------
    # Save evidence table.
    # ------------------------------------------------------------
    table = pd.DataFrame({
        "z": z,
        "alpha_DK": np.full_like(z, alpha_dk_best),
        "rd_DK_Mpc": np.full_like(z, rd_dk_best),
        "Omega_m_ref_LCDM": np.full_like(z, Core_Omega_m_LCDM),

        "E_LCDM": E_lcdm,
        "E_DK_raw_from_Gab_max": E_dk_raw,
        "E_DK_eff_projected_by_alpha_DK": E_dk_eff,

        "Omega_DE_eff_LCDM": Omega_DE_lcdm,
        "Omega_DE_eff_DKRD2_diagnostic": Omega_DE_dk,

        "w_LCDM": w_lcdm,
        "w_DKRD2_reconstructed": w_dk,
        "delta_w_DK_minus_LCDM": delta_w,
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="FIG04_EXPLICIT_W_RECONSTRUCTION",
        figure_id=4,
        strict=False,
        index=False,
        meta={
            "alpha_DK": alpha_dk_best,
            "rd_DK_Mpc": rd_dk_best,
            "Omega_m_ref_LCDM": Core_Omega_m_LCDM,
            "z_range": f"{z_min} to {z_max}",
            "alpha_source": "Global DESI BAO + Cosmic Chronometer calibration",
            "coupling_interpretation": (
                "Gab_max is the full thermodynamic-relativistic coupling envelope; "
                "Gab_eff applies alpha_DK as the observable large-scale geometric projection."
            ),
            "expansion_projection": (
                "E_DK_eff^2 = 1 + alpha_DK * (E_DK_raw^2 - 1)."
            ),
            "interpretation": (
                "Diagnostic w(z) reconstruction using a ΛCDM matter reference; "
                "not a fundamental dark-energy component in DK-RD2."
            ),
        },
    )

    # ------------------------------------------------------------
    # Save stats.
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "ΛCDM reference",
            "alpha_DK": np.nan,
            "rd_mpc": np.nan,
            "Omega_m_ref": float(Core_Omega_m_LCDM),
            "w_z_0p01": -1.0,
            "w_z_0p5": -1.0,
            "w_z_1": -1.0,
            "w_z_2": -1.0,
            "Omega_DE_z_0p01": float(Omega_L_LCDM),
            "Omega_DE_z_0p5": float(Omega_L_LCDM),
            "Omega_DE_z_1": float(Omega_L_LCDM),
            "Omega_DE_z_2": float(Omega_L_LCDM),
            "w_min": -1.0,
            "w_max": -1.0,
            "mean_delta_w_vs_minus1": 0.0,
            "N_finite": int(len(z)),
            "fit_mode": "reference_constant",
        },
        {
            "model": "DK-RD2 explicit diagnostic w(z)",
            "alpha_DK": float(alpha_dk_best),
            "rd_mpc": float(rd_dk_best),
            "Omega_m_ref": float(Core_Omega_m_LCDM),
            "w_z_0p01": w0,
            "w_z_0p5": w05,
            "w_z_1": w1,
            "w_z_2": w2,
            "Omega_DE_z_0p01": omega0,
            "Omega_DE_z_0p5": omega05,
            "Omega_DE_z_1": omega1,
            "Omega_DE_z_2": omega2,
            "w_min": float(np.nanmin(w_dk)) if np.any(finite_w) else np.nan,
            "w_max": float(np.nanmax(w_dk)) if np.any(finite_w) else np.nan,
            "mean_delta_w_vs_minus1": float(np.nanmean(w_dk + 1.0)) if np.any(finite_w) else np.nan,
            "Omega_DE_min": float(np.nanmin(Omega_DE_dk)) if np.any(finite_omega) else np.nan,
            "Omega_DE_max": float(np.nanmax(Omega_DE_dk)) if np.any(finite_omega) else np.nan,
            "N_finite": int(np.count_nonzero(finite_w)),
            "fit_mode": "diagnostic_LCDM_reference",
        },
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=4,
        fit_mode="EXPLICIT_W_RECONSTRUCTION",
        index=False,
    )

    # ------------------------------------------------------------
    # Plot configuration.
    # ------------------------------------------------------------
    fig, (ax_omega, ax_w, ax_delta) = plt.subplots(
        3,
        1,
        figsize=(11.8, 10.4),
        dpi=130,
        sharex=True,
        gridspec_kw={
            "height_ratios": [1.15, 2.25, 1.05],
            "hspace": 0.13,
        },
    )

    title_fs = 14
    label_fs = 9.5
    tick_fs = 8.2
    legend_fs = 7.4
    box_fs = 7.0
    footer_fs = 5.4

    # ------------------------------------------------------------
    # Panel A — Omega_DE_eff(z).
    # ------------------------------------------------------------
    ax_omega.plot(
        z,
        Omega_DE_lcdm,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=1.8,
        label=r"$\Lambda$CDM: $\Omega_{\Lambda}=\mathrm{const.}$",
    )

    ax_omega.plot(
        z,
        Omega_DE_dk,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=r"DK-RD2 diagnostic $\Omega_{\rm DE}^{eff}(z)$",
    )

    ax_omega.axhline(
        0.0,
        color="black",
        linewidth=0.8,
        alpha=0.45,
    )

    ax_omega.set_ylabel(r"$\Omega_{\rm DE}^{eff}(z)$", fontsize=label_fs)
    ax_omega.set_title(
        rf"Figure 04 — Explicit Diagnostic Reconstruction of $w(z)$"
        "\n"
        rf"$\alpha_{{\rm DK}}={alpha_dk_best:.4f}$, "
        rf"$r_d={rd_dk_best:.2f}\,\mathrm{{Mpc}}$ "
        rf"(global DESI BAO + CC calibration)",
        fontsize=title_fs,
        pad=8,
    )

    ax_omega.tick_params(axis="both", labelsize=tick_fs)
    ax_omega.grid(alpha=0.28)
    ax_omega.legend(fontsize=legend_fs, loc="best", framealpha=0.88)

    # ------------------------------------------------------------
    # Panel B — Explicit reconstructed w(z).
    # ------------------------------------------------------------
    ax_w.plot(
        z,
        w_lcdm,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=2.0,
        label=r"$\Lambda$CDM: $w=-1$",
    )

    ax_w.plot(
        z,
        w_dk,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.4,
        label=(
            rf"DK-RD2 reconstructed $w(z)$ "
            rf"$(\Omega_{{m,\rm ref}}={Core_Omega_m_LCDM:.3f})$"
        ),
    )

    ax_w.axhline(
        -1.0,
        color="black",
        linewidth=0.9,
        alpha=0.35,
    )

    ax_w.set_ylabel(r"Reconstructed $w(z)$", fontsize=label_fs)
    ax_w.tick_params(axis="both", labelsize=tick_fs)
    ax_w.grid(alpha=0.28)
    ax_w.legend(fontsize=legend_fs, loc="upper right", framealpha=0.88)

    eq_text = (
        r"$E_{\rm DK,eff}^2=1+\alpha_{\rm DK}(E_{\rm DK,raw}^2-1)$" "\n"
        r"$\Omega_{\rm DE}^{eff}=E_{\rm DK,eff}^2-\Omega_{m,\rm ref}(1+z)^3$" "\n"
        r"$w(z)=-1+\frac{1+z}{3}\frac{d\ln\Omega_{\rm DE}^{eff}}{dz}$"
    )

    ax_w.text(
        0.035,
        0.94,
        eq_text,
        transform=ax_w.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs,
        linespacing=1.08,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="black",
            linewidth=0.8,
            alpha=0.88,
        ),
    )

    result_text = (
        r"DK-RD2 reconstructed:" "\n"
        rf"$w(0.01)={w0:.3f}$" "\n"
        rf"$w(0.5)={w05:.3f}$" "\n"
        rf"$w(1.0)={w1:.3f}$"
    )

    ax_w.text(
        0.745,
        0.43,
        result_text,
        transform=ax_w.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs,
        linespacing=1.15,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.8,
            alpha=0.86,
        ),
    )

    note_text = (
        "Diagnostic reconstruction only:\n"
        rf"$\Omega_{{m,\rm ref}}={Core_Omega_m_LCDM:.3f}$ is the ΛCDM matter reference."
        "\nThis is not a fundamental dark-energy component in DK-RD2."
    )

    ax_w.text(
        0.035,
        0.055,
        note_text,
        transform=ax_w.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.4,
        linespacing=1.08,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.86,
        ),
    )

    highz_note = (
        "Diagnostic FRW mapping only.\n"
        "The high-z divergence is not a physical singularity of DK-RD2."
    )

    ax_w.text(
        0.035,
        0.235,
        highz_note,
        transform=ax_w.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.1,
        linespacing=1.08,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.86,
        ),
    )

    # ------------------------------------------------------------
    # Panel C — delta w.
    # ------------------------------------------------------------
    ax_delta.plot(
        z,
        delta_w,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=2.0,
        label=r"$\Delta w(z)=w_{\rm DK}(z)+1$",
    )

    ax_delta.axhline(
        0.0,
        color="black",
        linewidth=0.9,
        alpha=0.70,
    )

    ax_delta.set_xlabel("Redshift z", fontsize=label_fs)
    ax_delta.set_ylabel(r"$\Delta w$", fontsize=label_fs)
    ax_delta.tick_params(axis="both", labelsize=tick_fs)
    ax_delta.grid(alpha=0.28)
    ax_delta.legend(fontsize=legend_fs, loc="best", framealpha=0.88)

    # ------------------------------------------------------------
    # Footer.
    # ------------------------------------------------------------
    fig.text(
        0.5,
        0.018,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=footer_fs,
        color=DK_RD2_color,
    )

    fig.subplots_adjust(
        top=0.900,
        bottom=0.105,
        left=0.080,
        right=0.965,
        hspace=0.17,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats

def generate_figure05(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    z_min: float = 0.001,

    # Limit the diagnostic reconstruction to z ≤ 2.0.
    # Beyond this range Ω_DE_eff(z) approaches zero and the
    # reconstructed FRW equation-of-state becomes ill-defined.
    # This is a limitation of the diagnostic mapping, not of DK-RD2.
    z_max = 2.0,
    n_z: int = 500,
    Core_H0: float | None = None,
    Core_Omega_m_LCDM: float | None = None,
    Omega_L_LCDM: float | None = None,
    sigma8_0: float | None = None,
    a_min: float = 1.0e-3,
    n_grid: int = 5000,
):
    """
    Figure 05 — Growth of Cosmic Structure

    Conservative paper version:
    - Main plotted DK-RD2 growth curve uses background-only growth.
    - Modified-Poisson mu_DK version is computed only as diagnostic.
    - alpha_DK is obtained from the global DESI BAO + CC calibration.
    """

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure05 requires desi_bao_mean_path and desi_bao_cov_path.")

    if Core_H0 is None:
        Core_H0 = float(Core_Hubble_H0)

    if Core_Omega_m_LCDM is None:
        Core_Omega_m_LCDM = float(Core_OMEGA_M_LCDM)

    if Omega_L_LCDM is None:
        Omega_L_LCDM = float(Core_OMEGA_L_LCDM)

    if sigma8_0 is None:
        sigma8_0 = float(Core_SIGMA8_PLANCK_2018)

    # ------------------------------------------------------------
    # Global calibrated alpha_DK access from main()
    alpha_dk_best = alpha_dk_best_global
    rd_dk_best = rd_dk_best_global
    # ------------------------------------------------------------

    file_fig = generate_evidence("image", 5)
    file_table = generate_evidence("table", 5)
    file_stats = file_table.replace(".csv", "_stats.csv")

    z_eval = np.linspace(float(z_min), float(z_max), int(n_z))

    # ------------------------------------------------------------
    # Growth solver.
    # ------------------------------------------------------------
    def solve_growth_from_background(
        z_output,
        E_of_z,
        source_S_of_z,
        *,
        a_min_local: float = 1.0e-3,
        n_grid_local: int = 5000,
    ):
        z_output = np.asarray(z_output, dtype=float)

        x_i = np.log(float(a_min_local))
        x_f = 0.0

        x_grid = np.linspace(x_i, x_f, int(n_grid_local))
        a_grid = np.exp(x_grid)
        z_grid = (1.0 / a_grid) - 1.0

        E_grid = np.asarray(E_of_z(z_grid), dtype=float)
        E_grid = np.clip(E_grid, 1.0e-300, None)

        lnE = np.log(E_grid)
        dlnH_dx = np.gradient(lnE, x_grid)

        S_grid = np.asarray(source_S_of_z(z_grid), dtype=float)
        S_grid = np.where(np.isfinite(S_grid), S_grid, 0.0)

        D = np.zeros_like(x_grid, dtype=float)
        Dp = np.zeros_like(x_grid, dtype=float)

        D[0] = float(a_min_local)
        Dp[0] = D[0]

        dx = float(x_grid[1] - x_grid[0])

        def rhs(i, y0, y1):
            A = 2.0 + float(dlnH_dx[i])
            B = 1.5 * float(S_grid[i])
            return y1, -A * y1 + B * y0

        for i in range(len(x_grid) - 1):
            y0 = D[i]
            y1 = Dp[i]

            k1_0, k1_1 = rhs(i, y0, y1)
            k2_0, k2_1 = rhs(i, y0 + 0.5 * dx * k1_0, y1 + 0.5 * dx * k1_1)
            k3_0, k3_1 = rhs(i, y0 + 0.5 * dx * k2_0, y1 + 0.5 * dx * k2_1)
            k4_0, k4_1 = rhs(i, y0 + dx * k3_0, y1 + dx * k3_1)

            D[i + 1] = y0 + (dx / 6.0) * (k1_0 + 2.0 * k2_0 + 2.0 * k3_0 + k4_0)
            Dp[i + 1] = y1 + (dx / 6.0) * (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1)

        D = D / np.clip(D[-1], 1.0e-300, None)

        z_rev = z_grid[::-1]
        D_rev = D[::-1]

        D_out = np.interp(z_output, z_rev, D_rev)

        a_out = 1.0 / (1.0 + z_output)
        ln_a_out = np.log(np.clip(a_out, 1.0e-300, None))
        ln_D_out = np.log(np.clip(D_out, 1.0e-300, None))

        idx = np.argsort(ln_a_out)
        f_sorted = np.gradient(ln_D_out[idx], ln_a_out[idx])

        f_out = np.empty_like(f_sorted)
        f_out[idx] = f_sorted

        fsigma8_out = f_out * float(sigma8_0) * D_out

        return D_out, f_out, fsigma8_out

    # ------------------------------------------------------------
    # LCDM reference.
    # ------------------------------------------------------------
    def E_lcdm_func(z_in):
        return E_LCDM(
            z_in,
            Core_Omega_m=Core_Omega_m_LCDM,
            Core_Omega_L=Omega_L_LCDM,
        )

    def S_lcdm_func(z_in):
        z_arr = np.asarray(z_in, dtype=float)
        E_vals = np.asarray(E_lcdm_func(z_arr), dtype=float)

        return (
            float(Core_Omega_m_LCDM)
            * (1.0 + z_arr) ** 3
            / np.clip(E_vals ** 2, 1.0e-300, None)
        )

    D_lcdm, f_lcdm, fs8_lcdm = solve_growth_from_background(
        z_eval,
        E_lcdm_func,
        S_lcdm_func,
        a_min_local=a_min,
        n_grid_local=n_grid,
    )

    # ------------------------------------------------------------
    # DK-RD2 projected background.
    # ------------------------------------------------------------
    def E_dk_raw_func(z_in):
        return E_Relativistic(
            z_in,
            Core_Omega_m=None,
            Omega_L_value=None,
            v_model=None,
            T_model=None,
        )

    def E_dk_eff_func(z_in):
        z_arr = np.asarray(z_in, dtype=float)
        E_raw = np.asarray(E_dk_raw_func(z_arr), dtype=float)
        E2_eff = 1.0 + alpha_dk_best * (E_raw ** 2 - 1.0)

        return np.sqrt(np.clip(E2_eff, 1.0e-300, None))

    def mu_dk_projected_func(z_in):
        z_arr = np.asarray(z_in, dtype=float)

        return 1.0 + alpha_dk_best / (1.0 + z_arr)

    # Conservative paper source: background-only.
    def S_dk_background_only_func(z_in):
        z_arr = np.asarray(z_in, dtype=float)
        E_vals = np.asarray(E_dk_eff_func(z_arr), dtype=float)

        return (
            float(Core_Omega_m_LCDM)
            * (1.0 + z_arr) ** 3
            / np.clip(E_vals ** 2, 1.0e-300, None)
        )

    # Diagnostic only: background plus modified-Poisson response.
    # This projected mu_DK is not used in the paper curve.
    # The published Figure05 uses the background-only source term.
    def S_dk_mu_func(z_in):
        z_arr = np.asarray(z_in, dtype=float)
        S_bg = S_dk_background_only_func(z_arr)
        mu_vals = mu_dk_projected_func(z_arr)

        return S_bg * mu_vals

    # ------------------------------------------------------------
    # Solve DK-RD2 growth.
    # ------------------------------------------------------------
    D_dk_bg, f_dk_bg, fs8_dk_bg = solve_growth_from_background(
        z_eval,
        E_dk_eff_func,
        S_dk_background_only_func,
        a_min_local=a_min,
        n_grid_local=n_grid,
    )

    D_dk_mu, f_dk_mu, fs8_dk_mu = solve_growth_from_background(
        z_eval,
        E_dk_eff_func,
        S_dk_mu_func,
        a_min_local=a_min,
        n_grid_local=n_grid,
    )

    # ------------------------------------------------------------
    # Paper curve selection: conservative background-only.
    # ------------------------------------------------------------
    D_dk_plot = D_dk_bg
    f_dk_plot = f_dk_bg
    fs8_dk_plot = fs8_dk_bg
    S_dk_plot = S_dk_background_only_func(z_eval)

    # Diagnostic modified-Poisson curve.
    S_dk_mu = S_dk_mu_func(z_eval)

    mu_dk = mu_dk_projected_func(z_eval)

    S_lcdm = S_lcdm_func(z_eval)
    S_dk_bg = S_dk_background_only_func(z_eval)

    E_lcdm_eval = E_lcdm_func(z_eval)
    E_dk_eff_eval = E_dk_eff_func(z_eval)

    double_count_check = S_dk_mu / np.clip(S_dk_bg * mu_dk, 1.0e-300, None)
    double_count_error = np.nanmax(np.abs(double_count_check - 1.0))

    print("============================================================")
    print("Figure 05 DK-RD2 growth-source consistency check")
    print(f"max |S_DK_mu / (S_DK_bg * mu_DK) - 1| = {double_count_error:.3e}")
    print("If this is ~0, the diagnostic modified-Poisson source is algebraically consistent.")
    print("The plotted paper curve uses the conservative background-only source.")
    print("============================================================")

    # ------------------------------------------------------------
    # Diagnostics based on paper curve.
    # ------------------------------------------------------------
    delta_fs8 = fs8_dk_plot - fs8_lcdm
    frac_delta_fs8 = delta_fs8 / np.clip(fs8_lcdm, 1.0e-300, None)

    delta_D = D_dk_plot - D_lcdm
    frac_delta_D = delta_D / np.clip(D_lcdm, 1.0e-300, None)

    def interp_safe(x_grid, y_grid, x_value):
        mask = np.isfinite(x_grid) & np.isfinite(y_grid)
        if np.count_nonzero(mask) < 2:
            return np.nan
        return float(np.interp(x_value, x_grid[mask], y_grid[mask]))

    fs8_lcdm_0 = interp_safe(z_eval, fs8_lcdm, 0.0)
    fs8_lcdm_05 = interp_safe(z_eval, fs8_lcdm, 0.5)
    fs8_lcdm_1 = interp_safe(z_eval, fs8_lcdm, 1.0)

    fs8_dk_0 = interp_safe(z_eval, fs8_dk_plot, 0.0)
    fs8_dk_05 = interp_safe(z_eval, fs8_dk_plot, 0.5)
    fs8_dk_1 = interp_safe(z_eval, fs8_dk_plot, 1.0)

    frac_05 = interp_safe(z_eval, frac_delta_fs8, 0.5)
    frac_1 = interp_safe(z_eval, frac_delta_fs8, 1.0)
    frac_2 = interp_safe(z_eval, frac_delta_fs8, 2.0)

    # ------------------------------------------------------------
    # Save evidence table.
    # ------------------------------------------------------------
    table = pd.DataFrame({
        "z": z_eval,
        "alpha_DK": np.full_like(z_eval, alpha_dk_best),
        "rd_DK_Mpc": np.full_like(z_eval, rd_dk_best),
        "Omega_m_ref_LCDM": np.full_like(z_eval, Core_Omega_m_LCDM),
        "sigma8_0": np.full_like(z_eval, sigma8_0),

        "E_LCDM": E_lcdm_eval,
        "E_DK_eff_projected": E_dk_eff_eval,
        "mu_DK_projected": mu_dk,

        "source_S_LCDM": S_lcdm,
        "source_S_DKRD2_background_only": S_dk_bg,
        "source_S_DKRD2_modified_poisson": S_dk_mu,

        "D_LCDM": D_lcdm,
        "D_DKRD2_background_only": D_dk_bg,
        "D_DKRD2_modified_poisson": D_dk_mu,
        "D_DKRD2_paper_curve": D_dk_plot,

        "f_LCDM": f_lcdm,
        "f_DKRD2_background_only": f_dk_bg,
        "f_DKRD2_modified_poisson": f_dk_mu,
        "f_DKRD2_paper_curve": f_dk_plot,

        "fsigma8_LCDM": fs8_lcdm,
        "fsigma8_DKRD2_background_only": fs8_dk_bg,
        "fsigma8_DKRD2_modified_poisson": fs8_dk_mu,
        "fsigma8_DKRD2_paper_curve": fs8_dk_plot,

        "delta_fsigma8_DK_minus_LCDM": delta_fs8,
        "fractional_delta_fsigma8": frac_delta_fs8,

        "delta_D_DK_minus_LCDM": delta_D,
        "fractional_delta_D": frac_delta_D,

        "double_count_check_SDK_mu_over_Sbg_mu": double_count_check,
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="FIG05_GROWTH_OF_COSMIC_STRUCTURE",
        figure_id=5,
        strict=False,
        index=False,
        meta={
            "alpha_DK": alpha_dk_best,
            "rd_DK_Mpc": rd_dk_best,
            "Omega_m_ref_LCDM": Core_Omega_m_LCDM,
            "sigma8_0": sigma8_0,
            "z_range": f"{z_min} to {z_max}",
            "alpha_source": "Global DESI BAO + Cosmic Chronometer calibration",
            "growth_equation": "D_xx + [2 + dlnH/dx]D_x - 3/2 S(x)D = 0",
            "dk_background_projection": "E_DK_eff^2 = 1 + alpha_DK * (E_DK_raw^2 - 1)",
            "paper_curve": "DK-RD2 background-only growth",
            "diagnostic_curve": "DK-RD2 modified-Poisson mu_DK curve stored in CSV only",
            "double_count_check_max_error": double_count_error,
            "interpretation": (
                "Conservative linear growth diagnostic; no additional dark-energy sector "
                "or phenomenological growth parameter is fitted."
            ),
        },
    )

    # ------------------------------------------------------------
    # Save stats.
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "LambdaCDM growth reference",
            "alpha_DK": np.nan,
            "rd_mpc": np.nan,
            "Omega_m_ref": float(Core_Omega_m_LCDM),
            "sigma8_0": float(sigma8_0),
            "fsigma8_z0": fs8_lcdm_0,
            "fsigma8_z0p5": fs8_lcdm_05,
            "fsigma8_z1": fs8_lcdm_1,
            "fit_mode": "LCDM_LINEAR_GROWTH_REFERENCE",
        },
        {
            "model": "DK-RD2 background-only linear growth",
            "alpha_DK": float(alpha_dk_best),
            "rd_mpc": float(rd_dk_best),
            "Omega_m_ref": float(Core_Omega_m_LCDM),
            "sigma8_0": float(sigma8_0),
            "fsigma8_z0": fs8_dk_0,
            "fsigma8_z0p5": fs8_dk_05,
            "fsigma8_z1": fs8_dk_1,
            "fractional_delta_fsigma8_z0p5": frac_05,
            "fractional_delta_fsigma8_z1": frac_1,
            "fractional_delta_fsigma8_z2": frac_2,
            "max_abs_fractional_delta_fsigma8": (
                float(np.nanmax(np.abs(frac_delta_fs8)))
                if np.any(np.isfinite(frac_delta_fs8))
                else np.nan
            ),
            "mean_fractional_delta_fsigma8": (
                float(np.nanmean(frac_delta_fs8))
                if np.any(np.isfinite(frac_delta_fs8))
                else np.nan
            ),
            "double_count_check_max_error": double_count_error,
            "fit_mode": "DKRD2_BACKGROUND_ONLY_LINEAR_GROWTH_DIAGNOSTIC",
        },
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=5,
        fit_mode="GROWTH_OF_COSMIC_STRUCTURE",
        index=False,
    )

    # ------------------------------------------------------------
    # Plot.
    # ------------------------------------------------------------
    fig, (ax_fs8, ax_frac, ax_source) = plt.subplots(
        3,
        1,
        figsize=(11.8, 10.2),
        dpi=130,
        sharex=True,
        gridspec_kw={
            "height_ratios": [2.05, 1.25, 1.15],
            "hspace": 0.13,
        },
    )

    title_fs = 14
    label_fs = 9.5
    tick_fs = 8.2
    legend_fs = 7.4
    box_fs = 6.8
    footer_fs = 5.4

    # Panel A — f sigma8.
    ax_fs8.plot(
        z_eval,
        fs8_lcdm,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=2.0,
        label=r"$\Lambda$CDM reference $f\sigma_8(z)$",
    )

    ax_fs8.plot(
        z_eval,
        fs8_dk_plot,
        color=DK_RD2_color,
        linestyle=":",
        linewidth=2.0,
        alpha=0.95,
        label=r"DK-RD2 background-only growth",
    )

    ax_fs8.set_ylabel(r"$f\sigma_8(z)$", fontsize=label_fs)
    ax_fs8.set_title(
        rf"Figure 05 — Growth of Cosmic Structure"
        "\n"
        rf"$\alpha_{{\rm DK}}={alpha_dk_best:.4f}$, "
        rf"$r_d={rd_dk_best:.2f}\,\mathrm{{Mpc}}$ "
        rf"(global DESI BAO + CC calibration)",
        fontsize=title_fs,
        pad=8,
    )

    ax_fs8.tick_params(axis="both", labelsize=tick_fs)
    ax_fs8.grid(alpha=0.28)
    ax_fs8.legend(fontsize=legend_fs, loc="best", framealpha=0.88)

    result_text = (
        r"Background-only growth:" "\n"
        r"Representative values:" "\n"
        rf"$f\sigma_8^{{\mathrm{{DK,bg}}}}(0.5)={fs8_dk_05:.3f}$" "\n"
        rf"$f\sigma_8^{{\Lambda\mathrm{{CDM}}}}(0.5)={fs8_lcdm_05:.3f}$" "\n"
        rf"$\Delta_{{\mathrm{{frac}}}}(0.5)={100.0 * frac_05:.2f}\%$" "\n"
        rf"$\Delta_{{\mathrm{{frac}}}}(1.0)={100.0 * frac_1:.2f}\%$"
    )

    ax_fs8.text(
        0.035,
        0.52,
        result_text,
        transform=ax_fs8.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs,
        linespacing=1.08,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.8,
            alpha=0.88,
        ),
    )

    # Panel B — fractional deviation.
    ax_frac.plot(
        z_eval,
        100.0 * frac_delta_fs8,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=2.0,
        label=r"$100\times[f\sigma_8^{DK,bg}/f\sigma_8^{\Lambda CDM}-1]$",
    )

    ax_frac.axhline(
        0.0,
        color="black",
        linewidth=0.9,
        alpha=0.70,
    )

    ax_frac.set_ylabel(r"$\Delta f\sigma_8$ [%]", fontsize=label_fs)
    ax_frac.tick_params(axis="both", labelsize=tick_fs)
    ax_frac.grid(alpha=0.28)
    ax_frac.legend(fontsize=legend_fs, loc="best", framealpha=0.88)

    note_text = (
        "Growth diagnostic:\n"
        "The DK-RD2 curve uses the projected background from Figures 01–04.\n"
        "No independent growth parameter is fitted."
    )

    ax_frac.text(
        0.165,
        0.50,
        note_text,
        transform=ax_frac.transAxes,
        ha="left",
        va="top",
        fontsize=6.3,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.86,
        ),
    )

    # Panel C — growth source term.
    ax_source.plot(
        z_eval,
        S_lcdm,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=1.9,
        label=r"$S_{\Lambda CDM}(z)$",
    )

    ax_source.plot(
        z_eval,
        S_dk_plot,
        color=DK_RD2_color,
        linestyle=":",
        linewidth=1.9,
        alpha=0.95,
        label=r"$S_{\rm DK,bg}(z)$",
    )

    ax_source.set_xlabel("Redshift z", fontsize=label_fs)
    ax_source.set_ylabel(r"Growth source $S(z)$", fontsize=label_fs)
    ax_source.tick_params(axis="both", labelsize=tick_fs)
    ax_source.grid(alpha=0.28)
    ax_source.legend(fontsize=legend_fs, loc="best", framealpha=0.88)

    source_text = (
        r"$S_{\rm DK,bg}(z)=\Omega_{m,\rm ref}(1+z)^3/E_{\rm DK,eff}^2$" "\n"
        r"Modified-Poisson $\mu_{\rm DK}$ test stored only as diagnostic."
    )

    ax_source.text(
        0.035,
        0.92,
        source_text,
        transform=ax_source.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.86,
        ),
    )

    # Footer.
    fig.text(
        0.5,
        0.018,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=footer_fs,
        color=DK_RD2_color,
    )

    fig.subplots_adjust(
        top=0.900,
        bottom=0.105,
        left=0.080,
        right=0.965,
        hspace=0.17,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats

def generate_figure06(
    sidm_notebook_path: str = "data.nb",
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    n_r: int = 650,
):
    """
    # ============================================================
    # Figure 06 — Compact-Halo Prediction and Effective Gravitational Amplification
    # ============================================================
    #
    # Purpose
    # -------
    # Test whether compact-halo phenomenology commonly interpreted through
    # dense dark-matter structures can also be represented as an effective
    # gravitational amplification within the DK-RD2 framework.
    #
    # This figure uses published core-collapsed SIDM halo reconstructions as
    # an external observational benchmark and compares them with a DK-RD2
    # effective-gravity interpretation.
    #
    # Scientific role in this paper
    # -----------------------------
    # Figures 01–05 test DK-RD2 through cosmological observables:
    #
    #   • BAO + Cosmic Chronometer calibration
    #   • alpha_DK likelihood profiling
    #   • thermodynamic coupling projection
    #   • effective w(z) reconstruction
    #   • large-scale structure growth
    #
    # Figure 06 extends the falsification program to compact galactic and
    # sub-galactic systems.
    #
    # Important interpretation
    # ------------------------
    # This figure does NOT claim that SIDM or ΛCDM are wrong.
    # It evaluates whether part of the same compact-halo phenomenology can be
    # reproduced through an effective thermodynamic-relativistic gravitational
    # amplification:
    #
    #       mu_eff(r) = G_eff(r) / G0
    #
    # In standard ΛCDM / standard gravity:
    #
    #       mu_eff(r) = 1
    #
    # Therefore, any bounded radial amplification mu_eff(r) != 1 represents
    # a direct observational discriminator between constant-coupling gravity
    # and an effective-gravity interpretation.
    #
    # Global calibration
    # ------------------
    # alpha_DK is obtained from the same global DESI BAO + Cosmic Chronometer
    # calibration used in Figures 01–05.
    #
    # alpha_DK is NOT fitted to the compact-halo profiles.
    # It is stored and displayed only as the global cosmological calibration
    # for consistency across the paper.
    #
    # The compact-halo amplification profiles are phenomenological local
    # effective-gravity fits:
    #
    #       mu_DK(r) = 1 + A0 / [1 + (r / rc)^n]
    #
    # where A0, rc and n are fitted separately to each compact benchmark.
    # ============================================================
    """

    import re
    from pathlib import Path
    from scipy.interpolate import interp1d
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    file_fig = generate_evidence("image", 6)
    file_table = generate_evidence("table", 6)
    file_stats = file_table.replace(".csv", "_stats.csv")

    # ------------------------------------------------------------
    # Global calibrated alpha_DK access from main()
    alpha_dk_best = alpha_dk_best_global
    rd_dk_best = rd_dk_best_global
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Load external SIDM benchmark notebook.
    # ------------------------------------------------------------
    nb_path = Path(sidm_notebook_path)
    if not nb_path.exists():
        raise FileNotFoundError(
            f"Figure06 requires the Zenodo Mathematica notebook data.nb. "
            f"File not found: {sidm_notebook_path}"
        )

    nb_text = nb_path.read_text(encoding="utf-8", errors="ignore")

    # ------------------------------------------------------------
    # Local color aliases.
    # ------------------------------------------------------------
    dk_color = DK_RD2_color
    lcdm_color = LCDM_color
    dk_light = DK_LIGHT_color
    lcdm_light = LCDM_LIGHT_color
    error_color = ERROR_color
    secondary_color = SECONDARY_color

    cdm_color = "black"
    sidm30_color = lcdm_color
    sidm50_color = lcdm_light
    sidm100_color = secondary_color

    b1938_color = "#4a4a4a"
    fornax_color = "#8a8a8a"
    gd1_color = "0.25"
    gd1_band_color = error_color

    # ------------------------------------------------------------
    # Helper: extract Mathematica RowBox point arrays from data.nb.
    # ------------------------------------------------------------
    def _extract_nb_rowbox_pairs(var_name: str) -> np.ndarray:
        start = nb_text.find(var_name)
        if start < 0:
            raise ValueError(f"Variable {var_name} not found in {sidm_notebook_path}")

        candidate_vars = [
            "DMrhodataCDM",
            "DMrhodataS30",
            "DMrhodataS50",
            "DMrhodataS100",
            "LogLogPlot",
            "ListLogLogPlot",
        ]

        next_positions = [
            nb_text.find(v, start + 1)
            for v in candidate_vars
            if nb_text.find(v, start + 1) > start
        ]

        end = min(next_positions) if next_positions else start + 20000
        block = nb_text[start:end]

        pair_pattern = re.compile(
            r'RowBox\[\{"\{",\s*RowBox\[\{"([0-9.]+)"\s*,\s*",",\s*'
            r'(?:RowBox\[\{"([0-9.]+)"\s*,\s*"\*"\s*,\s*'
            r'RowBox\[\{"10"\s*,\s*"\^"\s*,\s*"(-?\d+)"\}\]\}\]|"([0-9.]+)")'
            r'\}\]\s*,\s*"\}"\}\]',
            re.S,
        )

        pairs = []
        for match in pair_pattern.finditer(block):
            x_val = float(match.group(1))

            if match.group(4) is not None:
                y_val = float(match.group(4))
            else:
                y_val = float(match.group(2)) * 10.0 ** int(match.group(3))

            pairs.append((x_val, y_val))

        if len(pairs) < 5:
            raise ValueError(
                f"Could not extract enough data points for {var_name}. "
                f"Check data.nb format."
            )

        arr = np.array(pairs, dtype=float)
        arr = arr[np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])]
        arr = arr[arr[:, 1] > 0.0]
        arr = arr[np.argsort(arr[:, 0])]
        return arr

    # ------------------------------------------------------------
    # 1. Extract published CDM / SIDM simulation curves.
    # ------------------------------------------------------------
    cdm_data = _extract_nb_rowbox_pairs("DMrhodataCDM")
    sidm30_data = _extract_nb_rowbox_pairs("DMrhodataS30")
    sidm50_data = _extract_nb_rowbox_pairs("DMrhodataS50")
    sidm100_data = _extract_nb_rowbox_pairs("DMrhodataS100")

    # ------------------------------------------------------------
    # 2. Observational analytic compact-halo benchmarks.
    # ------------------------------------------------------------
    r = np.logspace(np.log10(0.0029), np.log10(0.50), int(n_r))

    # B1938+666 pseudo-Jaffe-like benchmark.
    rho0_pj = 4.3e7
    rt_pj = 0.149
    rho_powell = rho0_pj * rt_pj**4 / (r**2 * (r**2 + rt_pj**2))

    # Fornax 6 Hernquist-like benchmark.
    M_fornax = 1.0e6
    a_fornax = 0.020
    rho_fornax = (M_fornax / (2.0 * np.pi)) * a_fornax / (
        r * (r + a_fornax) ** 3
    )

    # GD-1 compact perturber envelope.
    gd1_enc_params = np.array([
        [1.0000e-3, 8.7947e7],
        [1.2248e-2, 1.5250e8],
        [2.8573e-2, 9.6397e7],
        [3.1435e-2, 1.0663e7],
        [1.5736e-2, 4.6062e5],
        [1.0000e-3, 1.7580e5],
    ], dtype=float)

    gd1_profiles = []
    for a_i, M_i in gd1_enc_params:
        rho_i = (M_i / (2.0 * np.pi)) * a_i / (r * (r + a_i) ** 3)
        gd1_profiles.append(rho_i)

    gd1_profiles = np.array(gd1_profiles)
    rho_gd1_low = np.nanmin(gd1_profiles, axis=0)
    rho_gd1_high = np.nanmax(gd1_profiles, axis=0)
    rho_gd1_rep = np.sqrt(rho_gd1_low * rho_gd1_high)

    # ------------------------------------------------------------
    # 3. Interpolate CDM baseline from the external notebook.
    # ------------------------------------------------------------
    cdm_interp = interp1d(
        np.log(cdm_data[:, 0]),
        np.log(cdm_data[:, 1]),
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
    )

    rho_cdm_base = np.exp(cdm_interp(np.log(r)))
    rho_cdm_base = np.clip(rho_cdm_base, 1e-300, None)

    # ------------------------------------------------------------
    # 4. Required effective amplification relative to the CDM baseline.
    #
    # If a compact profile is interpreted as CDM baseline times an
    # effective gravitational response, then:
    #
    #       mu_required(r) = rho_target(r) / rho_CDM(r)
    #
    # In standard gravity / ΛCDM:
    #
    #       mu_eff(r) = 1
    #
    # ------------------------------------------------------------
    mu_powell_required = rho_powell / rho_cdm_base
    mu_fornax_required = rho_fornax / rho_cdm_base
    mu_gd1_rep_required = rho_gd1_rep / rho_cdm_base
    mu_gd1_low_required = rho_gd1_low / rho_cdm_base
    mu_gd1_high_required = rho_gd1_high / rho_cdm_base

    mu_lcdm_standard = np.ones_like(r)

    def mu_DK_profile(r_in, A0, r_c, n):
        """
        Phenomenological compact-system effective-gravity amplification.

            mu_DK(r) = 1 + A0 / [1 + (r/rc)^n]

        This is a local compact-system test and is not the cosmological
        alpha_DK calibration used in Figures 01–05.
        """
        r_in = np.asarray(r_in, dtype=float)
        return 1.0 + A0 / (1.0 + (r_in / r_c) ** n)

    def fit_mu_profile(target_mu, label, r_min=0.003, r_max=0.12):
        fit_mask = (
            (r >= r_min)
            & (r <= r_max)
            & np.isfinite(target_mu)
            & (target_mu > 0.0)
        )

        if np.count_nonzero(fit_mask) < 10:
            raise RuntimeError(f"Not enough valid points to fit {label}")

        def objective(theta):
            A0, r_c, n = theta

            if A0 <= 0.0 or r_c <= 0.0 or n <= 0.0:
                return 1.0e99

            mu_model = mu_DK_profile(r[fit_mask], A0, r_c, n)

            if not np.all(np.isfinite(mu_model)):
                return 1.0e99

            residual = np.log10(mu_model) - np.log10(target_mu[fit_mask])
            return float(np.mean(residual**2))

        res_fit = minimize(
            objective,
            x0=np.array([20.0, 0.010, 2.0]),
            bounds=[
                (0.01, 5000.0),
                (0.001, 0.250),
                (0.20, 10.0),
            ],
            method="L-BFGS-B",
        )

        if not res_fit.success:
            raise RuntimeError(f"DK-RD2 fit failed for {label}: {res_fit.message}")

        A0_best, rc_best, n_best = [float(x) for x in res_fit.x]
        mu_fit = mu_DK_profile(r, A0_best, rc_best, n_best)
        rho_eff = rho_cdm_base * mu_fit

        return {
            "label": label,
            "A0": A0_best,
            "rc": rc_best,
            "n": n_best,
            "mu": mu_fit,
            "rho_eff": rho_eff,
            "log10_mse": float(objective(res_fit.x)),
            "r_min": r_min,
            "r_max": r_max,
        }

    fit_powell = fit_mu_profile(
        mu_powell_required,
        "B1938+666",
        r_min=0.003,
        r_max=0.12,
    )

    fit_fornax = fit_mu_profile(
        mu_fornax_required,
        "Fornax 6",
        r_min=0.003,
        r_max=0.08,
    )

    fit_gd1 = fit_mu_profile(
        mu_gd1_rep_required,
        "GD-1 representative",
        r_min=0.003,
        r_max=0.12,
    )

    # ------------------------------------------------------------
    # Robustness test: sensitivity to the fitted radial window.
    # ------------------------------------------------------------
    robustness_windows = [
        (0.003, 0.120, "fiducial"),
        (0.004, 0.100, "medium"),
        (0.005, 0.080, "inner_conservative"),
    ]

    robustness_targets = [
        ("B1938+666", mu_powell_required),
        ("Fornax 6", mu_fornax_required),
        ("GD-1 representative", mu_gd1_rep_required),
    ]

    robustness_rows = []

    for target_label, target_mu in robustness_targets:
        for rmin_i, rmax_i, window_label in robustness_windows:
            fit_i = fit_mu_profile(
                target_mu,
                f"{target_label} robustness {window_label}",
                r_min=rmin_i,
                r_max=rmax_i,
            )

            robustness_rows.append({
                "target": target_label,
                "window_label": window_label,
                "r_min_kpc": rmin_i,
                "r_max_kpc": rmax_i,
                "A0_inner_excess": fit_i["A0"],
                "mu0_inner_limit_Geff_over_G0": 1.0 + fit_i["A0"],
                "r_c_kpc": fit_i["rc"],
                "r_c_pc": fit_i["rc"] * 1000.0,
                "transition_index_n": fit_i["n"],
                "log10_space_mse": fit_i["log10_mse"],
                "alpha_DK_global": alpha_dk_best,
                "rd_DK_global_Mpc": rd_dk_best,
            })

    robustness_df = pd.DataFrame(robustness_rows)

    file_robustness = file_table.replace(".csv", "_robustness.csv")

    dkrd2_to_csv(
        robustness_df,
        file_robustness,
        table_kind="FIG06_COMPACT_HALO_FIT_WINDOW_ROBUSTNESS",
        figure_id=6,
        strict=False,
        index=False,
        meta={
            "purpose": "Check sensitivity of compact-halo amplification parameters to the fitted radial window.",
            "mu_DK_form": "mu_DK(r) = 1 + A0 / (1 + (r/r_c)^n)",
            "interpretation": (
                "If mu0 remains within the same order of magnitude across fitting windows, "
                "the compact amplification result is robust. Large order-of-magnitude changes "
                "would indicate parameter degeneracy or strong window dependence."
            ),
            "alpha_DK_global": alpha_dk_best,
            "rd_DK_global_Mpc": rd_dk_best,
        },
    )

    # ------------------------------------------------------------
    # 5. Save evidence table.
    # ------------------------------------------------------------
    table = pd.DataFrame({
        "r_kpc": r,

        "alpha_DK_global": np.full_like(r, alpha_dk_best),
        "rd_DK_global_Mpc": np.full_like(r, rd_dk_best),
        "mu_LCDM_standard_gravity": mu_lcdm_standard,

        "rho_CDM_interpolated_from_data_nb_Msun_kpc3": rho_cdm_base,

        "rho_Powell_B1938_PseudoJaffe_Msun_kpc3": rho_powell,
        "rho_Fornax6_Hernquist_Msun_kpc3": rho_fornax,
        "rho_GD1_low_Msun_kpc3": rho_gd1_low,
        "rho_GD1_high_Msun_kpc3": rho_gd1_high,
        "rho_GD1_representative_geometric_mean_Msun_kpc3": rho_gd1_rep,

        "mu_required_Powell_Geff_over_G0": mu_powell_required,
        "mu_required_Fornax_Geff_over_G0": mu_fornax_required,
        "mu_required_GD1_low_Geff_over_G0": mu_gd1_low_required,
        "mu_required_GD1_high_Geff_over_G0": mu_gd1_high_required,
        "mu_required_GD1_representative_Geff_over_G0": mu_gd1_rep_required,

        "mu_DKRD2_B1938_Geff_over_G0": fit_powell["mu"],
        "mu_DKRD2_Fornax_Geff_over_G0": fit_fornax["mu"],
        "mu_DKRD2_GD1_rep_Geff_over_G0": fit_gd1["mu"],

        "rho_DKRD2_B1938_effective_Msun_kpc3": fit_powell["rho_eff"],
        "rho_DKRD2_Fornax_effective_Msun_kpc3": fit_fornax["rho_eff"],
        "rho_DKRD2_GD1_rep_effective_Msun_kpc3": fit_gd1["rho_eff"],
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="FIG06_COMPACT_HALO_EFFECTIVE_GRAVITY_AMPLIFICATION",
        figure_id=6,
        strict=False,
        index=False,
        meta={
            "notebook_source": str(nb_path),
            "alpha_DK_global": alpha_dk_best,
            "rd_DK_global_Mpc": rd_dk_best,
            "alpha_source": (
                "Global DESI BAO + Cosmic Chronometer calibration; "
                "not fitted to compact halos."
            ),
            "LCDM_baseline": "mu_eff(r)=1, equivalent to G_eff=G0.",
            "mu_DK_form": "1 + A0 / (1 + (r/r_c)^n)",
            "G0_source": "Core_G0 from DK_RD2_Core.py",
            "interpretation": (
                "Three separate DK-RD2 phenomenological compact amplification "
                "profiles are fitted to B1938+666, Fornax 6, and a GD-1 "
                "representative curve. The SIDM interpretation is preserved, "
                "while DK-RD2 is tested as an alternative effective-gravity "
                "description relative to the standard mu_eff=1 baseline."
            ),
        },
    )

    # ------------------------------------------------------------
    # 6. Save stats.
    # ------------------------------------------------------------
    stats_rows = []
    for fit in [fit_powell, fit_fornax, fit_gd1]:
        stats_rows.append({
            "model": f"DK-RD2 {fit['label']} compact amplification",
            "alpha_DK_global": alpha_dk_best,
            "rd_DK_global_Mpc": rd_dk_best,
            "LCDM_standard_gravity_mu": 1.0,
            "A0_inner_excess": fit["A0"],
            "mu_inner_limit_Geff_over_G0": 1.0 + fit["A0"],
            "r_c_kpc": fit["rc"],
            "r_c_pc": fit["rc"] * 1000.0,
            "transition_index_n": fit["n"],
            "log10_space_mse": fit["log10_mse"],
            "fit_radial_window_kpc": f"{fit['r_min']} to {fit['r_max']}",
            "fit_mode": "separate_phenomenological_effective_G_amplification",
        })

    stats = pd.DataFrame(stats_rows)

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=6,
        fit_mode="SEPARATE_PHENOMENOLOGICAL_EFFECTIVE_G_AMPLIFICATION",
        index=False,
    )

    # ------------------------------------------------------------
    # 7. Plot Figure 06.
    # ------------------------------------------------------------
    fig = plt.figure(figsize=(13.5, 8.9), dpi=130)

    gs = gridspec.GridSpec(
        3,
        1,
        height_ratios=[1.55, 1.0, 1.55],
        hspace=0.18,
    )

    ax_den = fig.add_subplot(gs[0])
    ax_mu = fig.add_subplot(gs[1], sharex=ax_den)
    ax_eff = fig.add_subplot(gs[2], sharex=ax_den)

    # ------------------------------------------------------------
    # Panel A — published SIDM / observed compact profiles.
    # ------------------------------------------------------------
    ax_den.loglog(
        cdm_data[:, 0],
        cdm_data[:, 1],
        color=cdm_color,
        linestyle="--",
        linewidth=1.7,
        label="CDM baseline",
    )

    ax_den.loglog(
        sidm30_data[:, 0],
        sidm30_data[:, 1],
        color=sidm30_color,
        linestyle="--",
        linewidth=1.45,
        label="SIDM30",
    )

    ax_den.loglog(
        sidm50_data[:, 0],
        sidm50_data[:, 1],
        color=sidm50_color,
        linestyle="--",
        linewidth=1.45,
        label="SIDM50",
    )

    ax_den.loglog(
        sidm100_data[:, 0],
        sidm100_data[:, 1],
        color=sidm100_color,
        linestyle="--",
        linewidth=1.45,
        label="SIDM100",
    )

    ax_den.loglog(
        r,
        rho_powell,
        color=b1938_color,
        linewidth=1.8,
        label="B1938+666 target",
    )

    ax_den.loglog(
        r,
        rho_fornax,
        color=fornax_color,
        linewidth=1.8,
        label="Fornax 6 target",
    )

    ax_den.fill_between(
        r,
        rho_gd1_low,
        rho_gd1_high,
        color=gd1_band_color,
        alpha=0.32,
        label="GD-1 envelope",
    )

    ax_den.loglog(
        r,
        rho_gd1_rep,
        color=gd1_color,
        linestyle="-.",
        linewidth=1.5,
        alpha=0.70,
        label="GD-1 representative",
    )

    ax_den.set_xlim(0.003, 0.50)
    ax_den.set_ylim(1e7, 2e11)
    ax_den.set_ylabel(
        r"$\rho(r)$ [$M_\odot\,{\rm kpc}^{-3}$]",
        fontsize=6.0,
    )

    ax_den.set_title(
        "Figure 06 — Compact-Halo Prediction and Effective Gravitational Amplification",
        fontsize=11.8,
        pad=12,
    )

    ax_den.text(
        0.5,
        1.02,
        (
            "Published SIDM profiles are preserved; DK-RD2 is tested as a "
            "phenomenological effective-gravity amplification. "
            rf"Global $\alpha_{{\rm DK}}={alpha_dk_best:.4f}$ is shown for consistency."
        ),
        transform=ax_den.transAxes,
        ha="center",
        va="bottom",
        fontsize=7.2,
    )

    ax_den.grid(alpha=0.23, which="both")
    ax_den.legend(
        fontsize=6.4,
        loc="upper right",
        ncols=4,
        framealpha=0.92,
    )

    # ------------------------------------------------------------
    # Panel B — required and fitted effective amplification.
    # ------------------------------------------------------------
    ax_mu.semilogx(
        r,
        mu_powell_required,
        color=b1938_color,
        linewidth=1.25,
        alpha=0.62,
        label=r"Required $\mu$: B1938+666",
    )

    ax_mu.semilogx(
        r,
        mu_fornax_required,
        color=fornax_color,
        linewidth=1.25,
        alpha=0.62,
        label=r"Required $\mu$: Fornax 6",
    )

    ax_mu.fill_between(
        r,
        mu_gd1_low_required,
        mu_gd1_high_required,
        color=gd1_band_color,
        alpha=0.28,
        label=r"Required $\mu$: GD-1 envelope",
    )

    ax_mu.semilogx(
        r,
        mu_gd1_rep_required,
        color=gd1_color,
        linestyle=":",
        linewidth=1.2,
        alpha=0.70,
        label=r"Required $\mu$: GD-1 representative",
    )

    ax_mu.semilogx(
        r,
        fit_powell["mu"],
        color=dk_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.15,
        label="DK-RD2 effective-gravity fit: B1938+666",
    )

    ax_mu.semilogx(
        r,
        fit_fornax["mu"],
        color=dk_light,
        linestyle=DK_RD2_linestyle,
        linewidth=2.15,
        label="DK-RD2 fit: Fornax 6",
    )

    ax_mu.semilogx(
        r,
        fit_gd1["mu"],
        color=secondary_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.15,
        label="DK-RD2 fit: GD-1 rep.",
    )

    ax_mu.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.2,
        alpha=0.90,
        label=r"$\Lambda$CDM / standard gravity: $\mu_{\rm eff}=1$",
    )

    ax_mu.set_yscale("log")
    ax_mu.set_ylim(5e-3, 5e4)
    ax_mu.set_ylabel(
        r"$\mu_{\rm eff}(r)=G_{\rm eff}/G_0$",
        fontsize=6.0,
    )

    ax_mu.grid(alpha=0.23, which="both")
    ax_mu.legend(
        fontsize=5.2,
        loc="upper right",
        ncols=3,
        framealpha=0.92,
    )

    fit_text = (
        "Separate DK-RD2 compact fits:\n"
        rf"B1938: $\mu_0={1+fit_powell['A0']:.1f}$, "
        rf"$r_c={fit_powell['rc']*1000:.2f}$ pc, "
        rf"$n={fit_powell['n']:.2f}$\n"
        rf"Fornax: $\mu_0={1+fit_fornax['A0']:.1f}$, "
        rf"$r_c={fit_fornax['rc']*1000:.2f}$ pc, "
        rf"$n={fit_fornax['n']:.2f}$\n"
        rf"GD-1 rep.: $\mu_0={1+fit_gd1['A0']:.1f}$, "
        rf"$r_c={fit_gd1['rc']*1000:.2f}$ pc, "
        rf"$n={fit_gd1['n']:.2f}$"
    )

    ax_mu.text(
        0.025,
        0.93,
        fit_text,
        transform=ax_mu.transAxes,
        fontsize=5.4,
        va="top",
        ha="left",
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="0.45",
            alpha=0.92,
        ),
    )

    # ------------------------------------------------------------
    # Panel C — effective density proxies.
    # ------------------------------------------------------------
    ax_eff.loglog(
        cdm_data[:, 0],
        cdm_data[:, 1],
        color=cdm_color,
        linestyle="--",
        linewidth=1.5,
        label="CDM baseline",
    )

    ax_eff.loglog(
        r,
        rho_powell,
        color=b1938_color,
        linewidth=1.2,
        alpha=0.55,
        label="B1938+666 target",
    )

    ax_eff.loglog(
        r,
        rho_fornax,
        color=fornax_color,
        linewidth=1.2,
        alpha=0.55,
        label="Fornax 6 target",
    )

    ax_eff.fill_between(
        r,
        rho_gd1_low,
        rho_gd1_high,
        color=gd1_band_color,
        alpha=0.28,
        label="GD-1 envelope",
    )

    ax_eff.loglog(
        r,
        fit_powell["rho_eff"],
        color=dk_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=r"CDM $\times\mu_{\rm DK}$: B1938",
    )

    ax_eff.loglog(
        r,
        fit_fornax["rho_eff"],
        color=dk_light,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=r"CDM $\times\mu_{\rm DK}$: Fornax",
    )

    ax_eff.loglog(
        r,
        fit_gd1["rho_eff"],
        color=secondary_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=r"CDM $\times\mu_{\rm DK}$: GD-1 rep.",
    )

    ax_eff.set_xlim(0.003, 0.50)
    ax_eff.set_ylim(1e7, 2e11)
    ax_eff.set_xlabel(
        "Radius r [kpc]",
        fontsize=5.6,
        labelpad=2,
    )
    ax_eff.set_ylabel(
        r"Effective density proxy [$M_\odot\,{\rm kpc}^{-3}$]",
        fontsize=6.0,
    )

    ax_eff.grid(alpha=0.23, which="both")
    ax_eff.legend(
        fontsize=5.4,
        loc="upper right",
        ncols=3,
        framealpha=0.92,
    )

    ax_eff.text(
        0.025,
        0.065,
        (
            "Interpretation: SIDM / ΛCDM compact-halo readings are preserved.\n"
            "DK-RD2 provides an alternative effective-gravity description using \n"
            "without modifying the published SIDM profiles.\n"
            "system-dependent compact amplification."
        ),
        transform=ax_eff.transAxes,
        fontsize=5.4,
        va="bottom",
        ha="left",
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="0.55",
            alpha=0.90,
        ),
    )

    for ax in [ax_den, ax_mu, ax_eff]:
        ax.tick_params(
            axis="both",
            labelsize=6.0,
        )

    # ------------------------------------------------------------
    # Footer.
    # ------------------------------------------------------------
    fig.text(
        0.5,
        0.012,
        f"Figure: {file_fig} | Source: {Core_autor_text} |\n"
        f"External comparison data: Zenodo Mathematica notebook from http://dx.doi.org/10.5281/zenodo.19116269",
        ha="center",
        fontsize=5.2,
        color=dk_color,
    )

    fig.subplots_adjust(
        top=0.91,
        bottom=0.125,
        left=0.105,
        right=0.985,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats, file_robustness

def generate_figure07(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    z_min: float = 0.0,
    z_max: float = 5.0,
    n_z: int = 900,

    # ============================================================
    ### The selected H0 values are representative of the currently reported
    ### early- and late-Universe determinations and are used only to illustrate
    ### the DK-RD2 thermodynamic projection mechanism.
    # ============================================================
    H0_late_ref: float = 73.0,
    H0_early_ref: float = 67.4,
    f_H0: float = 1.0,
    # ============================================================

):
    """
    # ============================================================
    # Figure 07 — Hubble-Tension Projection Test
    # ============================================================
    #
    # Purpose
    # -------
    # Test whether the thermodynamic projection structure of DK-RD2 can
    # generate distinct effective reconstructions of H0 between early-
    # and late-Universe observational regimes.
    #
    # Physical idea
    # -------------
    # DK-RD2 predicts that the observable thermodynamic correction scales as:
    #
    #     T0 / T(z) = 1 / (1 + z)
    #
    # Therefore, the projected gravitational contribution is suppressed at
    # high redshift and becomes maximal toward the present epoch.
    #
    # This figure does not claim to solve the Hubble tension by fitting H0.
    # Instead, it defines a falsifiable diagnostic:
    #
    #     H0_eff(z) = H0_early + (H0_late - H0_early) * P_DK(z)
    #
    # where P_DK(z) is the normalized DK thermodynamic projection:
    #
    #     P_DK(z) = [alpha_DK/(1+z)] / alpha_DK = 1/(1+z)
    #
    # The global DESI BAO + Cosmic Chronometer calibration alpha_DK is used
    # for consistency with Figures 01–06.
    #
    # Interpretation
    # --------------
    # If future H0 determinations converge to a single thermodynamic-
    # independent value, this DK-RD2 projection interpretation is disfavored.
    # If systematic differences persist and correlate with the effective
    # thermodynamic projection regime, DK-RD2 gains observational support.
    #
    # alpha_DK is not hardcoded.
    # It is obtained from get_alpha_DK(), which reuses the global
    # DESI BAO + Cosmic Chronometer calibration.
    # ============================================================
    """

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure07 requires desi_bao_mean_path and desi_bao_cov_path.")

    # ------------------------------------------------------------
    # Global calibrated alpha_DK access from main()
    alpha_dk_best = alpha_dk_best_global
    rd_dk_best = rd_dk_best_global
    # ------------------------------------------------------------

    file_fig = generate_evidence("image", 7)
    file_table = generate_evidence("table", 7)
    file_stats = file_table.replace(".csv", "_stats.csv")

    # ------------------------------------------------------------
    # Redshift grid and thermodynamic projection.
    # ------------------------------------------------------------
    z = np.linspace(float(z_min), float(z_max), int(n_z))

    T_z = Core_TCMB_K * (1.0 + z)
    thermal_factor = Core_TCMB_K / T_z

    # DK projected excess for beta = 1 reference normalization.
    delta_G_eff_over_G0 = alpha_dk_best * thermal_factor
    G_eff_over_G0 = 1.0 + delta_G_eff_over_G0

    # Normalized projection from 1 at z=0 to suppressed at high z.
    projection_norm = thermal_factor

    # ------------------------------------------------------------
    # Effective H0 projection diagnostic.
    # ------------------------------------------------------------
    H0_late_ref = float(H0_late_ref)
    H0_early_ref = float(H0_early_ref)

    H0_gap = H0_late_ref - H0_early_ref

    H0_eff_DK = H0_early_ref + f_H0 * H0_gap * projection_norm

    H0_late_line = np.full_like(z, H0_late_ref)
    H0_early_line = np.full_like(z, H0_early_ref)

    delta_H0_from_early = H0_eff_DK - H0_early_ref
    delta_H0_from_late = H0_eff_DK - H0_late_ref

    frac_gap_remaining = (H0_eff_DK - H0_early_ref) / np.clip(H0_gap, 1.0e-300, None)

    # ------------------------------------------------------------
    # Representative regimes.
    # ------------------------------------------------------------
    z_samples = np.array([0.0, 0.5, 1.0, 2.0, 5.0], dtype=float)

    def interp_safe(x_grid, y_grid, x_value):
        mask = np.isfinite(x_grid) & np.isfinite(y_grid)
        if np.count_nonzero(mask) < 2:
            return np.nan
        return float(np.interp(float(x_value), x_grid[mask], y_grid[mask]))

    sample_rows = []
    for zi in z_samples:
        sample_rows.append({
            "z": zi,
            "T0_over_Tz": interp_safe(z, thermal_factor, zi),
            "G_eff_over_G0": interp_safe(z, G_eff_over_G0, zi),
            "delta_G_eff_over_G0": interp_safe(z, delta_G_eff_over_G0, zi),
            "H0_eff_DK_km_s_Mpc": interp_safe(z, H0_eff_DK, zi),
            "delta_H0_from_early": interp_safe(z, delta_H0_from_early, zi),
            "delta_H0_from_late": interp_safe(z, delta_H0_from_late, zi),
        })

    sample_df = pd.DataFrame(sample_rows)

    # ------------------------------------------------------------
    # Save evidence table.
    # ------------------------------------------------------------
    table = pd.DataFrame({
        "z": z,
        "alpha_DK": np.full_like(z, alpha_dk_best),
        "rd_DK_Mpc": np.full_like(z, rd_dk_best),

        "T_CMB_z_K": T_z,
        "T0_over_Tz": thermal_factor,

        "delta_G_eff_over_G0_beta1": delta_G_eff_over_G0,
        "G_eff_over_G0_beta1": G_eff_over_G0,
        "projection_normalized": projection_norm,

        "H0_early_reference": H0_early_line,
        "H0_late_reference": H0_late_line,
        "H0_eff_DK_projection": H0_eff_DK,

        "delta_H0_from_early": delta_H0_from_early,
        "delta_H0_from_late": delta_H0_from_late,
        "fraction_of_H0_gap_projected": frac_gap_remaining,
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="FIG07_HUBBLE_TENSION_PROJECTION_TEST",
        figure_id=7,
        strict=False,
        index=False,
        meta={
            "alpha_DK": alpha_dk_best,
            "rd_DK_Mpc": rd_dk_best,
            "H0_early_reference": H0_early_ref,
            "H0_late_reference": H0_late_ref,
            "H0_gap": H0_gap,
            "z_range": f"{z_min} to {z_max}",
            "alpha_source": "Global DESI BAO + Cosmic Chronometer calibration",
            "thermal_projection": "T(z)=T0(1+z), T0/T(z)=1/(1+z)",
            "H0_projection_model": (
                "H0_eff(z)=H0_early+(H0_late-H0_early)/(1+z). "
                "This is a diagnostic projection, not an H0 fit."
            ),
            "H0_projection_fraction_f_H0": f_H0,
            "H0_projection_model": (
                "H0_eff(z)=H0_early+f_H0*(H0_late-H0_early)/(1+z). "
                "This is a diagnostic projection, not an H0 fit."
            ),
            "interpretation": (
                "Tests whether early- and late-Universe H0 determinations "
                "can correspond to different thermodynamic projection regimes."
            ),
        },
    )

    # ------------------------------------------------------------
    # Save stats.
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "DK-RD2 Hubble-tension projection diagnostic",
            "alpha_DK": float(alpha_dk_best),
            "rd_mpc": float(rd_dk_best),
            "H0_early_reference": H0_early_ref,
            "H0_late_reference": H0_late_ref,
            "f_H0_projection_fraction": float(f_H0),
            "H0_gap": H0_gap,
            "H0_eff_z0": interp_safe(z, H0_eff_DK, 0.0),
            "H0_eff_z0p5": interp_safe(z, H0_eff_DK, 0.5),
            "H0_eff_z1": interp_safe(z, H0_eff_DK, 1.0),
            "H0_eff_z2": interp_safe(z, H0_eff_DK, 2.0),
            "H0_eff_z5": interp_safe(z, H0_eff_DK, 5.0),
            "T0_over_T_z0": interp_safe(z, thermal_factor, 0.0),
            "T0_over_T_z1": interp_safe(z, thermal_factor, 1.0),
            "T0_over_T_z5": interp_safe(z, thermal_factor, 5.0),
            "delta_G_eff_over_G0_z0": interp_safe(z, delta_G_eff_over_G0, 0.0),
            "delta_G_eff_over_G0_z1": interp_safe(z, delta_G_eff_over_G0, 1.0),
            "delta_G_eff_over_G0_z5": interp_safe(z, delta_G_eff_over_G0, 5.0),
            "fit_mode": "HUBBLE_TENSION_THERMODYNAMIC_PROJECTION_DIAGNOSTIC",
        }
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=7,
        fit_mode="HUBBLE_TENSION_PROJECTION_TEST",
        index=False,
    )

    # ------------------------------------------------------------
    # Plot configuration.
    # ------------------------------------------------------------
    fig, (ax_h0, ax_proj, ax_g) = plt.subplots(
        3,
        1,
        figsize=(11.8, 10.2),
        dpi=130,
        sharex=True,
        gridspec_kw={
            "height_ratios": [2.0, 1.15, 1.15],
            "hspace": 0.15,
        },
    )

    title_fs = 14
    label_fs = 9.5
    tick_fs = 8.2
    legend_fs = 7.4
    box_fs = 6.8
    footer_fs = 5.4

    # ------------------------------------------------------------
    # Panel A — effective H0 projection.
    # ------------------------------------------------------------
    ax_h0.plot(
        z,
        H0_late_line,
        color="0.25",
        linestyle=":",
        linewidth=1.8,
        label=rf"Late-Universe reference $H_0={H0_late_ref:.1f}$",
    )

    ax_h0.plot(
        z,
        H0_early_line,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=2.0,
        label=rf"Early-Universe reference $H_0={H0_early_ref:.1f}$",
    )

    ax_h0.plot(
        z,
        H0_eff_DK,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.4,
        label=r"DK-RD2 projected $H_0^{eff}(z)$",
    )

    ax_h0.set_ylabel(
        r"$H_0^{eff}$ [km s$^{-1}$ Mpc$^{-1}$]",
        fontsize=label_fs,
    )

    ax_h0.set_title(
        rf"Figure 07 — Hubble-Tension Projection Test"
        "\n"
        rf"$\alpha_{{\rm DK}}={alpha_dk_best:.4f}$, "
        rf"$r_d={rd_dk_best:.2f}\,\mathrm{{Mpc}}$ "
        rf"(global DESI BAO + CC calibration)",
        fontsize=title_fs,
        pad=8,
    )

    ax_h0.tick_params(axis="both", labelsize=tick_fs)
    ax_h0.grid(alpha=0.28)
    ax_h0.legend(fontsize=legend_fs, loc="best", framealpha=0.88)

    delta_H0_z1 = interp_safe(z, delta_H0_from_early, 1.0)

    result_text = (
        r"Projection diagnostic:" "\n"
        r"Maximal projection case" "\n"
        rf"$H_0^{{early}}={H0_early_ref:.1f}$" "\n"
        rf"$H_0^{{late}}={H0_late_ref:.1f}$" "\n"
        rf"$\Delta H_0={H0_gap:.1f}$" "\n"        
        rf"$H_0^{{eff}}(z)=H_0^{{early}}+f_{{H0}}\Delta H_0/(1+z)$" "\n"
        rf"$f_{{H0}}={f_H0:.2f}$" "\n"
        rf"$\Delta H_0^{{DK}}(z=1)={delta_H0_z1:.1f}$"
    )

    ax_h0.text(
        0.035,
        0.55,
        result_text,
        transform=ax_h0.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs-2,
        linespacing=1.1,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.8,
            alpha=0.88,
        ),
    )

    # ------------------------------------------------------------
    # Panel B — thermal projection.
    # ------------------------------------------------------------
    ax_proj.plot(
        z,
        projection_norm,
        color="black",
        linestyle="--",
        linewidth=1.9,
        label=r"$T_0/T(z)=1/(1+z)$",
    )

    ax_proj.plot(
        z,
        frac_gap_remaining,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.1,
        label=r"Projected fraction of the $H_0$ gap",
    )

    ax_proj.set_ylabel("Projection factor", fontsize=label_fs)
    ax_proj.tick_params(axis="both", labelsize=tick_fs)
    ax_proj.grid(alpha=0.28)
    ax_proj.legend(fontsize=legend_fs, loc="best", framealpha=0.88)

    note_text = (
        r"$T(z)=T_0(1+z)$" "\n"
        r"High-$z$ observations sample suppressed DK projection." "\n"
        r"Low-$z$ observations sample the maximal projected regime."
    )

    ax_proj.text(
        0.2,
        0.90,
        note_text,
        transform=ax_proj.transAxes,
        ha="left",
        va="top",
        fontsize=6.3,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.86,
        ),
    )

    # ------------------------------------------------------------
    # Panel C — effective coupling projection.
    # ------------------------------------------------------------
    ax_g.plot(
        z,
        G_eff_over_G0,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=2.1,
        label=r"$G_{eff}(z)/G_0=1+\alpha_{\rm DK}/(1+z)$",
    )

    ax_g.plot(
        z,
        delta_G_eff_over_G0,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.0,
        alpha=0.85,
        label=r"$\Delta G_{eff}/G_0=\alpha_{\rm DK}/(1+z)$",
    )

    ax_g.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.0,
        alpha=0.75,
        label=r"Standard gravity $G_{eff}=G_0$",
    )

    ax_g.set_xlabel("Redshift z", fontsize=label_fs)
    ax_g.set_ylabel(r"Coupling projection", fontsize=label_fs)
    ax_g.tick_params(axis="both", labelsize=tick_fs)
    ax_g.grid(alpha=0.28)
    ax_g.legend(fontsize=legend_fs, loc="best", framealpha=0.88)

    source_text = (
        r"$\Delta G_{eff}/G_0=\alpha_{\rm DK}\,T_0/T(z)$" "\n"
        rf"$\alpha_{{\rm DK}}={alpha_dk_best:.4f}$ is fixed by BAO + CC."
    )

    ax_g.text(
        0.035,
        0.7,
        source_text,
        transform=ax_g.transAxes,
        ha="left",
        va="top",
        fontsize=box_fs,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.86,
        ),
    )

    # ------------------------------------------------------------
    # Footer.
    # ------------------------------------------------------------
    fig.text(
        0.5,
        0.018,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=footer_fs,
        color=DK_RD2_color,
    )

    fig.subplots_adjust(
        top=0.900,
        bottom=0.105,
        left=0.085,
        right=0.965,
        hspace=0.17,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats

def generate_figure08(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    a_min: float = 1.0e-10,
    a_max: float = 1.0e0,
    n_a: int = 1200,
):
    """
    # ============================================================
    # Figure 08 — Competing Physical Mechanisms:
    # Dark-Sector Explanations versus DK-RD2 Thermodynamic Projection
    # ============================================================
    #
    # Purpose
    # -------
    # Provide a conceptual model-discrimination figure comparing:
    #
    #   1. ΛCDM plus additional dark-sector mechanisms often invoked
    #      to address recent DESI / Hubble-tension-related anomalies.
    #
    #   2. DK-RD2, where the same class of observational signatures is
    #      generated by the thermodynamic evolution of the effective
    #      gravitational coupling.
    #
    # Scientific role
    # ---------------
    # This figure does NOT reproduce any published Hot NEDE / DAO figure.
    # It creates an original schematic inspired by the same physical logic:
    # extra dark-sector energy components, phase-transition epochs, ΔNeff,
    # and dark-radiation interactions.
    #
    # The middle panel shows the DK-RD2 alternative:
    #
    #       G_eff(z)/G0 = 1 + alpha_DK / (1 + z)
    #
    # with alpha_DK fixed by the global DESI BAO + Cosmic Chronometer
    # calibration used throughout Figures 01–07.
    #
    # The lower panel summarizes the falsifiable observational roadmap:
    # which measurements can discriminate between DK-RD2 and dark-sector
    # explanations.
    #
    # Interpretation
    # --------------
    # The comparison emphasizes physical mechanism and falsifiability,
    # not merely goodness of fit.
    # ============================================================
    """

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure08 requires desi_bao_mean_path and desi_bao_cov_path.")

    # ------------------------------------------------------------
    # Global calibrated alpha_DK access from main()
    alpha_dk_best = alpha_dk_best_global
    rd_dk_best = rd_dk_best_global
    # ------------------------------------------------------------

    file_fig = generate_evidence("image", 8)
    file_table = generate_evidence("table", 8)
    file_stats = file_table.replace(".csv", "_stats.csv")
    file_signatures = file_table.replace(".csv", "_signatures.csv")

    # ------------------------------------------------------------
    # Scale-factor grid and redshift conversion.
    # ------------------------------------------------------------
    a = np.logspace(np.log10(float(a_min)), np.log10(float(a_max)), int(n_a))
    z = (1.0 / a) - 1.0
    z_safe = np.clip(z, 0.0, 1.0e12)

    # ------------------------------------------------------------
    # Panel A: original schematic dark-sector ingredients.
    #
    # These curves are qualitative templates, not a fit to external data.
    # They illustrate extra components typically invoked in early-dark-energy,
    # dark-radiation, Hot NEDE, DAO, or DRMD scenarios.
    # ------------------------------------------------------------
    loga = np.log10(a)

    def smooth_step(x, x0, width):
        return 1.0 / (1.0 + np.exp(-(x - x0) / width))

    a_bbn_min = 1.0e-10
    a_bbn_max = 2.0e-9
    a_phase = 1.0e-6
    a_ir = 1.8e-4

    rho_dr = 0.045 + 0.50 * smooth_step(loga, np.log10(a_phase), 0.055)
    rho_nede = 0.012 + 0.52 * smooth_step(loga, np.log10(a_phase), 0.030)
    rho_gauge = 0.035 + 0.45 * smooth_step(loga, np.log10(a_phase), 0.040)
    rho_higgs = 0.012 + 0.065 * (1.0 - smooth_step(loga, np.log10(a_ir), 0.060))
    delta_neff = 0.045 + 0.45 * smooth_step(loga, np.log10(a_phase), 0.040)

    rho_dr += 0.10 * smooth_step(loga, np.log10(a_ir), 0.10)
    rho_gauge += 0.07 * smooth_step(loga, np.log10(a_ir), 0.10)
    delta_neff += 0.08 * smooth_step(loga, np.log10(a_ir), 0.10)

    # ------------------------------------------------------------
    # Panel B: DK-RD2 thermodynamic projection.
    # ------------------------------------------------------------
    thermal_projection = a.copy()  # T0/T(z)=1/(1+z)=a
    delta_G_over_G0 = alpha_dk_best * thermal_projection
    G_eff_over_G0 = 1.0 + delta_G_over_G0

    # DK-RD2 reference curve shown in Panel A
    # normalized to the dark-sector schematic scale

    Geff_over_G0_norm = (
                                G_eff_over_G0 - np.min(G_eff_over_G0)
                        ) / (
                                np.max(G_eff_over_G0) - np.min(G_eff_over_G0)
                        )

    Geff_over_G0_norm = 0.04 + 0.50 * Geff_over_G0_norm

    # ------------------------------------------------------------
    # Panel C: falsification roadmap.
    #
    # This is not a likelihood matrix. It summarizes which observables
    # distinguish DK-RD2 from dark-sector explanations.
    # ------------------------------------------------------------
    roadmap_rows = [
        {
            "observable": "H(z), distances",
            "DK_RD2_prediction": "Projected thermodynamic expansion response",
            "Dark_sector_prediction": "Expansion modified by Λ, w(z), or early-sector energy injection",
            "primary_DKRD2_figure": "Figures 01, 03, 04",
            "falsification_test": "Redshift-dependent expansion reconstruction",
        },
        {
            "observable": "w_eff(z)",
            "DK_RD2_prediction": "Emergent diagnostic w(z), no fundamental dark-energy fluid",
            "Dark_sector_prediction": "Physical evolving dark-energy equation of state",
            "primary_DKRD2_figure": "Figure 04",
            "falsification_test": "FRW reconstruction versus physical dark-energy dynamics",
        },
        {
            "observable": "fσ8(z)",
            "DK_RD2_prediction": "Small growth deviation from the same projected background",
            "Dark_sector_prediction": "Growth altered by dark-sector dynamics or expansion changes",
            "primary_DKRD2_figure": "Figure 05",
            "falsification_test": "Independent structure-growth measurements",
        },
        {
            "observable": "Lensing geometry",
            "DK_RD2_prediction": "Distance/lensing deviations correlated with effective coupling",
            "Dark_sector_prediction": "Standard geometry plus dark-sector parameter shifts",
            "primary_DKRD2_figure": "Figures 04, future lensing",
            "falsification_test": "Weak/strong lensing distance reconstruction",
        },
        {
            "observable": "Compact halos",
            "DK_RD2_prediction": "Local μ_eff(r)=G_eff/G0 amplification",
            "Dark_sector_prediction": "Dense dark-matter or SIDM halo structure",
            "primary_DKRD2_figure": "Figure 06",
            "falsification_test": "Compact-halo gravitational response",
        },
        {
            "observable": "H0 reconstruction",
            "DK_RD2_prediction": "Thermodynamic projection between early and late regimes",
            "Dark_sector_prediction": "Early-energy or dark-radiation solution",
            "primary_DKRD2_figure": "Figure 07",
            "falsification_test": "Correlation of H0 inference with thermodynamic regime",
        },
    ]

    roadmap_df = pd.DataFrame(roadmap_rows)

    # ------------------------------------------------------------
    # Save evidence table.
    # ------------------------------------------------------------
    table = pd.DataFrame({
        "a": a,
        "z": z_safe,
        "alpha_DK_global": np.full_like(a, alpha_dk_best),
        "rd_DK_global_Mpc": np.full_like(a, rd_dk_best),

        "dark_sector_template_rho_dr_over_rho_nu": rho_dr,
        "dark_sector_template_rho_nede_over_rho_nu": rho_nede,
        "dark_sector_template_rho_gauge_over_rho_nu": rho_gauge,
        "dark_sector_template_rho_higgs_over_rho_nu": rho_higgs,
        "dark_sector_template_delta_neff": delta_neff,

        "DKRD2_T0_over_Tz": thermal_projection,
        "DKRD2_delta_Geff_over_G0": delta_G_over_G0,
        "DKRD2_Geff_over_G0": G_eff_over_G0,
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="FIG08_DARK_SECTOR_VS_DKRD2_MECHANISM",
        figure_id=8,
        strict=False,
        index=False,
        meta={
            "alpha_DK_global": alpha_dk_best,
            "rd_DK_global_Mpc": rd_dk_best_global,
            "alpha_source": "Global DESI BAO + Cosmic Chronometer calibration",
            "dark_sector_panel": (
                "Qualitative original schematic inspired by Hot NEDE / DAO / DRMD-type "
                "mechanisms; not a reproduction of any external figure."
            ),
            "DKRD2_panel": (
                "Thermodynamic projection Geff/G0=1+alpha_DK/(1+z), "
                "using the same alpha_DK as Figures 01-07."
            ),
            "roadmap_panel": (
                "Falsification roadmap summarizing observational domains able to "
                "discriminate DK-RD2 from dark-sector explanations."
            ),
            "interpretation": (
                "Figure compares physical mechanisms and falsifiable observational "
                "targets rather than statistical fit quality."
            ),
        },
    )

    dkrd2_to_csv(
        roadmap_df,
        file_signatures,
        table_kind="FIG08_FALSIFICATION_ROADMAP",
        figure_id=8,
        strict=False,
        index=False,
        meta={
            "alpha_DK_global": alpha_dk_best,
            "rd_DK_global_Mpc": rd_dk_best_global,
            "purpose": (
                "Summarize the observational tests that distinguish DK-RD2 "
                "from dark-sector mechanisms."
            ),
        },
    )

    # ------------------------------------------------------------
    # Save diagnostic stats.
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "Dark-sector schematic",
            "primary_mechanism": "Additional early dark-sector components",
            "requires_extra_dark_energy_sector": True,
            "requires_dark_radiation_or_delta_neff": True,
            "requires_phase_transition_template": True,
            "requires_dark_sector_interactions": True,
            "uses_global_alpha_DK": False,
            "alpha_DK_global": np.nan,
            "rd_DK_global_Mpc": np.nan,
            "fit_mode": "QUALITATIVE_MECHANISM_COMPARISON",
        },
        {
            "model": "DK-RD2 effective-gravity template",
            "primary_mechanism": "Thermodynamic evolution of effective gravitational coupling",
            "requires_extra_dark_energy_sector": False,
            "requires_dark_radiation_or_delta_neff": False,
            "requires_phase_transition_template": False,
            "requires_dark_sector_interactions": False,
            "uses_global_alpha_DK": True,
            "alpha_DK_global": alpha_dk_best,
            "rd_DK_global_Mpc": rd_dk_best_global,
            "Geff_over_G0_z0": float(1.0 + alpha_dk_best),
            "Geff_over_G0_z1": float(1.0 + alpha_dk_best / 2.0),
            "Geff_over_G0_z3": float(1.0 + alpha_dk_best / 4.0),
            "fit_mode": "QUALITATIVE_MECHANISM_COMPARISON",
        },
        {
            "model": "DK-RD2 falsification roadmap",
            "primary_mechanism": "Multi-observable discrimination",
            "number_of_observable_tests": int(len(roadmap_df)),
            "observable_tests": ", ".join(roadmap_df["observable"].tolist()),
            "alpha_DK_global": alpha_dk_best,
            "rd_DK_global_Mpc": rd_dk_best,
            "fit_mode": "FALSIFICATION_ROADMAP",
        },
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=8,
        fit_mode="MODEL_DISCRIMINATION_MECHANISM_COMPARISON",
        index=False,
    )

    # ------------------------------------------------------------
    # Plot Figure 08.
    # ------------------------------------------------------------
    fig, (ax_dark, ax_dk, ax_road) = plt.subplots(
        3,
        1,
        figsize=(12.8, 10.2),
        dpi=130,
        gridspec_kw={
            "height_ratios": [1.55, 1.25, 1.35],
            "hspace": 0.25,
        },
    )

    title_fs = 13.5
    label_fs = 9.0
    tick_fs = 7.6
    legend_fs = 6.6
    note_fs = 6.4
    table_fs = 6.2
    footer_fs = 5.4

    # ------------------------------------------------------------
    # Panel A — dark-sector mechanism schematic.
    # ------------------------------------------------------------
    ax_dark.set_xscale("log")
    ax_dark.set_yscale("log")

    ax_dark.plot(
        a,
        rho_dr,
        color= LCDM_color,
        linewidth=1.8,
        label=r"dark radiation / $\Delta N_{\rm eff}$ template",
    )

    ax_dark.plot(
        a,
        rho_nede,
        color="#cc8500",
        linewidth=1.8,
        label=r"early dark-energy injection template",
    )

    ax_dark.plot(
        a,
        rho_gauge,
        color="#996300",
        linewidth=1.8,
        label=r"dark gauge-sector template",
    )

    ax_dark.plot(
        a,
        rho_higgs,
        color="#664200",
        linewidth=1.8,
        label=r"phase-transition field template",
    )

    ax_dark.plot(
        a,
        delta_neff,
        color="#332100",
        linewidth=1.6,
        linestyle="--",
        label=r"$\Delta N_{\rm eff}$ effective template",
    )

    # DK-RD2 reference mechanism
    ax_dark.fill_between(
        a,
        0.95 * Geff_over_G0_norm,
        1.05 * Geff_over_G0_norm,
        color="#6666ff",
        alpha=0.12,
    )

    ax_dark.plot(
        a,
        Geff_over_G0_norm,
        color=DK_RD2_color,
        linewidth=2.8,
        label="DK-RD2 projected H0eff(z)",
    )

    ax_dark.axvspan(
        a_bbn_min,
        a_bbn_max,
        color="#ffb833",
        alpha=0.35,
        label="BBN constrained region",
    )

    ax_dark.axvline(
        a_phase,
        color="0.35",
        linestyle=":",
        linewidth=1.1,
    )

    ax_dark.axvline(
        a_ir,
        color="0.35",
        linestyle=":",
        linewidth=1.1,
    )

    ax_dark.text(
        1.6e-10,
        0.18,
        "BBN\nconstraint",
        fontsize=8.0,
        ha="left",
        va="center",
    )

    ax_dark.text(
        1.15e-6,
        0.018,
        "early-sector\ntransition",
        fontsize=7.2,
        ha="left",
        va="bottom",
    )

    ax_dark.text(
        2.1e-4,
        0.018,
        "IR / late\nsector",
        fontsize=7.2,
        ha="left",
        va="bottom",
    )

    ax_dark.text(
        0.012,
        0.92,
        (
            "Panel A — Dark-sector explanations\n"
            "DESI / Hubble-related deviations are modeled using extra sectors:\n"
            r"$\rho_{\rm EDE}$, dark radiation, phase transitions, $\Delta N_{\rm eff}$, or DAO/DRMD interactions."
        ),
        transform=ax_dark.transAxes,
        ha="left",
        va="top",
        fontsize=note_fs,
        linespacing=1.15,
        bbox=dict(
            boxstyle="round,pad=0.30",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.8,
            alpha=0.90,
        ),
    )

    ax_dark.set_xlim(a_min, 1.0e-3)
    ax_dark.set_ylim(8.0e-3, 1.6)
    ax_dark.set_ylabel(r"relative dark-sector scale", fontsize=label_fs)
    ax_dark.set_title(
        "Figure 08 — Competing Physical Mechanisms and Falsification Roadmap",
        fontsize=title_fs,
        pad=10,
    )

    ax_dark.grid(alpha=0.28, which="both")
    ax_dark.tick_params(axis="both", labelsize=tick_fs)
    ax_dark.legend(
        fontsize=legend_fs-3,
        loc="upper center",
        framealpha=0.88,
        ncols=2,
    )

    # ------------------------------------------------------------
    # Panel B — DK-RD2 thermodynamic projection.
    # ------------------------------------------------------------
    ax_dk.set_xscale("log")

    ax_dk.plot(
        a,
        thermal_projection,
        color="black",
        linestyle="--",
        linewidth=1.8,
        label=r"$T_0/T(z)=a=1/(1+z)$",
    )

    ax_dk.plot(
        a,
        delta_G_over_G0,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.1,
        label=r"$\Delta G_{\rm eff}/G_0=\alpha_{\rm DK}a$",
    )

    ax_dk.plot(
        a,
        G_eff_over_G0,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=2.2,
        label=r"$G_{\rm eff}/G_0=1+\alpha_{\rm DK}a$",
    )

    ax_dk.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.0,
        alpha=0.7,
        label=r"standard gravity $G_{\rm eff}=G_0$",
    )

    ax_dk.text(
        0.012,
        0.92,
        (
            "Panel B — DK-RD2 effective-gravity projection\n"
            "No additional dark-energy sector, dark radiation, phase transition, or primordial-spectrum modification is introduced.\n"
            rf"One global calibration fixes the projection: $\alpha_{{\rm DK}}={alpha_dk_best:.4f}$."
        ),
        transform=ax_dk.transAxes,
        ha="left",
        va="center",
        fontsize=note_fs,
        linespacing=1.15,
        bbox=dict(
            boxstyle="round,pad=0.30",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.8,
            alpha=0.90,
        ),
    )

    ax_dk.set_xlim(a_min, 1.0)
    ax_dk.set_ylim(-0.04, 1.38)
    ax_dk.set_ylabel("DK-RD2 projection", fontsize=label_fs)
    ax_dk.grid(alpha=0.28, which="both")
    ax_dk.tick_params(axis="both", labelsize=tick_fs)
    ax_dk.legend(
        fontsize=legend_fs,
        loc="center",
        framealpha=0.88,
        ncols=2,
    )

    # ------------------------------------------------------------
    # Panel C — falsification roadmap table.
    # ------------------------------------------------------------
    ax_road.axis("off")

    ax_road.text(
        0.012,
        0.98,
        (
            "Panel C — Observational falsification roadmap\n"
            "DK-RD2 and dark-sector models may fit similar datasets, but they predict different correlated signatures."
        ),
        transform=ax_road.transAxes,
        ha="left",
        va="top",
        fontsize=note_fs,
        linespacing=1.15,
        bbox=dict(
            boxstyle="round,pad=0.30",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.8,
            alpha=0.90,
        ),
    )

    table_columns = [
        "Observable",
        "DK-RD2 signature",
        "Typical alternative explanation",
        "Test",
    ]

    table_cells = [
        [
            "H(z)",
            "Projected expansion",
            "Λ / w(z) / early energy",
            "DESI, Euclid",
        ],
        [
            r"$w_{\rm eff}(z)$",
            "Emergent diagnostic",
            "Physical DE fluid",
            "FRW reconstruction",
        ],
        [
            r"$f\sigma_8(z)$",
            "Background-linked growth",
            "Growth from sector dynamics",
            "LSS surveys",
        ],
        [
            "Lensing",
            "Coupling-correlated geometry",
            "Parameter-shift geometry",
            "Weak/strong lensing",
        ],
        [
            "Compact halos",
            r"local $\mu_{\rm eff}(r)$",
            "DM/SIDM density structure",
            "sub-galactic systems",
        ],
        [
            r"$H_0$",
            "thermal projection",
            "early-energy solution",
            "early vs late probes",
        ],
    ]

    roadmap_table = ax_road.table(
        cellText=table_cells,
        colLabels=table_columns,
        cellLoc="center",
        loc="center",
        bbox=[0.02, 0.02, 0.96, 0.78],
    )

    roadmap_table.auto_set_font_size(False)
    roadmap_table.set_fontsize(table_fs)

    for (row, col), cell in roadmap_table.get_celld().items():
        cell.set_edgecolor("0.60")
        cell.set_linewidth(0.6)

        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor("#3b4cc0")
        else:
            if col == 1:
                cell.set_facecolor("#eaf2ff")
            elif col == 2:
                cell.set_facecolor("#fff2e0")
            else:
                cell.set_facecolor("white")

    # ------------------------------------------------------------
    # Footer.
    # ------------------------------------------------------------
    fig.text(
        0.5,
        0.018,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=footer_fs,
        color=DK_RD2_color,
    )

    fig.subplots_adjust(
        top=0.925,
        bottom=0.090,
        left=0.085,
        right=0.965,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats

def set_plot_text_defaults():
    """
    Ensure consistent, English-only plot styling across the project.
    """
    import matplotlib as mpl
    mpl.rcParams["axes.titlesize"] = 12
    mpl.rcParams["axes.labelsize"] = 11
    mpl.rcParams["legend.fontsize"] = 9
    mpl.rcParams["figure.titlesize"] = 13
    mpl.rcParams["font.family"] = "DejaVu Sans"  # or another available family
    # All labels/titles you set elsewhere should be English strings.


if __name__ == "__main__":
    import os

    Core_out_dir_path = "evidence"
    os.makedirs(Core_out_dir_path, exist_ok=True)

    set_plot_text_defaults()

    # ==============================================================
    # DK-RD2 ThermoGravity Framework — Figure selector
    #
    # Examples:
    #   exec_figs = ""        # run all paper figures
    #   exec_figs = "all"     # run all paper figures
    #   exec_figs = "1,2"     # run Figure 01 and Figure 02 only
    #   exec_figs = "1-4"     # run Figure 01..04
    #
    # Figure IDs for this paper:
    #   1 — Figure 01: Stellar thermodynamic activity vs gravitational confinement
    #   2 — Figure 02: Effective equation of state w(z): DK-RD2 vs ΛCDM
    #   3 — Figure 03: Type Ia SN distance modulus μ(z): DK-RD2 vs ΛCDM
    #   4 — Figure 04: Degravification / BOAT-like transient diagnostic
    #   5 — Figure 05: Falsifiability summary table: DK-RD2 vs ΛCDM
    #
    # Optional / inherited diagnostics:
    #   8 — Optional: Linear growth observable fσ8(z)
    #   9 — Optional: CLASS/CMB geometric consistency, if retained as appendix
    #
    # Removed from this paper main sequence:
    #   old Fig.04 H(z)        -> already covered in previous DK-RD2 paper
    #   old Fig.05 Lensing SIS -> already covered in previous DK-RD2 paper
    #   old SPHEREx proxy      -> not required here
    # ==============================================================

    banner = r"""
░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░      ░▒▓███████▓▒░░▒▓███████▓▒░░▒▓███████▓▒░                                                          
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░                                                         
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░                                                         
░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░  ▒█▒  ░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░                                                          
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░                                                                
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░                                                                
░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓████████▓▒░
"""

    print("Falsifiable Predictions of the DK-RD2")
    print("DK-RD2 ThermoGravity Framework")
    print(Core_autor_text)

    exec_figs = "all"

    print(
        "====================== Figure IDs for this paper ======================\n"
        "0  - Run all paper figures\n"
        "\n"
        "1  — Figure 01: Statistical calibration of alpha_DK and r_d\n"
        "2  — Figure 02: Likelihood constraint on alpha_DK\n"
        "3  — Figure 03: Effective expansion history and dynamical dark-energy reconstruction\n"
        "4  — Figure 04: Geometric gravitational-lensing signatures\n"
        "5  — Figure 05: Growth of cosmic structure through f_sigma8(z)\n"
        "6  — Figure 06: Compact-halo prediction and effective gravitational amplification\n"
        "7  — Figure 07: Hubble-tension projection test\n"
        "8  — Figure 08: Model discrimination versus DESI alternatives\n"
        "\n"
        "X  — Exit or Ctrl+C"
    )

    print("selector examples:")
    print('  exec_figs = "" or "0" # run all paper figures')
    print('  exec_figs = "all"     # run all paper figures')
    print('  exec_figs = "1,2"     # run Figure 01 and Figure 02 only')
    print('  exec_figs = "1-5"     # run all main paper figures\n')

    print(f"Input selector value: {exec_figs!r}")
    user_spec = input(
        "Type a new selection (e.g. 1,2 or 1-5; 0 = all; ENTER = keep current): "
    ).strip()

    if user_spec == "":
        selection_spec = exec_figs
    elif user_spec == "0":
        selection_spec = "all"
    else:
        selection_spec = user_spec

    ALL_FIGS = {1, 2, 3, 4, 5, 6, 7, 8}
    OPTIONAL_FIGS = set()
    AVAILABLE_FIGS = ALL_FIGS.union(OPTIONAL_FIGS)

    def _parse_exec_list(spec: str | None, available: set[int]) -> set[int]:
        """
        Parse execution specification into a set of figure numbers.

        Supported:
            None / "" / "all"  -> all main paper figures
            "1,2,5"            -> explicit list
            "1-5,8"            -> ranges plus optional diagnostics
        """
        if spec is None:
            return set(ALL_FIGS)

        s = str(spec).strip()

        if s == "" or s.lower() == "all":
            return set(ALL_FIGS)

        out: set[int] = set()
        parts = [p.strip() for p in s.split(",") if p.strip()]

        for p in parts:
            if "-" in p:
                a, b = [x.strip() for x in p.split("-", 1)]
                if a.isdigit() and b.isdigit():
                    lo, hi = int(a), int(b)
                    if lo > hi:
                        lo, hi = hi, lo
                    out.update(range(lo, hi + 1))
            else:
                if p.isdigit():
                    out.add(int(p))

        return out.intersection(available)

    RUN = _parse_exec_list(selection_spec, AVAILABLE_FIGS)

    print(banner)
    print(f"▶ Exec selector = {selection_spec!r} → running: {sorted(RUN)}")
    if selection_spec == "all":
        print ("Generating Déjà vu figures...patience, the matrix... had been changed.")

    path = "data/DESI/bao_data/desi_bao_dr2/"
    desi_bao_mean_path = path + "desi_gaussian_bao_ALL_GCcomb_mean.txt"
    desi_bao_cov_path = path + "desi_gaussian_bao_ALL_GCcomb_cov.txt"

    # ------------------------------------------------------------
    # Official global DK-RD2 calibration.
    # This must be computed once and reused by all figures.
    # ------------------------------------------------------------
    alpha_dk_best, rd_dk_best = calibrate_alpha_DK_global(
        desi_bao_mean_path=desi_bao_mean_path,
        desi_bao_cov_path=desi_bao_cov_path,
        force=False,
        verbose=True,
    )

    if 1 in RUN:
        # ============================================================
        # Figure 01 — DESI Expansion and BAO Calibration of α_DK
        # ============================================================
        print("Figure01  DESI Expansion and BAO Calibration of α_DK")
        print("Generating figure... hang tight, this may take a while.")
        fig1_png, fig1_csv, fig1_stats = generate_figure01(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 01 saved: {fig1_png}")
        print(f"✔ Table saved:     {fig1_csv}")

    # ============================================================
    # Figure 02 — Statistical Constraint on the DK Projection Factor
    # ============================================================
    if 2 in RUN:
        print("Figure 02: Statistical Constraint on the DK Projection Factor")
        print("Generating figure... hang tight, this may take a while.")
        print("... or more than a while...")
        fig2_img, fig2_table, fig2_stats = generate_figure02(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 02 saved: {fig2_img}")
        print(f"✔ Table saved:     {fig2_table}")
        print(f"✔ Stats saved:     {fig2_stats}")

    if 3 in RUN:
        print("Figure 03: DK-RD2 Coupling Structure")
        print("Generating figure... hang tight, this may take a while.")
        fig3_img, fig3_table, fig3_stats = generate_figure03(desi_bao_mean_path,desi_bao_cov_path)
        print(f"✔ Figure 03 saved: {fig3_img}")
        print(f"✔ Table saved:     {fig3_table}")
        print(f"✔ Stats saved:     {fig3_stats}")

    # ============================================================
    # Figure 04 — Effective Equation of State w_eff(z)
    # ============================================================
    if 4 in RUN:
        print("Figure 04: Explicit Diagnostic Reconstruction of w(z)")
        print("Generating figure... hang tight, this may take a while.")
        fig4_img, fig4_table, fig4_stats = generate_figure04(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 04 saved: {fig4_img}")
        print(f"✔ Table saved:     {fig4_table}")
        print(f"✔ Stats saved:     {fig4_stats}")

    # ------------------------------------------------------------------
    # Figure 05 — Stellar thermodynamic activity
    # ------------------------------------------------------------------
    if 5 in RUN:
        print("Figure 05: Stellar thermodynamic activity")
        print("Generating figure... hang tight, this may take a while.")
        print("A little more than a while...")
        fig5_img, fig5_table, fig5_stats = generate_figure05(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 05 saved: {fig5_img}")
        print(f"✔ Table saved:     {fig5_table}")
        print(f"✔ Stats saved:     {fig5_stats}")

    # ============================================================
    # Figure 06 — Effective Gravitational Interpretation of Core-Collapsed SIDM Halos in DK-RD2
    # ============================================================
    sidm_notebook_path = "data/data.nb" # from http://dx.doi.org/10.5281/zenodo.19116269
    if 6 in RUN:
        print("Figure 06: Effective Gravitational Interpretation of Core-Collapsed SIDM Halos in DK-RD2")
        print("External comparison data: Zenodo Mathematica notebook data.nb")
        print("data.nb download from http://dx.doi.org/10.5281/zenodo.19116269")
        print("Generating figure... hang tight, this may take a while.")
        fig6_img, fig6_table, fig6_stats, fig6_robustness = generate_figure06(
            sidm_notebook_path=sidm_notebook_path,
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 06 saved: {fig6_img}")
        print(f"✔ Table saved:     {fig6_table}")
        print(f"✔ Stats saved:     {fig6_stats}")
        print("Figure 06 robustness table:", fig6_robustness)

    if 7 in RUN:
        print("Figure 07: Hubble-Tension Projection Test")
        print("Generating figure... hang tight, this may take a while.")

        if "generate_figure07" not in globals():
            print("⚠ Figure 07 routine is not implemented yet.")
        else:
            fig7_img, fig7_table, fig7_stats = generate_figure07(
                desi_bao_mean_path=desi_bao_mean_path,
                desi_bao_cov_path=desi_bao_cov_path,
            )
            print(f"✔ Figure 07 saved: {fig7_img}")
            print(f"✔ Table saved:     {fig7_table}")
            print(f"✔ Stats saved:     {fig7_stats}")

    if 8 in RUN:
        print("Figure 08: Model Discrimination")
        print("Generating figure... hang tight, this may take a while.")

        if "generate_figure08" not in globals():
            print("⚠ Figure 08 routine is not implemented yet.")
        else:
            fig8_img, fig8_table, fig8_stats = generate_figure08(
                desi_bao_mean_path=desi_bao_mean_path,
                desi_bao_cov_path=desi_bao_cov_path,
            )
            print(f"✔ Figure 08 saved: {fig8_img}")
            print(f"✔ Table saved:     {fig8_table}")
            print(f"✔ Stats saved:     {fig8_stats}")

    print("\nDone. Falsifiable Predictions of the DK-RD2 figure generation completed.")
    banner = r"""
    ░▒▓████████▓▒░▒▓██████▓▒░░▒▓█▓▒░       ░▒▓███████▓▒░▒▓█▓▒░░▒▓██████▓▒░░▒▓███████▓▒░░▒▓█▓▒░      ░▒▓████████▓▒░
    ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░
    ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░
    ░▒▓██████▓▒░░▒▓████████▓▒░▒▓█▓▒░       ░▒▓██████▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓███████▓▒░░▒▓█▓▒░      ░▒▓██████▓▒░
    ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░             ░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░
    ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░             ░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░
    ░▒▓█▓▒░     ░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓███████▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓████████▓▒░▒▓████████▓▒░

    ░▒▓███████▓▒░░▒▓███████▓▒░░▒▓████████▓▒░▒▓███████▓▒░░▒▓█▓▒░░▒▓██████▓▒░▒▓████████▓▒░▒▓█▓▒░░▒▓██████▓▒░░▒▓███████▓▒░ ░▒▓███████▓▒░
    ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
    ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░        ░▒▓█▓▒░   ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
    ░▒▓███████▓▒░░▒▓███████▓▒░░▒▓██████▓▒░ ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░        ░▒▓█▓▒░   ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░
    ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░        ░▒▓█▓▒░   ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░
    ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ ░▒▓█▓▒░   ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░
    ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓███████▓▒░░▒▓█▓▒░░▒▓██████▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░
    """
    print(banner)
    print("— GabE=mc²+JzzOZ & Luludns -> ∞Ψ -")
    print(Core_autor_text)
    print(f"Output Evidence files in: '{Core_out_dir_path}/'")
