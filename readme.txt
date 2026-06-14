______    _     _  __ _       _     _       ______             _ _      _   _                 
|  ___|  | |   (_)/ _(_)     | |   | |      | ___ \           | (_)    | | (_)                
| |_ __ _| |___ _| |_ _  __ _| |__ | | ___  | |_/ / __ ___  __| |_  ___| |_ _  ___  _ __  ___ 
|  _/ _` | / __| |  _| |/ _` | '_ \| |/ _ \ |  __/ '__/ _ \/ _` | |/ __| __| |/ _ \| '_ \/ __|
| || (_| | \__ \ | | | | (_| | |_) | |  __/ | |  | | |  __/ (_| | | (__| |_| | (_) | | | \__ \
\_| \__,_|_|___/_|_| |_|\__,_|_.__/|_|\___| \_|  |_|  \___|\__,_|_|\___|\__|_|\___/|_| |_|___/
                                                                                              
                                                                                              
______ _   __     ____________  _____                                                         
|  _  \ | / /     | ___ \  _  \/ __  \                                                        
| | | | |/ /______| |_/ / | | |`' / /'                                                        
| | | |    \______|    /| | | |  / /                                                          
| |/ /| |\  \     | |\ \| |/ / ./ /___                                                        
|___/ \_| \_/     \_| \_|___/  \_____/                                                        
                                                                                              
                                                                                              

DK-RD2 Falsifiable Predictions Pipeline
=======================================
Release Date 6-7 June 2026

Title
-----
Falsifiable Predictions of the DK-RD2 Framework:
Observational Tests Beyond ΛCDM

Authors
-------
Gabriel Martín del Campo Flores
Jazmín Olvera Zacarías

Contact
-------
Gabriel Martín del Campo Flores
Email: gabemdelc@gmail.com

GitHub https://github.com/gabemdelc/DK-RD2-Falsifiable-Predictions
Zenodo: https://doi.org/10.5281/zenodo.20637136

Overview
--------
This package contains the reproducible Python pipeline used to generate in /evidence directory the
figures, tables, and statistical diagnostics for the manuscript (the :

    Falsifiable Predictions of the DK-RD2 Framework:
    Observational Tests Beyond ΛCDM

DK-RD2 is a thermodynamic-relativistic framework in which the effective
gravitational coupling is modeled as

    Gab(T, v) ≈ G0 · [1 + α_DK · (v²/c²) · (T0/T)]

where α_DK is the geometric projection factor calibrated from DESI BAO +
Cosmic Chronometer data. For cosmological applications the thermal background
is approximated as

    T(z) ≈ T0 · (1 + z)

The purpose of this code is not to prove DK-RD2 a priori, but to provide a
reproducible falsification pipeline: every figure is designed to expose a
specific observable that can support, constrain, or rule out the framework.


Package files
-------------
The minimal runnable package should contain:

    Falsifiable_ Predictions_DK-RD2.py
        Main driver script for the June 2026 falsifiable-predictions paper.
        It loads data, performs the global α_DK calibration, generates figures,
        and writes all tables/statistics under ./evidence/.

    DK_RD2_Core.py
        Core computational engine. It includes constants, Gab(T,v), DK-RD2 and
        ΛCDM expansion functions, distances, BAO observables, likelihood tools,
        growth/lensing helpers, and CSV output utilities.

    readme.txt or README_Falsifiable_Predictions_DK_RD2.txt
        This installation and usage guide.

Recommended optional files:

    requirements.txt
        Python dependency list, if you create one from the install command below.

    LICENSE
        MIT License text.


Required directory structure
----------------------------
Run the code from the project root directory. The expected structure is:

    project_root/
    ├── Falsifiable_ Predictions_DK-RD2.py
    ├── DK_RD2_Core.py
    ├── readme.txt
    ├── data/
    │   ├── DESI/
    │   │   └── bao_data/
    │   │       └── desi_bao_dr2/
    │   │           ├── desi_gaussian_bao_ALL_GCcomb_mean.txt
    │   │           └── desi_gaussian_bao_ALL_GCcomb_cov.txt
    │   └── data.nb
    └── evidence/

The script creates ./evidence/ automatically if it does not already exist.

Important data notes
--------------------
The DESI BAO files are required for the global α_DK and r_d calibration. The
current driver computes this global calibration before generating any selected
figure, so the following files must exist:

    data/DESI/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt
    data/DESI/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt

Figure 06 additionally requires the external SIDM benchmark file:

    data/data.nb

The code comments identify this file as the Mathematica notebook benchmark from
Zenodo:

    http://dx.doi.org/10.5281/zenodo.19116269

If data.nb is missing, Figure 06 will fail, but the other figures can still be
run by selecting only the desired figure numbers.


Installation
------------
Python 3.12 or newer is recommended.

Linux / macOS
~~~~~~~~~~~~~
From inside the project directory:

    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install --upgrade pip
    python -m pip install numpy scipy matplotlib pandas astropy openpyxl

Windows PowerShell
~~~~~~~~~~~~~~~~~~
From inside the project directory:

    py -3.12 -m venv .venv
    .\.venv\Scripts\Activate.ps1
    python -m pip install --upgrade pip
    python -m pip install numpy scipy matplotlib pandas astropy openpyxl

Minimal dependencies
~~~~~~~~~~~~~~~~~~~~
Required for the main pipeline:

    numpy
    scipy
    matplotlib
    pandas

Recommended / optional:

    astropy    Optional FITS support used by auxiliary SPHEREx helpers in the core.
    openpyxl   Optional Excel support for generic table loaders.

CLASS is not required by the current Falsifiable Predictions driver. The current
paper pipeline is focused on background-level observables, BAO calibration,
statistical likelihoods, effective equation-of-state diagnostics, growth,
compact-halo tests, Hubble projection, and model discrimination.


Running the code
----------------
Because the main script filename contains spaces, run it with quotes:

    python "Falsifiable_ Predictions_DK-RD2.py"

The script will:

1. Create ./evidence/ if needed.
2. Print the DK-RD2 banner.
3. Print the available figure IDs.
4. Ask for a selector.
5. Compute the official global DK-RD2 calibration:

       α_DK from DESI BAO + Cosmic Chronometers
       r_d from the same joint likelihood

6. Generate the selected figures and CSV files.

Selector examples
-----------------
At the prompt, use:

    ENTER       Keep default selector, usually all.
    0           Run all figures.
    all         Run all figures.
    1           Run only Figure 01.
    1,2         Run Figures 01 and 02.
    1-4         Run Figures 01 through 04.
    1-4,7,8     Run Figures 01, 02, 03, 04, 07, and 08.
    X           Exit.

If you prefer a filename without spaces, you may copy the script:

    cp "Falsifiable_ Predictions_DK-RD2.py" Falsifiable_Predictions_DK_RD2.py
    python Falsifiable_Predictions_DK_RD2.py


Figure overview
---------------
The current driver generates eight main figures for the falsifiable-predictions
paper.

Figure 01 — Statistical calibration of α_DK and r_d
    Calibrates the DK-RD2 projection factor α_DK and the BAO sound horizon r_d
    using DESI DR2 Gaussian BAO data plus a Cosmic Chronometer anchor. It also
    compares DK-RD2 and ΛCDM distance observables.

Figure 02 — Likelihood constraint on α_DK
    Profiles Δχ²(α_DK) while independently optimizing r_d for each α_DK value.
    This tests whether α_DK is observationally constrained rather than arbitrary.

Figure 03 — Thermodynamic-relativistic coupling structure
    Visualizes the maximal and projected DK-RD2 gravitational coupling,
    including the DESI-calibrated projection factor and the bounded saturation
    behavior of the coupling.

Figure 04 — Effective dark-energy diagnostic reconstruction
    Reconstructs Ω_DE,eff(z), w_eff(z), and Δw(z) from the DK-RD2 expansion
    history. This is a diagnostic FRW reinterpretation only; DK-RD2 does not
    introduce a fundamental dark-energy fluid.

Figure 05 — Growth of cosmic structure / fσ8(z)
    Computes the DK-RD2 growth prediction from the same globally calibrated
    background and compares it against ΛCDM-style growth behavior.

Figure 06 — Compact-halo prediction and effective gravitational amplification
    Uses published SIDM benchmark profiles from data.nb as an external
    comparison and tests whether compact gravitational behavior can also be
    represented through a bounded DK-RD2 effective-gravity amplification.

Figure 07 — Hubble-tension projection test
    Evaluates whether thermodynamic projection of the effective coupling can
    generate distinct early- and late-Universe H0 reconstructions.

Figure 08 — Model discrimination versus DESI alternatives
    Compares DK-RD2 against alternative explanations of DESI-scale deviations,
    including evolving dark energy, Hot NEDE, and DAO/DRMD scenarios.


Outputs
-------
All generated outputs are written to ./evidence/.

Typical outputs include:

    Falsifiable_ Predictions_DK-RD2_image_XX.png
    Falsifiable_ Predictions_DK-RD2_table_XX.csv
    Falsifiable_ Predictions_DK-RD2_table_XX_stats.csv

Exact filenames are generated by DK_RD2_Core.py through the generate_evidence()
helper and may vary depending on the configured output naming convention.

The CSV outputs include numerical tables, residuals, likelihood summaries,
calibration metadata, χ² diagnostics, and robustness checks where applicable.


Scientific interpretation
-------------------------
This code evaluates DK-RD2 as a falsifiable observational framework.

The central calibrated relation is:

    Gab(T, v) ≈ G0 · [1 + α_DK · (v²/c²) · (T0/T)]

with α_DK calibrated from DESI BAO + Cosmic Chronometer data and reused across
all figures. The framework is challenged if future observations show any of the
following:

    α_DK → 0 within uncertainties.
    The DESI BAO + CC likelihood minimum disappears.
    The reconstructed expansion remains fully ΛCDM-compatible.
    The predicted growth history fσ8(z) is not observed.
    Compact-halo amplification signatures are absent.
    Early- and late-Universe H0 reconstructions converge without DK-RD2-like
    thermodynamic projection effects.
    Future lensing/growth/expansion correlations favor dark-sector extensions
    rather than DK-RD2 thermodynamic-gravity signatures.


Limitations and scope
---------------------
The current pipeline should be interpreted as a background-level and
phenomenological falsification program. A complete validation of DK-RD2 will
require a fully covariant formulation of the coupling and a first-principles
implementation of modified perturbation dynamics in a Boltzmann solver.

The code is intended to make the current predictions reproducible and testable,
not to claim final theoretical closure.


Troubleshooting
---------------
1. FileNotFoundError for DESI BAO files
   Check that the following files exist relative to the project root:

       data/DESI/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt
       data/DESI/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt

2. Figure 06 fails with missing data.nb
   Download or place the SIDM benchmark notebook at:

       data/data.nb

   Or rerun the script and select figures excluding 6, for example:

       1-5,7,8

3. ImportError for numpy/scipy/matplotlib/pandas
   Activate your virtual environment and reinstall dependencies:

       python -m pip install numpy scipy matplotlib pandas astropy openpyxl

4. Problems caused by the script filename containing spaces
   Run with quotes:

       python "Falsifiable_ Predictions_DK-RD2.py"

   Or copy it to a simpler local filename:

       cp "Falsifiable_ Predictions_DK-RD2.py" Falsifiable_Predictions_DK_RD2.py
       python Falsifiable_Predictions_DK_RD2.py


Citation / credit
-----------------
If you use this code or its generated figures, please cite the associated DK-RD2
manuscripts and include credit to the authors.

Primary framework:

    Martín del Campo Flores, G. (2026).
    What Is Gravity? From Emergent Mass to Thermodynamic-Geometric Dynamics
    in the DK-RD2 Framework.
    Zenodo. https://doi.org/10.5281/zenodo.20078175

Falsifiable-predictions manuscript:

    Olvera Zacarías, J.; Martín del Campo Flores, G. (2026).
    Falsifiable Predictions of the DK-RD2 Framework:
    Observational Tests Beyond ΛCDM.
    GitHub https://github.com/gabemdelc/DK-RD2-Falsifiable-Predictions
    Zenodo: https://doi.org/10.5281/zenodo.20637136    


License
-------
MIT License.

You are free to use, modify, and distribute the code for research and educational
purposes, provided that proper credit is given and the license terms are respected.


Motto
-----
Gravity is not assumed — it emerges.
GabE=mc²+JzzOZ & Luludns => ∞Ψ
