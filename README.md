# DK-RD2-Falsifiable-Predictions
DK-RD2: Thermodynamic-relativistic cosmology framework with falsifiable tests beyond ΛCDM.

**Reference implementation for:**

**Falsifiable Predictions of the DK-RD2 Framework: Observational Tests Beyond ΛCDM**

Authors:

* Gabriel Martín del Campo Flores
* Jazmín Olvera Zacarías

---

## Overview

DK-RD2 (Duo-Kinetic Relativistic Dynamics) is a thermodynamic-relativistic cosmology framework designed to investigate falsifiable alternatives to the standard ΛCDM model.

The framework explores whether a temperature- and velocity-dependent effective gravitational coupling can reproduce key cosmological observations without introducing a fundamental cosmological constant.

The effective coupling is modeled as

$$
G_{ab}(T,v)
\approx
G_0
\left(
1+\alpha_{DK}
\frac{v^2}{c^2}
\frac{T_0}{T}
\right)
$$

where:

- $G_0$ is Newton's gravitational constant
- $T_0$ is the present-day CMB temperature
- $T$ is the effective thermal background
- $v$ is the characteristic relativistic velocity
- $\alpha_{DK}$ is the DESI-calibrated projection factor
The goal of this repository is not to prove DK-RD2, but to provide a fully reproducible observational falsification framework.

---

## Features

* DESI BAO calibration
* Cosmic Chronometer integration
* αDK likelihood reconstruction
* Effective dark-energy diagnostics
* Growth-history predictions (f\sigma_8(z))
* Compact-halo amplification tests
* Hubble-tension projection analysis
* Model discrimination against:

  * ΛCDM
  * Dynamic Dark Energy
  * Hot NEDE
  * DAO/DRMD

---

## Repository Structure

```text
.
├── Falsifiable_Predictions_DK_RD2.py
├── DK_RD2_Core.py
├── requirements.txt
├── data/
│   ├── DESI/
│   └── data.nb
└── evidence/
```

---

## Installation

Create a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Running

```bash
python "Falsifiable_Predictions_DK_RD2.py"
```

Generated figures, tables, and diagnostics are written to:

```text
./evidence/
```

---

## Main Outputs

The current pipeline generates:

1. αDK calibration
2. αDK likelihood constraints
3. Thermodynamic-relativistic coupling structure
4. Effective dark-energy reconstruction
5. Growth-history predictions
6. Compact-halo amplification tests
7. Hubble-tension projection diagnostics
8. Model discrimination tests

---

## Data Sources

The repository uses:

* DESI BAO DR2 public likelihood data
* Cosmic Chronometer measurements
* Published SIDM benchmark profiles (Figure 06)

The SIDM benchmark notebook used in Figure 06 is associated with:

DOI: 10.5281/zenodo.19116269

---

## Scientific Scope

The current implementation should be interpreted as a background-level and phenomenological falsification program.

A complete validation of DK-RD2 will require:

* A fully covariant formulation of the coupling
* First-principles perturbation dynamics
* Boltzmann-solver implementations

---

## Citation

If you use this code, please cite:

Martín del Campo Flores, G. (2026)

*Falsifiable Predictions of the DK-RD2 Framework: Observational Tests Beyond ΛCDM.*

DOI: 10.5281/zenodo.20637136

and

Olvera Zacarías, J.; Martín del Campo Flores, G. (2026)

*Falsifiable Predictions of the DK-RD2 Framework: Observational Tests Beyond ΛCDM.*

DOI: 10.5281/zenodo.20637136

---

## License

MIT License
