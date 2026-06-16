# brinkgrad

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20538833.svg)](https://doi.org/10.5281/zenodo.20538833)
[![CI](https://github.com/NisongMonyimba/brinkgrad/actions/workflows/ci.yml/badge.svg)](https://github.com/NisongMonyimba/brinkgrad/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0000--7558--8580-green)](https://orcid.org/0009-0000-7558-8580)

**An open-source FEniCSx framework for adjoint-based topology optimisation of coupled Brinkman–convection-diffusion systems in porous media.**

Nisong Monyimba · Vincent Pizziconi · Aurel Coza  
School of Biological and Health Systems Engineering, Arizona State University, Tempe AZ USA

> **Contribution:** `brinkgrad` provides an openly documented, CI-verified, and Docker-reproducible FEniCSx framework for coupled Brinkman–convection–diffusion topology optimisation with SUPG-stabilised adjoint transport.

---

## What brinkgrad does

`brinkgrad` finds a spatially-varying permeability field inside a microfluidic channel such that the outlet concentration profile matches a user-specified target. It solves a density-based topology optimisation problem on a coupled Brinkman–convection-diffusion system using:

- **Continuous adjoint** — gradient computation in two PDE solves regardless of mesh size
- **SUPG stabilisation** of both forward and adjoint equations at Pe ~ 10²–10⁵
- **OC optimiser** with Helmholtz PDE filtering and Heaviside β-continuation
- **End-to-end reproducibility** via Docker, GitHub Actions CI, and Zenodo archival

**Scope:** `brinkgrad` v1.0 targets 2D density optimisation on structured meshes within the Brinkman–SIMP surrogate model. Post-processing manufacturability analysis is available via `robust_projection()`.

---

## Reproduce all results (one command)

```bash
git clone https://github.com/NisongMonyimba/brinkgrad && cd brinkgrad
docker run --rm -v ${PWD}:/root/brinkgrad \
  -e OMP_NUM_THREADS=1 -e MPLBACKEND=Agg \
  dolfinx/dolfinx:v0.7.3 bash -c \
  "cd /root/brinkgrad && pip install -e . --quiet && bash run_all.sh"
```

Expected output: RMSE = 0.058 after ~1 hour on 80×20 mesh (standard laptop, no GPU).  
All results permanently archived: **DOI: [10.5281/zenodo.20538833](https://doi.org/10.5281/zenodo.20538833)**

---

## Quick start

```python
from brinkgrad import GradientGeneratorOptimizer

opt = GradientGeneratorOptimizer(
    Lx=2000e-6, Ly=500e-6, nx=80, ny=20,
    target_expr=lambda x: x[1]/500e-6,   # linear gradient target
    w_f=1e-3, w_c=50.0, V_star=0.5)
opt.run(max_iter=600)
print(f"RMSE: {opt.rmse:.4f}")
```

**Change the target profile** by modifying one argument:

```python
import numpy as np
step = lambda x: np.where(x[1] > 250e-6, 1.0, 0.0)
opt = GradientGeneratorOptimizer(..., target_expr=step, ...)
```

See `examples/` for seven ready-to-run scripts.

---

## Installation

```bash
# Inside Docker (recommended for reproducibility)
docker run --rm -it dolfinx/dolfinx:v0.7.3 bash
pip install git+https://github.com/NisongMonyimba/brinkgrad.git

# Local (requires FEniCSx v0.7.3)
git clone https://github.com/NisongMonyimba/brinkgrad
cd brinkgrad && pip install -e .
```

---

## Package structure
brinkgrad/
├── adjoint.py          # Continuous adjoint solver + sensitivity assembly
├── solver.py           # Forward Brinkman + convection-diffusion solver
├── optimizer.py        # OC bisection + MMA (nlopt) update rules
├── gradient_optimizer.py  # Top-level optimisation loop
├── utilities.py        # Helmholtz filter, Heaviside projection, alpha(ρ)
├── mesh.py             # Structured triangular mesh generation
├── manufacturability.py   # robust_projection() post-processing
├── postprocess.py      # RMSE, gray-zone fraction, outlet profiles
├── binary_validation.py   # Binary design validation (threshold study)
├── taylor_test.py      # Gradient verification (Taylor remainder test)
└── experimental/       # Auxiliary modules (not part of core API)
examples/
├── linear_target.py           # Primary benchmark (linear gradient)
├── step_profile_target.py     # Step profile target
├── double_peak_target.py      # Double-peak target
├── robust_projection_demo.py  # Binary design post-processing
├── gallery_targets.py         # All targets in one script
├── run_convergence_study.py   # Mesh/iteration sensitivity
└── christmas_tree_comparison.py  # Comparison run
tests/
├── test_import.py   # Import integrity (8 tests)
└── test_core.py     # Unit tests for core functions (10 tests)
notebooks/
└── quickstart.ipynb  # Interactive tutorial

---

## Verification

| Layer | Test | Result |
|-------|------|--------|
| Code health | 18 automated tests (CI) | All pass on every push |
| Taylor remainder | Slope ≈ 1.0 (Riesz representative) | Pass |
| Pipeline | Mini-optimisation (20 iter) | J decreases monotonically |
| Benchmark | Linear gradient, 80×20 mesh | RMSE = 0.058 |

---

## Reproducibility

| Item | Detail |
|------|--------|
| Docker image | `dolfinx/dolfinx:v0.7.3` (pinned) |
| CI | GitHub Actions, every push |
| Archive | Zenodo DOI: 10.5281/zenodo.20538833 |
| Licence | MIT |
| Single-command reproduction | `bash run_all.sh` inside Docker |

---

## Citation

```bibtex
@software{brinkgrad2026,
  author  = {Monyimba, Nisong and Pizziconi, Vincent and Coza, Aurel},
  title   = {brinkgrad: An open-source FEniCSx framework for adjoint-based
             topology optimisation of coupled flow-transport in porous media},
  year    = {2026},
  doi     = {10.5281/zenodo.20538833},
  url     = {https://github.com/NisongMonyimba/brinkgrad},
  license = {MIT}
}
```

See also `CITATION.cff` for machine-readable citation metadata.
