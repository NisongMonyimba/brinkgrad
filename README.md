# micrograd

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20479523.svg)](https://doi.org/10.5281/zenodo.20479523)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0000--7558--8580-green)](https://orcid.org/0009-0000-7558-8580)
[![CI](https://github.com/NisongMonyimba/micrograd/actions/workflows/ci.yml/badge.svg)](https://github.com/NisongMonyimba/micrograd/actions)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**An open-source FEniCSx framework for adjoint-based topology optimisation of porous microfluidic mixers.**

Nisong Monyimba · Vincent Pizziconi · Aurel Coza
School of Biological and Health Systems Engineering, Arizona State University, Tempe AZ USA

---

## What this does

`micrograd` finds a spatially-varying permeability field inside a microfluidic channel such that the outlet concentration profile matches a user-specified target. It solves a density-based topology optimisation problem on a coupled Brinkman–convection-diffusion system using:

- **Continuous adjoint** for gradient computation — two PDE solves regardless of mesh size
- **SUPG stabilisation** of both forward and adjoint equations at Pe ~ 10²–10⁵
- **Hybrid OC/MMA continuation schedule** with Helmholtz PDE filtering and Heaviside projection
- **End-to-end reproducibility** via Docker, GitHub Actions CI, and Zenodo archival

The primary benchmark achieves **RMSE = 0.058** for a linear gradient target on an 80×20 mesh in approximately one hour on a standard laptop (Intel Core i7, 16 GB RAM, single thread, no GPU).

> **Scope and limitation:** the optimised permeability field is a porous-medium surrogate operating within the Brinkman–SIMP model — it is **not a manufacturable binary design**. Hard-thresholding to binary permeability increases RMSE by +356% because the porous optimum exploits diffusive transport through intermediate-density regions that have no binary analogue. Users requiring binary geometries should apply robust projection ([Wang et al. 2011](https://doi.org/10.1007/s00158-010-0602-y)) during optimisation.

---

## Reproduce all results

```bash
git clone https://github.com/NisongMonyimba/micrograd
cd micrograd
docker run --rm -v ${PWD}:/root/micrograd \
  -e OMP_NUM_THREADS=1 -e MPLBACKEND=Agg \
  dolfinx/dolfinx:v0.7.3 bash -c \
  "cd /root/micrograd && bash run_all.sh"
```

This single command pulls the official FEniCSx Docker image, installs `micrograd`, and executes the complete pipeline. All figures and result files are written directly to the repository tree. Runtime: ~1 hour, no GPU required.

All results are permanently archived at Zenodo: **DOI: [10.5281/zenodo.20479523](https://doi.org/10.5281/zenodo.20479523)**

---

## Quick start

**Single optimisation run (linear gradient, 1400 iterations):**

```bash
docker run --rm -v ${PWD}:/root/micrograd \
  dolfinx/dolfinx:v0.7.3 bash -c \
  "cd /root/micrograd && python examples/linear_target.py"
```

**Interactive Jupyter notebook:**

```bash
docker run --rm -p 8888:8888 -v ${PWD}:/root/micrograd \
  dolfinx/dolfinx:v0.7.3 bash -c \
  "cd /root/micrograd && jupyter notebook --ip=0.0.0.0"
```

**Python API:**

```python
from micrograd import GradientGeneratorOptimizer

opt = GradientGeneratorOptimizer(
    target_expr=lambda x: x[1] / 500e-6,   # linear gradient
    Lx=2000e-6, Ly=500e-6,                  # 2 mm x 0.5 mm domain
    nx=80, ny=20                             # primary mesh
)
opt.run(max_iter=600)
print(f"RMSE = {opt.rmse:.4f}")
```

---

## Key results

| Schedule | Mesh | Iterations | RMSE | Gray fraction | Notes |
|----------|------|-----------|------|---------------|-------|
| Hybrid OC/MMA, β=128 | 80×20 | 1400 | **0.058** | 0.074 | Primary result |
| Hybrid OC/MMA, β=16 | 80×20 | 600 | 0.079 | 0.443 | Short run |
| Hybrid OC/MMA | 160×40 | 600 | 0.091 | 0.676 | Mesh sensitivity |
| Hybrid OC/MMA | 160×40 | 2400 | 0.107 | 0.557 | Multi-modal confirmed |
| Pure MMA | 80×20 | 600 | 0.304 | 0.412 | Validates hybrid schedule |

The hybrid OC/MMA schedule is **3.8× better RMSE** than pure MMA at the same iteration budget. Pure MMA stagnates near the initial objective because it cannot escape the uniform-density plateau at β=1; OC's bisection-based update provides the large initial displacement needed.

The optimisation is **not mesh-convergent**: finer meshes find different local minima (RMSE 0.091–0.107 on 160×40), a known characteristic of non-convex topology optimisation.

---

## Implementation contributions

### 1. Continuous adjoint for coupled Brinkman–convection-diffusion

The adjoint system is derived analytically, yielding:

- The coupling term **−λ∇c** in the flow adjoint equation
- The Dirichlet outlet BC **λ = c_target − c** derived directly from the misfit functional
- Explicit sensitivity expression usable by OC/MMA without solving Mĝ = g

This derivation is not previously published for this coupled system in FEniCSx.

### 2. SUPG-stabilised adjoint

The same element-wise stabilisation parameter τ is used for both forward and adjoint transport equations. Reduces outlet oscillation amplitude from 0.4 to < 10⁻³ and outlet slope change to < 2% across Pe ~ 10²–10⁵.

### 3. Modular continuation schedule

Five components, each addressing a diagnosed failure mode:

| Component | Failure mode | Effect |
|-----------|-------------|--------|
| Helmholtz PDE filter | Checkerboard instability | Smooth, mesh-independent ρ̃ |
| Heaviside β-continuation | Gray-zone lock-in | Drives ρ̄ → {0,1} gradually |
| α_max ramping | Early stagnation | Maintains OC step size |
| Hybrid OC/MMA | MMA plateau at β=1 | OC escapes; MMA stabilises |
| Sinusoidal initialisation | Symmetry trapping | Breaks left-right symmetry |

Gray fraction reduces from 0.443 (β=16, 600 iter) to 0.074 (β=128, 1400 iter).

### 4. End-to-end reproducibility at laptop scale

- Single Docker command reproduces all results
- GitHub Actions CI validates every commit inside `dolfinx/dolfinx:v0.7.3`
- Zenodo DOI for permanent archival
- Runtime: ~1 hour, no GPU required

---

## Why continuous adjoint instead of dolfin-adjoint?

The discrete adjoint via [dolfin-adjoint](https://www.dolfin-adjoint.org/) would provide exact gradient magnitudes automatically. The continuous adjoint was chosen for three reasons:

1. The derivation makes adjoint BCs **explicit**: λ = c_target − c on Γ_out, derived directly from the misfit functional.
2. The coupling term **−λ∇c** is physically interpretable as the sensitivity of concentration to velocity changes.
3. The L² Riesz representative is **directly usable by OC and MMA** without solving Mĝ = g, reducing implementation complexity.

The discrete adjoint and L-BFGS are identified as future work.

---

## Adjoint verification

| Test | Mesh | Result | Interpretation |
|------|------|--------|---------------|
| Taylor remainder (FD slope) | 80×20 | 2.07 | Fréchet differentiability confirmed |
| Pearson correlation with FD | 80×20 | 0.88 | Correct descent direction |

The large relative magnitude error (~3303%) is expected: the assembled sensitivity is the L² Riesz representative g, not the Euclidean gradient ĝ = M⁻¹g. Solving Mĝ = g with the consistent mass matrix would eliminate this error. For OC and MMA, only gradient direction matters.

---

## Repository structure

```
micrograd/
├── micrograd/                        # Core Python package
│   ├── solver.py                     # Brinkman + convection-diffusion forward solver
│   ├── adjoint.py                    # Continuous adjoint and sensitivity
│   ├── optimizer.py                  # OC and MMA updaters
│   ├── utilities.py                  # Helmholtz filter, Heaviside projection, RMSE
│   ├── gradient_optimizer.py         # Main topology optimisation class
│   ├── taylor_test.py                # Gradient verification (Taylor remainder)
│   ├── stabilization_validation.py   # SUPG verification
│   ├── convergence_study.py          # Mesh convergence study
│   ├── binary_validation.py          # Hard-thresholding and binary gap
│   ├── manufacturability.py          # Minimum feature size checks
│   ├── postprocess.py                # Figures and result export
│   ├── scalability.py                # Timing and scaling benchmarks
│   ├── validation.py                 # Navier-Stokes validation (optional)
│   ├── christmas_tree.py             # Christmas-tree reference generator
│   ├── multiobjective.py             # Multi-objective sweep utilities
│   ├── uncertainty_quantification.py # PCE-based UQ (optional)
│   ├── filter_sensitivity.py         # Filter radius sensitivity
│   ├── experimental_metrics.py       # Pe, Re, Da computation
│   ├── mesh.py                       # Mesh generation helpers
│   ├── compatibility.py              # FEniCSx version compatibility
│   └── __init__.py
├── examples/
│   ├── linear_target.py              # Linear gradient (primary benchmark)
│   ├── double_peak_target.py         # Double-peak concentration target
│   ├── gallery_targets.py            # Multiple target profiles
│   └── run_convergence_study.py      # 4-mesh convergence sweep
├── tests/
│   ├── test_import.py                # Import and smoke tests
│   └── __init__.py
├── manuscript/                       # LaTeX source (submit-ready)
│   ├── main.tex                      # Master document
│   ├── macros.tex                    # All numerical results as macros
│   ├── references.bib
│   ├── abstract.tex
│   ├── chapter1_introduction.tex
│   ├── chapter2_mathematical_model.tex
│   ├── chapter3_numerical_methods.tex
│   ├── chapter4_results.tex
│   ├── chapter5_discussion.tex
│   ├── chapter6_conclusion.tex
│   ├── chapter7_data_availability.tex
│   ├── chapter8_appendices.tex
│   ├── cover_letter.tex
│   └── figures/                      # Manuscript figures (PDF)
├── figures/                          # Generated output figures (PDF + PNG)
├── docs/
│   ├── si/                           # docs/si/ — supplementary figures (PDF)
│   ├── api.rst                       # API documentation source
│   └── conf.py                       # Sphinx config
├── .github/
│   └── workflows/ci.yml              # GitHub Actions CI pipeline
├── setup.py                          # pip-installable package
├── environment.yaml                  # Conda environment spec
├── run_all.sh                        # Full reproduction script
├── submission_checklist.sh           # Pre-submission manuscript checker
├── generate_macros.py                # Auto-generate macros.tex from results
├── CITATION.cff                      # Machine-readable citation
├── LICENSE                           # MIT
└── README.md
```

---

## Manuscript

**Title:** `micrograd`: An open-source FEniCSx framework for adjoint-based topology optimisation of porous microfluidic mixers using Brinkman–convection-diffusion equations

**Authors:** Nisong Monyimba, Vincent Pizziconi, Aurel Coza

**Target journal:** *Engineering with Computers* (Springer)

**Preprint:** arXiv cs.NA *(submitted)*

```bibtex
@software{micrograd2025zenodo,
  author    = {Monyimba, Nisong and Pizziconi, Vincent and Coza, Aurel},
  title     = {{micrograd: Topology optimisation of microfluidic
                concentration gradient generators (v1.0.0)}},
  year      = {2026},
  publisher = {Zenodo},
  version   = {1.0.0},
  doi       = {10.5281/zenodo.20479523},
  url       = {https://github.com/NisongMonyimba/micrograd},
  note      = {Zenodo: 10.5281/zenodo.20479523}
}
```

---

## Requirements

**Recommended:** [Docker Desktop](https://www.docker.com/products/docker-desktop/) — pulls `dolfinx/dolfinx:v0.7.3` automatically, no other setup needed.

**Native install:** FEniCSx 0.7.3, PETSc, petsc4py, numpy, scipy, matplotlib.

---

## License

MIT — see [LICENSE](LICENSE).

---

## Contact

**Nisong Monyimba**
School of Biological and Health Systems Engineering, Arizona State University
Email: nisongmonyimba278@gmail.com
ORCID: [0009-0000-7558-8580](https://orcid.org/0009-0000-7558-8580)
GitHub: [NisongMonyimba](https://github.com/NisongMonyimba)
