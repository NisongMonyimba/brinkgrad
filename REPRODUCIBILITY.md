# brinkgrad Reproducibility Checklist

## One-click reproduction (Docker)
```bash
git clone https://github.com/NisongMonyimba/brinkgrad.git
cd brinkgrad
docker run --rm -v $(pwd):/shared dolfinx/dolfinx:v0.7.3 \
    bash -c "cd /shared && pip install -e . --no-deps -q && bash run_all.sh"
```
All results are written to `results/` as JSON files.

## Environment
| Component | Version |
|-----------|---------|
| Docker image | `dolfinx/dolfinx:v0.7.3` (pinned) |
| FEniCSx | 0.7.3 |
| Python | 3.11 |
| PETSc | 3.20 |
| numpy | ≥1.24 |

## Software metrics
| Metric | Value |
|--------|-------|
| Lines of code | 1,892 (non-empty, non-comment) |
| Public functions/classes | 68 across 9 modules |
| Automated tests | 15 (8 unit + 7 CI integration) |
| Test coverage | 8.7% (import-level; full pipeline tested via CI) |
| Peak memory (20×5, 50 iter) | 29.5 MB |
| Runtime (80×20, 1 iter) | ~1 s (8-core laptop) |
| Example scripts | 6 |

## Archived outputs
- Zenodo DOI: [10.5281/zenodo.20538833](https://doi.org/10.5281/zenodo.20538833)
- All result JSON files included in archive

## CI badge
![CI](https://github.com/NisongMonyimba/brinkgrad/actions/workflows/ci.yml/badge.svg)

## Verification
The primary result (RMSE = 0.058, 80×20 mesh, 1400 iterations) is
reproduced by `bash run_all.sh` to floating-point precision within
the pinned Docker environment.
