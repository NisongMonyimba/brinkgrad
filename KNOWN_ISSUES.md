# Known Issues and Limitations

This file documents known limitations of `brinkgrad` v1.0.0.
See also Section 6 (Discussion) of the manuscript for full details.

---

## L1: Binary gap (+356% RMSE after hard-thresholding)

**What happens:** The optimised permeability field is a porous-medium surrogate.
Hard-thresholding to binary solid/fluid increases RMSE by +356%.

**Root cause:** The Brinkman–SIMP model allows diffusive transport through
intermediate-density elements that have no physical binary analogue.
The porous optimum exploits this gray-zone mixing; it is lost upon thresholding.

**Workaround:** Apply robust projection ([Wang et al. 2011](https://doi.org/10.1007/s00158-010-0602-y))
during optimisation to penalise intermediate densities continuously.
This is identified as future work.

---

## L2: Mesh non-convergence

**What happens:** Refining the mesh from 80×20 to 160×40 gives worse RMSE
(0.091 at 600 iter, 0.107 at 2400 iter vs 0.079 at 600 iter on 80×20).

**Root cause:** Topology optimisation is non-convex. Finer meshes find
different local minima, not better ones. This is a known property of
density-based methods, not a bug.

**Implication:** Results are mesh-dependent. The 80×20 mesh is the
validated primary configuration. Do not interpret finer-mesh results
as more accurate.

---

## L3: Continuous adjoint magnitude error (~3303%)

**What happens:** The assembled sensitivity g_i = ∫(∂J/∂ρ̄)ψ_i dx is the
L² Riesz representative, not the Euclidean gradient ĝ = M⁻¹g.
The lumped approximation g/m_i gives large relative errors vs finite differences.

**Impact:** None for OC and MMA — both depend only on gradient sign/direction,
confirmed by Pearson correlation ≥ 0.80.

**Workaround:** Solve Mĝ = g with the consistent mass matrix, or use
dolfin-adjoint for exact discrete gradients. Needed only for L-BFGS.

---

## L4: OC/MMA switch point is empirical

**What happens:** The switch from OC (iterations 0–299) to MMA (300+)
is validated by comparison (pure MMA gives RMSE=0.304 vs hybrid 0.079)
but the switch point is not theoretically proven optimal.

**Implication:** Different problems may benefit from different switch points.
A systematic schedule search is future work.

---

## L5: 2D only, Brinkman surrogate only

All results are for 2D Brinkman–convection-diffusion.
3D extrusion and full Navier–Stokes validation are not implemented
in the current release. Body-fitted NS validation code exists in
`brinkgrad/validation.py` but requires `pyvista`/`gmsh` not in
the standard Docker image.

---

## L6: DOLFINx version pinned to v0.7.3

The code is tested with `dolfinx/dolfinx:v0.7.3`.
Newer versions (v0.8.x+) are expected to work but have not been
systematically tested. If you encounter issues with newer versions,
please open a GitHub issue.

---

*For questions or bug reports, open an issue at
https://github.com/NisongMonyimba/brinkgrad/issues*
