#!/usr/bin/env python3
"""
examples/robust_projection_demo.py
Demonstrates three-field robust projection (Wang et al. 2011)
applied post-optimisation to improve design binarity.

Usage (inside Docker):
    python examples/robust_projection_demo.py
"""
import numpy as np
import sys; sys.path.insert(0, '.')
from brinkgrad import GradientGeneratorOptimizer
from brinkgrad.utilities import helmholtz_filter, heaviside_projection
from brinkgrad.solver import forward_solve
from brinkgrad.adjoint import adjoint_and_sensitivity
from brinkgrad.manufacturability import robust_projection
import brinkgrad.utilities as _ut
import dolfinx.fem as fem, ufl

Lx, Ly = 2000e-6, 500e-6
_ut.alpha_max = 1e3

print("Running 100-iteration optimisation (beta=4->16)...")
opt = GradientGeneratorOptimizer(
    Lx=Lx, Ly=Ly, nx=20, ny=5,
    target_expr=lambda x: x[1]/Ly,
    w_f=1e-3, w_c=5e1, V_star=0.5)
opt.run(max_iter=100, beta_continuation=[4, 8, 16], move=0.2)

# Nominal design at beta=16
helmholtz_filter(opt.rho, opt.rho_filt, opt.V_rho, opt.r_filter)
heaviside_projection(opt.rho_filt, opt.rho_phys, 16)
u0,_,c0 = forward_solve(opt.msh, opt.boundary_data, opt.rho_phys, P_in=1000.0)

ds_out = ufl.Measure("ds", domain=opt.msh,
                     subdomain_data=opt.boundary_data["facet_tag"],
                     subdomain_id=3)

def compute_metrics(rho_field, c_field, label):
    delta = fem.Function(fem.functionspace(opt.msh, ("Lagrange",1)))
    delta.interpolate(opt.target_expr)
    delta.x.array[:] -= c_field.x.array[:]; delta.x.scatter_forward()
    Jc = float(fem.assemble_scalar(fem.form(0.5*ufl.inner(delta,delta)*ds_out)))
    rho = rho_field.x.array
    gray   = float(np.mean((rho > 0.1) & (rho < 0.9)))
    binary = float(np.mean(rho <= 0.1) + np.mean(rho >= 0.9))
    rmse   = float(np.sqrt(2*Jc))
    print(f"  {label:12s}: RMSE={rmse:.4f}  gray={gray:.3f}  binary={binary:.3f}")
    return rmse, gray, binary

print("\n=== Results ===")
_, _, _ = compute_metrics(opt.rho_phys, c0, "nominal")

# Three-field robust projection (eta_d=0.3)
eroded, nominal_rp, dilated = robust_projection(
    opt.rho_filt, beta=16, eta_d=0.3, V_rho=opt.V_rho)

for field, label in [(eroded,"eroded"), (nominal_rp,"rp-nominal"), (dilated,"dilated")]:
    u2,_,c2 = forward_solve(opt.msh, opt.boundary_data, field, P_in=1000.0)
    compute_metrics(field, c2, label)

print("\nConclusion: robust_projection() is available in brinkgrad.manufacturability")
print("The eroded field provides the most binary-feasible design.")
print("For binary-convergent optimisation, call robust_projection() at each")
print("iteration and optimise against the worst-case (eroded) field.")
