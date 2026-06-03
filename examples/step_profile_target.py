#!/usr/bin/env python3
"""
examples/step_profile_target.py
Demonstrates software flexibility: only the target_expr changes.
Step profile: c*(y) = 1 if y > Ly/2 else 0

This example, together with linear_target.py and double_peak_target.py,
shows that micrograd supports arbitrary concentration targets by
changing a single parameter (target_expr).

Usage (inside Docker):
    python examples/step_profile_target.py
"""
import numpy as np
import sys; sys.path.insert(0, '.')
from micrograd import GradientGeneratorOptimizer
import micrograd.utilities as _ut

Lx, Ly = 2000e-6, 500e-6
_ut.alpha_max = 1e3

# Step profile: top half = 1, bottom half = 0
step_expr = lambda x: np.where(x[1] > Ly/2, 1.0, 0.0)

print("Running step profile optimisation (target: c* = 1 for y > Ly/2)...")
opt = GradientGeneratorOptimizer(
    Lx=Lx, Ly=Ly, nx=20, ny=5,
    target_expr=step_expr,
    w_f=1e-3, w_c=5e1, V_star=0.5)
opt.run(max_iter=200, beta_continuation=[4, 8, 16, 32], move=0.2)

print("\nStep profile optimisation complete.")
print("Only target_expr changed vs linear_target.py.")
print("All other solver, adjoint, and optimiser code is identical.")
print("This demonstrates the software flexibility of micrograd.")
