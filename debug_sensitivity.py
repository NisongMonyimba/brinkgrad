import numpy as np
from brinkgrad.gradient_optimizer import GradientGeneratorOptimizer

# Match the call from linear_target.py (adjust keyword if necessary)
opt = GradientGeneratorOptimizer(Lx=2000e-6, Ly=500e-6, target='linear')
# Run one iteration but manually inspect
opt.run(max_iter=1)
# Get the last sensitivity
sens = opt.sensitivity  # assuming the optimiser stores it
print(f"Sensitivity min: {sens.min():.6e}, max: {sens.max():.6e}, mean: {sens.mean():.6e}")
print(f"Design change max: {np.abs(opt.rho.x.array - opt.rho_old.x.array).max():.6e}")
