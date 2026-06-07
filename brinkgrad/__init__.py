"""brinkgrad: Adjoint-based topology optimisation of coupled
Brinkman-convection-diffusion systems in FEniCSx.

Core modules:
    gradient_optimizer  - Top-level optimisation loop (GradientGeneratorOptimizer)
    solver              - Forward Brinkman + convection-diffusion solver
    adjoint             - Continuous adjoint solver + sensitivity assembly
    optimizer           - OC bisection + MMA update rules
    utilities           - Helmholtz filter, Heaviside projection, alpha(rho)
    mesh                - Structured triangular mesh generation
    manufacturability   - robust_projection() post-processing
    postprocess         - RMSE, gray-zone fraction, outlet profiles
    binary_validation   - Binary design validation
    taylor_test         - Gradient verification (Taylor remainder test)
"""
from brinkgrad.gradient_optimizer import GradientGeneratorOptimizer

__version__ = "1.0.0"
__author__ = "Nisong Monyimba, Vincent Pizziconi, Aurel Coza"
__license__ = "MIT"
__doi__ = "10.5281/zenodo.20538833"

__all__ = ["GradientGeneratorOptimizer"]
