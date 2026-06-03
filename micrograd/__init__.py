# Compatibility shim — package renamed to brinkgrad
# This directory is retained for backward compatibility
# Import from brinkgrad instead
import warnings
warnings.warn('micrograd is renamed brinkgrad; update your imports', DeprecationWarning)
from brinkgrad import *
from brinkgrad import GradientGeneratorOptimizer
