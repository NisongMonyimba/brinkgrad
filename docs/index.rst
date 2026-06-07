brinkgrad
=========

An open-source FEniCSx framework for adjoint-based topology optimisation
of coupled Brinkman-convection-diffusion systems in porous media.

.. toctree::
   :maxdepth: 2

   api

Quick start
-----------

.. code-block:: python

   from brinkgrad import GradientGeneratorOptimizer

   opt = GradientGeneratorOptimizer(
       Lx=2000e-6, Ly=500e-6, nx=80, ny=20,
       target_expr=lambda x: x[1]/500e-6,
       w_f=1e-3, w_c=50.0, V_star=0.5)
   opt.run(max_iter=600)
   print(f"RMSE: {opt.rmse:.4f}")

Links
-----

* `GitHub <https://github.com/NisongMonyimba/brinkgrad>`_
* `Zenodo DOI <https://doi.org/10.5281/zenodo.20538833>`_
* `Paper <https://doi.org/10.5281/zenodo.20538833>`_
