from brinkgrad import GradientGeneratorOptimizer
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, csv
from dolfinx import fem

targets = {
    'linear': lambda x: x[1]/500e-6,
    'sigmoid': lambda x: 1/(1+np.exp(-10*(x[1]/500e-6-0.5))),
    'double_peak': lambda x: np.sin(np.pi*x[1]/500e-6)**2,
    'staircase': lambda x: np.clip((x[1]//(500e-6/3)+1)/3,0,1)
}

def main():
    os.makedirs('figures', exist_ok=True)
    fig, axes = plt.subplots(2,2,figsize=(10,8))
    for ax,(name,target) in zip(axes.flat, targets.items()):
        opt = GradientGeneratorOptimizer(Lx=2000e-6, Ly=500e-6, nx=80, ny=20,
                                         target_expr=target, w_f=1e-3, w_c=5e1, V_star=0.5)
        # rho initialised internally by optimizer
        rho_phys = opt.run(max_iter=400, beta_continuation=[1,2,4,8,16], move=0.2)
        c_h = opt.c_h; x = opt.msh.geometry.x
        outlet_facets = opt.boundary_data["outlet"]
        dofs = fem.locate_dofs_topological(c_h.function_space, 1, outlet_facets)
        y_out = x[dofs,1]; c_out = c_h.x.array[dofs]
        idx = np.argsort(y_out); y_s, c_s = y_out[idx], c_out[idx]
        ax.plot(y_s*1e6,c_s,'b-',label='TopOpt')
        ax.plot(y_s*1e6,target(np.array([y_s*0,y_s])),'r--',label='target')
        ax.set_title(name); ax.set_xlabel('y (µm)'); ax.set_ylabel('c')
    fig.tight_layout(); fig.savefig('figures/gallery_profiles.pdf')
    print("Gallery saved to figures/gallery_profiles.pdf")
if __name__=='__main__': main()