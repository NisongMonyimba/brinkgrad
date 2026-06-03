"""
Publication-ready topology figure for the microfluidic gradient generator.
Produces: figures/topology_optimised.pdf  (and .png at 300 dpi)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dolfinx import fem
from brinkgrad import GradientGeneratorOptimizer
from brinkgrad.utilities import helmholtz_filter, heaviside_projection

# ── 1. Run optimiser ────────────────────────────────────────────────────────
Lx, Ly = 2000e-6, 500e-6
opt = GradientGeneratorOptimizer(
    Lx=Lx, Ly=Ly, nx=20, ny=5,
    target_expr=lambda x: x[1] / Ly,
    w_f=1e-7, w_c=5e4, V_star=0.5)
opt.run(max_iter=400, beta_continuation=[1,2,4,8,16,32,64], move=0.05)

# ── 2. Extract density on mesh nodes ───────────────────────────────────────
coords   = opt.msh.geometry.x          # (n_nodes, 3)
rho_vals = opt.rho_phys.x.array.copy()  # physical density [0,1]
x_nodes  = coords[:, 0] * 1e6          # µm
y_nodes  = coords[:, 1] * 1e6

# ── 3. Build structured grid via scatter interpolation ─────────────────────
Nx, Ny = 200, 50                        # interpolation grid
xi = np.linspace(0, Lx*1e6, Nx)
yi = np.linspace(0, Ly*1e6, Ny)
Xi, Yi = np.meshgrid(xi, yi)

from scipy.interpolate import griddata
Zi = griddata((x_nodes, y_nodes), rho_vals,
              (Xi, Yi), method='linear')
# Fill NaNs at boundary
from scipy.ndimage import generic_filter
mask = np.isnan(Zi)
if mask.any():
    filled = generic_filter(np.where(mask, 0, Zi),
                            np.nanmean, size=5)
    Zi = np.where(mask, filled, Zi)
Zi = np.clip(Zi, 0, 1)

# ── 4. Publication figure ───────────────────────────────────────────────────
# Colour map: white = fluid, black = solid (standard in topology opt)
cmap_topo = mcolors.LinearSegmentedColormap.from_list(
    'topopt', ['#f5f0e8', '#2a2118'], N=256)

fig, axes = plt.subplots(
    2, 1, figsize=(7.0, 3.8),
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.35})

# ── Panel A: density heatmap ────────────────────────────────────────────────
ax = axes[0]
im = ax.imshow(Zi, origin='lower', aspect='auto',
               extent=[0, Lx*1e6, 0, Ly*1e6],
               cmap=cmap_topo, vmin=0, vmax=1,
               interpolation='bilinear')

# Binary threshold contour at rho=0.5
ax.contour(Xi, Yi, Zi, levels=[0.5],
           colors=['#e85d26'], linewidths=0.8, linestyles='-')

cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label(r'Density $\rho$ (0 = fluid, 1 = solid)',
               fontsize=7, labelpad=4)
cbar.ax.tick_params(labelsize=6)

ax.set_xlabel(r'$x$ (µm)', fontsize=8)
ax.set_ylabel(r'$y$ (µm)', fontsize=8)
ax.tick_params(labelsize=7)
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(50))
ax.set_title('(a) Optimised topology — linear concentration gradient',
             fontsize=8, loc='left', pad=4)

# Mark inlets / outlet
ax.annotate('Inlet 1\n(c = 0)', xy=(0, Ly*1e6*0.25),
            xytext=(-320, Ly*1e6*0.25),
            fontsize=6, color='#1a6fa8',
            arrowprops=dict(arrowstyle='->', color='#1a6fa8', lw=0.8),
            annotation_clip=False)
ax.annotate('Inlet 2\n(c = 1)', xy=(0, Ly*1e6*0.75),
            xytext=(-320, Ly*1e6*0.75),
            fontsize=6, color='#c0392b',
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=0.8),
            annotation_clip=False)
ax.annotate('Outlet', xy=(Lx*1e6, Ly*1e6*0.5),
            xytext=(Lx*1e6+80, Ly*1e6*0.5),
            fontsize=6, color='#27ae60',
            arrowprops=dict(arrowstyle='->', color='#27ae60', lw=0.8),
            annotation_clip=False)

# ── Panel B: outlet concentration profile ───────────────────────────────────
ax2 = axes[1]
out_facets = opt.boundary_data["outlet"]
dofs = fem.locate_dofs_topological(
    opt.c_h.function_space, opt.msh.topology.dim-1, out_facets)
y_out = opt.msh.geometry.x[dofs, 1] * 1e6
c_out = opt.c_h.x.array[dofs]
idx   = np.argsort(y_out)
y_s, c_s = y_out[idx], c_out[idx]
y_t   = np.linspace(0, Ly*1e6, 200)

ax2.plot(y_t, y_t/(Ly*1e6), 'r--', lw=1.2, label='Target (linear)')
ax2.plot(y_s, np.clip(c_s, 0, 1), 'b-', lw=1.4, label='TopOpt result')
ax2.set_xlabel(r'$y$ (µm)', fontsize=8)
ax2.set_ylabel(r'Concentration $c$', fontsize=8)
ax2.set_xlim(0, Ly*1e6)
ax2.set_ylim(-0.05, 1.05)
ax2.tick_params(labelsize=7)
ax2.legend(fontsize=6.5, framealpha=0.7, loc='upper left')
ax2.set_title('(b) Outlet concentration profile', fontsize=8, loc='left', pad=4)

rmse = float(np.sqrt(np.mean((np.clip(c_s, 0, 1) - y_s/(Ly*1e6))**2)))
ax2.text(0.98, 0.06, f'RMSE = {rmse:.3f}',
         transform=ax2.transAxes, fontsize=7,
         ha='right', color='#333333',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8))

# ── 5. Save ─────────────────────────────────────────────────────────────────
out_dir = os.path.dirname(os.path.abspath(__file__))
pdf_path = os.path.join(out_dir, 'topology_optimised.pdf')
png_path = os.path.join(out_dir, 'topology_optimised.png')

fig.savefig(pdf_path, dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
fig.savefig(png_path, dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {pdf_path}")
print(f"Saved: {png_path}")
print(f"RMSE = {rmse:.4e}")
