# brinkgrad/manufacturability.py (updated)
import numpy as np, matplotlib.pyplot as plt
try:
    from scipy.ndimage import distance_transform_edt
    from scipy.interpolate import griddata
except ImportError:
    distance_transform_edt = None
    griddata = None

def measure_min_feature_size(rho_phys, V_rho, threshold=0.5, resolution=1e-6,
                             plot_distribution=False, output_file=None):
    fab_limit_um = 10.0; fab_limit_m = fab_limit_um*1e-6
    pts = V_rho.tabulate_dof_coordinates()[:,:2]
    vals = (rho_phys.x.array > threshold).astype(np.float64)
    x_min, y_min = pts.min(axis=0); x_max, y_max = pts.max(axis=0)
    nx = int((x_max-x_min)/resolution)+1; ny = int((y_max-y_min)/resolution)+1
    x_grid = np.linspace(x_min, x_max, nx); y_grid = np.linspace(y_min, y_max, ny)
    X, Y = np.meshgrid(x_grid, y_grid)
    binary_grid = griddata(pts, vals, (X, Y), method='nearest', fill_value=0.0)
    dist_pixels = distance_transform_edt(binary_grid==1); dist_m = dist_pixels * resolution
    fluid_mask = binary_grid==1
    if not np.any(fluid_mask): return {"min_width_m":0.0, "avg_width_m":0.0, "min_width_um":0.0, "avg_width_um":0.0, "fab_ok":False, "fab_limit_um":fab_limit_um}
    min_width = 2 * dist_m[fluid_mask].min(); avg_width = 2 * np.mean(dist_m[fluid_mask])
    ok = min_width >= fab_limit_m
    if plot_distribution:
        plt.figure(); widths_um = 2 * dist_m[fluid_mask] * 1e6
        plt.hist(widths_um, bins=50, alpha=0.7, color='steelblue')
        plt.axvline(fab_limit_um, color='r', linestyle='--', label=f'Fab limit ({fab_limit_um} µm)')
        plt.axvline(min_width*1e6, color='k', linestyle=':', label=f'Min width = {min_width*1e6:.1f} µm')
        plt.xlabel('Channel width [µm]'); plt.ylabel('Count'); plt.legend(); plt.title('Channel width distribution')
        if output_file: plt.savefig(output_file, dpi=150)
        plt.close()
    return {"min_width_m":min_width, "avg_width_m":avg_width, "min_width_um":min_width*1e6, "avg_width_um":avg_width*1e6, "fab_ok":ok, "fab_limit_um":fab_limit_um}

def robust_projection(rho_filt, beta, eta_d, V_rho):
    """Three-field robust projection: eroded, nominal, dilated.
    Wang et al. (2011) doi:10.1007/s00158-010-0602-y
    eta_d: dilation threshold (e.g. 0.3); eroded uses 1-eta_d.
    """
    import ufl
    from dolfinx import fem
    try:
        _pts = V_rho.element.interpolation_points
        if callable(_pts): _pts = _pts()
    except Exception:
        _pts = V_rho.element.interpolation_points()

    def _heaviside_expr(rho, b, eta):
        return (ufl.tanh(b*eta) + ufl.tanh(b*(rho-eta))) / (
                ufl.tanh(b*eta) + ufl.tanh(b*(1.0-eta)))

    eroded  = fem.Function(V_rho)
    nominal = fem.Function(V_rho)
    dilated = fem.Function(V_rho)
    for fn, eta in [(eroded, 1.0-eta_d), (nominal, 0.5), (dilated, eta_d)]:
        expr = _heaviside_expr(rho_filt, float(beta), eta)
        fn.interpolate(fem.Expression(expr, _pts))
        fn.x.scatter_forward()
    return eroded, nominal, dilated