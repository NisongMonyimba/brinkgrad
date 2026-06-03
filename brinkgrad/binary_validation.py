"""
Binary-geometry validation on the Brinkman mesh.
No pyvista or Gmsh required.

Three forward solves after optimization:
  (0) Gray Brinkman  — the optimization result (smooth rho)
  (A) Binary Brinkman — hard-threshold rho to {0,1}, keep smooth alpha interp
  (B) Binary Stokes   — hard-threshold rho, alpha=0 in fluid (true free-flow)

Compares RMSE and hydraulic resistance to assess porous-leakage error.
"""
import numpy as np
from dolfinx import fem
from micrograd.solver import forward_solve
import micrograd.utilities as _ut


def compute_outlet_rmse(c_h, msh, outlet_facets, Ly):
    """Compute RMSE at outlet P1 nodes vs linear target."""
    V_c  = c_h.function_space
    dofs = fem.locate_dofs_topological(V_c, 1, outlet_facets)
    coords = msh.geometry.x
    dof_coords = V_c.tabulate_dof_coordinates()
    y   = dof_coords[dofs, 1]
    c   = c_h.x.array[dofs]
    idx = np.argsort(y)
    rmse = float(np.sqrt(np.mean((c[idx] - y[idx]/Ly)**2)))
    return rmse, y[idx], c[idx]


def compute_flow_rate(u_h, msh, outlet_facets):
    """Integrate u_x over outlet facet (m^2/s per unit depth)."""
    import ufl
    from dolfinx.fem import assemble_scalar, form as fem_form
    n   = ufl.FacetNormal(msh)
    ds  = ufl.Measure("ds", domain=msh, subdomain_data=None)
    # Simple: sum u_x * dof values at outlet nodes
    V_u = u_h.function_space
    dofs_x = fem.locate_dofs_topological(V_u.sub(0).collapse()[0], 1,
                                          outlet_facets)[0]
    u_x = u_h.sub(0).collapse().x.array
    # Use mean u_x * outlet width as approximation
    Ly_vals = msh.geometry.x[dofs_x // 2, 1] if len(dofs_x) > 0 else None
    return float(np.mean(np.abs(u_h.sub(0).collapse().x.array))) * \
           float(msh.geometry.x[:, 1].max())


def run_binary_validation(opt, threshold=0.5, P_in=1000.0, verbose=True):
    msh     = opt.msh
    bdata   = opt.boundary_data
    Ly      = float(msh.geometry.x[:, 1].max())
    Lx      = float(msh.geometry.x[:, 0].max())
    rho_arr = opt.rho_phys.x.array.copy()

    outlet_facets = bdata["outlet"]
    results = {}

    def rmse_and_R(u_h, c_h, label):
        rmse, y, c = compute_outlet_rmse(c_h, msh, outlet_facets, Ly)
        # Q = integral of u_x at outlet — approximate via DOF mean
        V_ux = u_h.sub(0).collapse()
        dofs_out = fem.locate_dofs_topological(
            V_ux.function_space, 1, outlet_facets)
        u_out = V_ux.x.array[dofs_out]
        # outlet width ≈ Ly, so Q ≈ mean(u_x) * Ly
        Q = float(np.mean(u_out)) * Ly
        R = P_in / max(abs(Q), 1e-30)
        if verbose:
            print(f"  {label:25s}: RMSE={rmse:.4f}  Q={Q:.3e} m²/s  "
                  f"R={R:.3e} Pa·s/m²")
        return rmse, R

    if verbose:
        print("=" * 65)
        print("BINARY VALIDATION — comparing gray vs binary designs")
        print(f"  Gray-zone fraction: "
              f"{np.mean((rho_arr>0.05)&(rho_arr<0.95))*100:.1f}%")
        print(f"  Fluid volume fraction: {np.mean(rho_arr>threshold)*100:.1f}%")
        print("=" * 65)

    # ── (0) Gray Brinkman (optimization result) ───────────────────────────────
    r, R = rmse_and_R(opt.u_h, opt.c_h, "Gray Brinkman (opt)")
    results.update(rmse_gray=r, R_gray=R)

    # ── Build binary rho field ────────────────────────────────────────────────
    rho_bin = fem.Function(opt.rho_phys.function_space)
    rho_bin.x.array[:] = np.where(rho_arr > threshold, 1.0, 0.0)
    rho_bin.x.scatter_forward()

    # ── (A) Binary Brinkman ───────────────────────────────────────────────────
    u_a, _, c_a = forward_solve(msh, bdata, rho_bin, P_in=P_in)
    r, R = rmse_and_R(u_a, c_a, "Binary Brinkman")
    results.update(rmse_binary_brinkman=r, R_binary_brinkman=R)

    # ── (B) Binary Stokes (alpha=0 in fluid) ──────────────────────────────────
    orig_min     = _ut.alpha_min
    _ut.alpha_min = 0.0          # pure Stokes in fluid regions
    u_b, _, c_b = forward_solve(msh, bdata, rho_bin, P_in=P_in)
    _ut.alpha_min = orig_min
    r, R = rmse_and_R(u_b, c_b, "Binary Stokes (alpha=0)")
    results.update(rmse_binary_stokes=r, R_binary_stokes=R)

    if verbose:
        print("=" * 65)
        g = results["rmse_gray"]
        bb = results["rmse_binary_brinkman"]
        bs = results["rmse_binary_stokes"]
        print(f"  Binary Brinkman ΔRMSE vs gray: {(bb-g)/g*100:+.1f}%")
        print(f"  Binary Stokes   ΔRMSE vs gray: {(bs-g)/g*100:+.1f}%")
        if abs((bs-g)/g) < 0.15:
            print("  ✓ Brinkman approximation validated: <15% RMSE change")
        else:
            print("  ✗ WARNING: >15% change — porous leakage may be significant")
        print("=" * 65)

    results["gray_frac"]  = float(np.mean((rho_arr>0.05)&(rho_arr<0.95)))
    results["binary_vol"] = float(np.mean(rho_arr > threshold))
    return results


if __name__ == "__main__":
    import json, pathlib
    from micrograd import GradientGeneratorOptimizer

    print("Running 400-iteration optimisation + binary validation...")
    opt = GradientGeneratorOptimizer(
        Lx=2000e-6, Ly=500e-6, nx=80, ny=20,
        target_expr=lambda x: x[1]/500e-6,
        w_f=1e-3, w_c=5e1, V_star=0.5
    )
    opt.run(max_iter=400,
            beta_continuation=[1,2,4,8,16],
            move=0.05)

    results = run_binary_validation(opt, verbose=True)
    pathlib.Path("/tmp/binary_validation.json").write_text(
        json.dumps(results, indent=2))
    print("\nSaved to /tmp/binary_validation.json")
