# solver.py — penalty method (gamma = 1e12)
import basix.ufl
from dolfinx import fem
from dolfinx.fem.petsc import LinearProblem
import ufl
from petsc4py import PETSc
from .utilities import alpha, D_eff

def forward_solve(msh, boundary_data, rho_phys, mu=1e-3, D_fluid=1e-9,
                  P_in=1000.0, solver_type="direct"):
    i1    = boundary_data["inlet1"]
    i2    = boundary_data["inlet2"]
    out   = boundary_data["outlet"]
    walls = boundary_data["walls"]
    fd    = msh.topology.dim - 1

    P2 = basix.ufl.element("Lagrange", msh.topology.cell_name(), 2,
                            shape=(msh.geometry.dim,))
    P1 = basix.ufl.element("Lagrange", msh.topology.cell_name(), 1)
    W  = fem.functionspace(msh, basix.ufl.mixed_element([P2, P1]))
    Vc = fem.functionspace(msh, ("Lagrange", 1))

    (u, p) = ufl.TrialFunctions(W)
    (v, q) = ufl.TestFunctions(W)

    a = (mu * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
         + alpha(rho_phys) * ufl.inner(u, v) * ufl.dx
         - p * ufl.div(v) * ufl.dx
         - q * ufl.div(u) * ufl.dx)
    L = (ufl.inner(fem.Constant(msh, PETSc.ScalarType((0.0, 0.0))), v) * ufl.dx
         + fem.Constant(msh, PETSc.ScalarType(0.0)) * q * ufl.dx)

    gamma = fem.Constant(msh, PETSc.ScalarType(1e12))
    ft = boundary_data["facet_tag"]

    ds_walls = ufl.Measure("ds", domain=msh, subdomain_data=ft, subdomain_id=0)
    ds_in1   = ufl.Measure("ds", domain=msh, subdomain_data=ft, subdomain_id=1)
    ds_in2   = ufl.Measure("ds", domain=msh, subdomain_data=ft, subdomain_id=2)
    ds_out   = ufl.Measure("ds", domain=msh, subdomain_data=ft, subdomain_id=3)

    a += gamma * ufl.inner(u, v) * ds_walls

    a += gamma * p * q * ds_in1
    L += gamma * fem.Constant(msh, PETSc.ScalarType(P_in)) * q * ds_in1

    a += gamma * p * q * ds_in2
    L += gamma * fem.Constant(msh, PETSc.ScalarType(P_in)) * q * ds_in2

    a += gamma * p * q * ds_out
    # outlet pressure zero → no extra L term

    wh = LinearProblem(a, L, bcs=[],
                       petsc_options={"ksp_type": "preonly", "pc_type": "lu"}).solve()
    uh = wh.sub(0).collapse()
    ph = wh.sub(1).collapse()

    # Convection‑diffusion
    c, d = ufl.TrialFunction(Vc), ufl.TestFunction(Vc)
    D   = D_eff(rho_phys, D_fluid)
    hc  = ufl.CellDiameter(msh)
    um  = ufl.sqrt(ufl.dot(uh, uh) + 1e-20)
    Pe  = um * hc / (2.0 * D)
    tau = hc / (2.0 * um) * (1.0 / ufl.tanh(Pe) - 1.0 / Pe)

    a_c = (ufl.dot(uh, ufl.grad(c)) * d * ufl.dx
           + D * ufl.inner(ufl.grad(c), ufl.grad(d)) * ufl.dx
           + tau * ufl.dot(uh, ufl.grad(c)) * ufl.dot(uh, ufl.grad(d)) * ufl.dx)
    L_c = fem.Constant(msh, PETSc.ScalarType(0.0)) * d * ufl.dx

    bc_c1 = fem.dirichletbc(PETSc.ScalarType(1.0),
                             fem.locate_dofs_topological(Vc, fd, i1), Vc)
    bc_c2 = fem.dirichletbc(PETSc.ScalarType(0.0),
                             fem.locate_dofs_topological(Vc, fd, i2), Vc)

    ch = LinearProblem(a_c, L_c, bcs=[bc_c1, bc_c2],
                       petsc_options={"ksp_type": "preonly", "pc_type": "lu"}).solve()
    return uh, ph, ch