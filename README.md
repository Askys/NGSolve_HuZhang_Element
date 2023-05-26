# Hu-Zhang elements for NGSolve

An implementation of the two-dimensional Hu-Zhang 
element in [NGSolve](https://ngsolve.org/) with 
two examples for the TFSRM and QFSRM plate formulations

## Usage with NGSolve

An example file with the QFSRM formulation of the Reissner-Mindlin plate:

```python
#Import the extension  
from ngsolve.fem import CompilePythonModule
from pathlib import Path

txt = Path('mymodule.cpp').read_text()
m = CompilePythonModule(Path('mymodule.cpp'), init_function_name='mymodule')

from ngsolve.meshes import MakeStructured2DMesh
from ngsolve import *
from ngsolve.webgui import Draw

# Define the polynomial order p >= 3
order = 3

# Set material constants
t = 1e-5
ks = 5/6
E = 1
nu = 0.3
mu = E / (2*(1+nu))

# Define mesh subdivision
h = 4

# Define the compliance operator
def A(mat):
    return 12*(1+nu)/E*(mat-nu/(nu+1)*Trace(mat)*Id(2))

# Define the force
def f0(a):
    return a*(a-1)

def f1(a):
    return 5 * a ** 2 - 5 * a + 1

def f2(a):
    return 2*a - 1

g = 200 * E / (1-nu**2) * ( f0(x)**3*f1(y) + f0(y)**3*f1(x) + f0(x)*f0(y)*f1(x)*f1(y) )

# Build the mesh
mesh = MakeStructured2DMesh(quads=False, nx=h, ny=h)

# Define the finite element spaces
fesL2 = L2(mesh, order=order-1)
fesVL2 = VectorL2(mesh, order=order-1)
fesHsD = m.HDivSymFESpace(mesh, order=order)
fesHd = HDiv(mesh, order=order-1, RT=True)
fes = fesL2*fesVL2*fesHsD*fesHd

# Set the bilinear form
(w, p, M, q), (dw, dp, dM, dq) = fes.TnT()
a = BilinearForm(fes, symmetric=True, symmetric_storage=True, condense=False)
a += (InnerProduct(dM, A(M))
      + InnerProduct(div(dM), p)
      + InnerProduct(dp, div(M))
      + (t**2/(ks*mu))*InnerProduct(dq, q)
      - InnerProduct(div(dq), w)
      - InnerProduct(dw, div(q))
      - InnerProduct(dq, p)
      - InnerProduct(dp, q)
      )*dx

# Define the right-hand side
f = LinearForm(fes)
f += (-dw * g)*dx

# Solve the system
sol = GridFunction(fes)

r = sol.vec.CreateVector()
w = sol.vec.CreateVector()

with TaskManager():
    f.Assemble()
    a.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(a.condense), inverse="umfpack")

    a.Apply(sol.vec, r)

    r.data -= f.vec
    if a.condense:
        r.data += a.harmonic_extension_trans * r
    w.data = inv * r
    if a.condense:
        w.data += a.harmonic_extension * w
        w.data += a.inner_solve * r
    sol.vec.data -= w

# Get and plot results
gfw,gfp,gfM,gfq = sol.components

Draw(gfw, mesh, "w", deformation=True, order=order)
Draw(gfp, mesh, "p", order=order)
Draw(Norm(gfM), mesh, "M", order=order)
Draw(gfq, mesh, "q", order=order)
```

