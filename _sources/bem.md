Boundary Element Methods
------------------------

**What is a standard BEM problem?**

One of the following pdes with given boundary condition:

- Laplace equation $\quad \Delta u = 0$
- Helmholtz equation $ \quad -\Delta u - \kappa^2\, u = 0$
- Maxwells equations  $\quad \mathbf{curl}\mathbf{curl}\, \boldsymbol u - \kappa^2 \,\boldsymbol u = \boldsymbol 0$


Consider for instance interior Laplace with Dirichlet condition:

$$ \left\{ \begin{array}{rcl l} \Delta u &=& 0, \quad &\Omega \subset \mathbb R^3\,,\\ \gamma_0 u&=& u_0, \quad &\Gamma = \partial \Omega\,.\end{array} \right. $$ 

**How to derive a BEM from a standard BEM problem?**

1. choose a representation in terms of layer potentials:

$$ u(x) = \underbrace{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{1}{\| x-y\|} } \, j(y)\, \mathrm{d}\sigma_y }_{\displaystyle{ \mathrm{SL}(j) } }$$


2. derive a boundary integral equation for unknown trace and consider a variational formulation:

$$ \forall \, v\in H^{-\frac12}(\Gamma): \quad \left\langle \gamma_0 \left(\mathrm{SL}(j)\right), v \right\rangle_{-\frac12} = \left\langle u_0, v\right\rangle_{-\frac12} \,. $$

3. discretize the variational formulation with conforming finite element spaces:

$$ \mathrm{V} \, \mathrm{j} =   \underbrace{ \mathrm{M} \,\mathrm{u}_0}_{\displaystyle{\mathrm{rhs} }}  $$ 

4. evaluate the solution with formula from 1. for any $ x \in\Omega$!

**Why is BEM beneficial?** 

- problem dimension is reduced
- exterior problems are not an issue
- the solution is very accurate 

**Why is BEM not everybodies darling?** 

- only linear, isotropic material 
- source terms cause Newton potentials 
- singular and hypersingular integration is difficult
- system matrices are dense

#### NG-BEM Matrix Approximation

**What is ACA approximation?** 
- short explanation

**What parameters define the ACA approximation?** 
- leafsize
- eta
- eps

**Where to set the ACA parameters?** 
- when calling the layer potential operator
 
#### NG-BEM Vision 

- compatibility to NGSolve 
- efficient and accurate BEM solver  
- small, user-friendly Python interface  
- extensible to individual needsneeds   

```
j = GridFunction(fesL2)
pre = BilinearForm(u*v*ds, diagonal=True).Assemble().mat.Inverse()
with TaskManager(pajetrace=1000*1000*1000):
    V = SingleLayerPotentialOperator(fesL2, intorder=10, leafsize=40, eta=3., eps=1e-4, method="aca")  
    CG(mat = V.mat, pre=pre, rhs = rhs.vec, sol=j.vec, tol=1e-8, maxsteps=200, initialize=False)

```

#### NG-BEM Next Steps

- general EM scattering (open PEC screens and dielectrics)
- address low-frequency EM scattering
- non-trivial FEM-BEM coupling with artificial transmission boundary
- set up demos for Stokes and Lamé equations
- implement multipole approximation (gold standard)
- fast point-wise evaluation of representation formula 
- thorough testing
