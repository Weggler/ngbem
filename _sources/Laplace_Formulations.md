Boundary Integral Equations for Laplace
=============================

#### Notations of trace operators:

$$ \begin{array}{r rcl } \textnormal{Dirichlet trace} \quad & \gamma_0 u &=& u  \\[1ex] \textnormal{Neumann trace} \quad & \gamma_1 u &=& \langle \boldsymbol n,   \nabla\, u \rangle \,. \end{array} $$

#### Properties of trace spaces:

$$ \begin{array}{r rcl l} \textnormal{Dirichlet trace} \quad & \gamma_0 u &\in& H^{\frac12}\left( \Gamma\right) \quad &\textnormal{weakly continuous}\\  \textnormal{Neumann trace} \quad & \gamma_1 u &\in& H^{-\frac12}\left( \Gamma\right) \quad & \textnormal{not continuous}\,. \end{array} $$


## Interior Laplace with Dirichlet Condition

Let $u$ denote the electrostatic potential that arises under given Dirichlet boundary condition inside a source-free domain $\Omega \in \mathbb R^3$. Thus, $u$ solves the interior boundary value problem 

|  |  |  |
| -|--|- |
|$$ \left\{ \begin{array}{rcl l} \Delta u &=& 0\,, \quad & \Omega \subset \mathbb R^3\,, \\ \gamma_0 u &=& u_0\,, \quad & \Gamma = \partial \Omega\,. \end{array} \right. $$ | $\quad\quad\quad$  | ![](demos/resources/BEM_interior.png)  |
 

From here we can choose an direct or an indirect ansatz.  

$$ \begin{array}{ll rcl } 1. & \underline{\textnormal{Direct ansatz (Dirichlet to Neumann map)}} &&& \\ 
&&&&\\ 
&\textnormal{ansatz} & u(x) &=& \underbrace{ \int_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{1}{\| x-y\|} } \, u_1(y)\, \mathrm{d}\sigma_y}_{\displaystyle{ \mathrm{ SL}(u_1) }} - \underbrace{ \int_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{\langle n(y) , x-y\rangle }{\| x-y\|^3} } \, u_0(y)\, \mathrm{d}\sigma_y}_{\displaystyle{ \mathrm{DL}(u_0) }} \\ 
&&&&\\   
&\textnormal{variational formulation }  &  \left\langle \gamma_0 \left(\mathrm{SL}(u_1)\right), v \right\rangle_{-\frac12} &=& \left\langle u_0, v\right\rangle_{-\frac12} + \left\langle \gamma_0 \left(\mathrm{DL}(u_0)\right), v\right\rangle_{-\frac12}, \quad v\in H^{-\frac12}(\Gamma) \\ 
&&&& \\ & \textnormal{discretisation (MoM)} & \mathrm{u}_1 &=& \left( \frac12 \,\mathrm{M} + \mathrm{K} \right) \, \mathrm{u}_0 \\ 
&&&& \\ 
2. & \underline{\textnormal{Indirect ansatz}} &&& \\ 
&&&& \\
& \textnormal{ansatz} & u(x) &=&  \underbrace{ \int_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{1}{\| x-y\|} } \, j(y)\, \mathrm{d}\sigma_y }_{\displaystyle{ \mathrm{SL}(j) } } \\ 
&&&&\\
&\textnormal{variational formulation }  & \left\langle \gamma_0 \left(\mathrm{SL}(j)\right), v \right\rangle_{-\frac12} &=& \left\langle u_0, v\right\rangle_{-\frac12}, \quad v\in H^{-\frac12}(\Gamma) \\ 
&&&&\\ 
\\ & \textnormal{discretisation (MoM)} & \mathrm{V} \, \mathbf{j} &=& \mathrm{M} \,\mathbf{u}_0  \end{array} $$ 


## Electrostatic potential with given Neumann data

Let $u$ denote the electrostatic potential that arises under given Neumann boundary condition inside a source-free domain $\Omega \in \mathbb R^3$. Thus, $u$ solves the boundary value problem (note the marginal difference in the given trace compared to the above stated problem)

|  |  |  |
| -|--|- |
|$$ \left\{ \begin{array}{rcl l} \Delta u &=& 0\,, \quad & \Omega \subset \mathbb R^3\,, \\ \gamma_1 u &=& u_1\,, \quad & \Gamma = \partial \Omega\,. \end{array} \right. $$ | $\quad\quad\quad$  | ![](demos/resources/BEM_interior.png)  |


From here we can choose an direct or an indirect ansatz. 

$$ \begin{array}{ll rcl } 1. & \underline{\textnormal{Direct ansatz (Neumann to Dirichlet map)}} &&& \\ 
&&&&\\ 
&\textnormal{ansatz} & u(x) &=& \underbrace{ \int_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{1}{\| x-y\|} } \, u_1(y)\, \mathrm{d}\sigma_y}_{\displaystyle{ \mathrm{SL}(u_1) }} + \underbrace{ \int_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{\langle n(y), x-y \rangle}{\| x-y\|^3} } \, u_1( y)\, \mathrm{d}\sigma_y}_{\displaystyle{ \mathrm{DL}(u_1) }} \\ 
&&&&\\   
&\textnormal{variational formulation }  &  \left\langle v, \gamma_1 \left(\mathrm{DL}(u_0)\right) \right\rangle_{-\frac12}  &=& \left\langle u_1, v\right\rangle_{-\frac12} - \left\langle v, \gamma_1 \left(\mathrm{SL}(u_1)\right) \right\rangle_{-\frac12}, \quad v\in H^{\frac12}(\Gamma) \\ 
&&&& \\ & \textnormal{discretisation with stabilization (MoM)} & \left( \mathrm{D} + \mathrm{S}\right) \mathrm{u}_0 &=& \left( \frac12 \mathrm{M} - \mathrm{K}^\intercal\right) \, \mathrm{u}_1 \\ 
&&&& \\ 
2. & \underline{\textnormal{Indirect ansatz}} &&& \\ 
&&&& \\
& \textnormal{ansatz} & u(x) &=&  \underbrace{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{ \langle n(y), x-y \rangle }{\| x-y\|^3} } \, m(y)\, \mathrm{d}\sigma_y}_{\displaystyle{ \mathrm{DL}(m) }} \\ 
&&&&\\
&\textnormal{variational formulation }  & \left\langle v, \gamma_1 \left(\mathrm{DL}(m)\right) \right\rangle_{-\frac12} &=& \left\langle u_1, v\right\rangle_{-\frac12}, \quad v\in H^{\frac12}(\Gamma) \\ 
&&&&\\ 
\\ & \textnormal{discretisation with stabilisation(MoM)} & \left( \mathrm{D} + S\right) \, \mathrm{m} &=&  \mathrm{M}\,\mathrm{u}_1  \end{array} $$ 

## NG-BEM Python Functions 

| Symbol | Operator | trial space | test space | NG-BEM | trial NG-Solve | test NG-Solve |   
|-|-|-|-|-|-|-|
| $\mathrm V $ | single layer | $H^{-\frac12}(\Gamma)$ | $H^{-\frac12}(\Gamma)$         | SingleLayerPotentialOperator | SurfaceL2 | SurfaceL2|
| $\mathrm K $ | double layer | $H^{\frac12}(\Gamma)$ | $H^{-\frac12}(\Gamma)$          | DoubleLayerPotentialOperator | H1 | SurfaceL2 |
| $\mathrm D$ | hypersingular  | $H^{\frac12}(\Gamma)$ | $H^{\frac12}(\Gamma)$          | HypersingularOperator | H1 | H1 |
| $\mathrm K' $ | adjoint double layer | $H^{-\frac12}(\Gamma)$ | $H^{\frac12}(\Gamma)$ | DoubleLayerPotentialOperator | SurfaceL2 | H1 |               










