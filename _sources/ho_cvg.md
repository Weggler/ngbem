# Convergence Rates for HO BEM 


NG-BEM implements high order boundary element methods. The accuracy of the method depends on the order of the finite element space which needs to be discretized. On this page some numerical results of NG-BEM are presented and 

1) Convergence rates known from theory are verified. To get there we consider analytically solvable problems with an unkown in one of the trace spaces

$$ H^{\frac12}(\Gamma), \quad \boldsymbol H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma) \quad \textnormal{and} \quad H^{-\frac12}(\Gamma)\,.$$

The unkown trace is computed with a BEM for different order of approximation $p$ and for a sequence of uniform refined meshes and $L_2$-error as function of degrees of freedom is analysed. 

2) High accuracy of the pde solution given by a boundary integral representation is exemplarily shown. Here, an electrostatic potential in $H^{1}(\Omega)$, solution of Laplace bvp, is considered.


**Test 1: Laplace Dirichlet trace in $H^{\frac12}(\Gamma)$** 

Consider the harmonic function $u$,   

$$ u: \boldsymbol x \mapsto \dfrac{1}{\| \boldsymbol x - \boldsymbol x_s\|}, \quad \boldsymbol x_s = (1,1,1)^\intercal \in \Omega^c\,, $$ 

as unique solution of the Neumann boundary value problem

|  |  |  |
| -|--|- |
|$\begin{array}{rcl l} \Delta u &=& 0, \quad  & \mathrm{in} \, \Omega\,, \\ \gamma_1 u &=& u_1, \quad & \mathrm{on} \, \Gamma\,. \end{array}$    | $\quad\quad\quad$  | ![](demos/resources/BEM.png)  |

Let $\Omega$ be the unit ball and $\Gamma$ the unit sphere, respectively. Using the NG-BEM solver we compute the Dirichlet data $u_0\in H^{\frac12}(\Gamma)$ for varying order $p$. The following diagram shows the $L_2$-error of the approximated Dirichlet data as function of the number of degrees of freedom for $p=1,...,4$.

![](demos/resources/plot_Laplace_Dirichlet.png) 


**Test 2: Laplace Neumann trace in $H^{-\frac12}(\Gamma)$** 

To analyse convergence rates for traces in $H^{-\frac12}(\Gamma)$ consider now the Dirichlet boundary value problem with unkonwn Neumann data, i.e., 


|  |  |  |
| -|--|- |
|$\begin{array}{rcl l} \Delta u &=& 0, \quad  & \mathrm{in} \, \Omega\,, \\ \gamma_0 u &=& u_0, \quad & \mathrm{on} \, \Gamma\,. \end{array}$    | $\quad\quad\quad$  | ![](demos/resources/BEM.png)  |

Thus, again $u$ from Test 1 is the unique solution of this problem. Using the NG-BEM solver we compute the Neumann data $u_1\in H^{-\frac12}(\Gamma)$ for varying order $p$. The following diagram shows the $L_2$-error of the approximated Neumann data as function of the number of degrees of freedom for $p=1,...,4$.     

![](demos/resources/plot_Laplace_Neumann.png) 


**Common notes for Test 1 and 2:** 

* We actually considered the Laplace operator with mixed boundary values and solved the bvp only once.
* Our manufactured solution $u$ does not depend on the mesh, i.e., it is a solution for all meshes. This is why the geometrical approximation order is $1$ for all test runs.
* Check out the repository's test folder if you like to run this test.



**Test 3: Laplace Represenation formula for potential in $H^{1}(\Omega)$** 

Given complete Cauchy data $u_0 \in H^{\frac12}(\Gamma)$ and $u_1 \in H^{-\frac12}(\Gamma)$, the solution $u$ of the boundary value problem 

$$\begin{array}{rcl l} \Delta u &=& 0, \quad  & \mathrm{in} \, \Omega\,, \\ \gamma_0 u &=& u_0, \quad & \mathrm{on} \, \Gamma_0\,,\\ \gamma_1 u &=& u_1, \quad & \mathrm{on} \, \Gamma_1\,, \end{array}$$

is given by the represenation formula 

$$u(x) = \int\limits_{\Gamma} G(x-y)\,u_1(y) \, \mathrm{d}\sigma_y - \int\limits_{\Gamma} \langle n(y), \nabla_y G(x-y) \rangle \, u_0(y) \, \mathrm{d} \sigma_y\,.$$


In order to illustrate the high accuracy of the high order BEM, we evaluate the representation formula for numerically computed Cauchy data. In our test setting $\Omega$ is again the unit ball and $u$ is evaluated on a rectangular screen placed inside the ball (see figure below). 

The following diagram shows the $L_2$-error of solution $u$ as function of degrees of freedom for $p=1,2,3,4$. Note that $p$ is the order of approximation we use for the BEM, i.e., to approximate the traces $u_0$ on $\Gamma_1$ and $u_1$ on $\Gamma_0$, respectively.


|  |  |  |
| -|--|- |
| ![](demos/resources/plot_Laplace_Solution_on_Screen.png)    | $\quad$  | ![](demos/resources/BEM_Screen.png)  |


**Test 4: Maxwell Neumann trace in $H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma)$** 

Besides electrostatics we can also solve boundary value problems from electromagnetics. The relevant trace spaces are $\boldsymbol H^{-\frac12}(\mathrm{curl}_\Gamma, \Gamma)$, and its dual, $\boldsymbol H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma)$. NG-BEM offers finite element spaces for both trace spaces. 

In order to verify convergence rates for high order bem, we consider the scattering of a plane wave $\boldsymbol E_{\mathrm{inc}}$ at a perfect electric conducting sphere. This scattering problem is analytically solvable by so-called Mie series. In order to keep the geometrical approximation error negligible, the geometrical approximation oder is set to $4$ for all tests runs. 

The scattered electric field $\boldsymbol E$ solves the following boundary value problem: 

|  |  |  |
| -|--|- |
|$\begin{array}{rcl l} \mathbf{\mathrm{curl}}\,\mathbf{\mathrm{curl}}\, \boldsymbol E - \kappa^2 \, \boldsymbol E&=& \boldsymbol 0, \quad  & \mathrm{in} \, \Omega^c\,, \\ \gamma_R \boldsymbol E &=& \boldsymbol E_{\mathrm{inc}}, \quad & \mathrm{on} \, \Gamma\,, \\ \textnormal{ + radiation  condition} & & &\|x\| \to \infty\,. \end{array}$    | $\quad$  | ![](demos/resources/BEM.png)  |

Using the NG-BEM solver we compute the Neumann data $\gamma_N \boldsymbol E \in \boldsymbol H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma)$ for varying orders of approximation $p$. The following diagram shows the $L_2$-error of the approximated Neumann data as function of the number of degrees of freedom for $p=0,...,4$::w
     

![](demos/resources/plot_Maxwell_Mie.png) 


**Notes:** 

* The Mie series solves the scattering at the sphere. A high order mesh is required to obtain high order convergence at least if the traces are approximated with $p>2$.   
* Check out the repository's test folder if you like to run this test.

