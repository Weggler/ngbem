Boundary Integral Equations for Maxwell 
=============================

#### Notations for relevant trace operators:

$$ \begin{array}{r rcl } \textnormal{Dirichlet trace} \quad & \gamma_R \boldsymbol u &=& \left( \boldsymbol n \times \boldsymbol u \right) \times \boldsymbol n \\[1ex] \textnormal{rotated Dirichlet trace} \quad & \gamma_D \boldsymbol u &=& \boldsymbol n \times \boldsymbol u \\ \textnormal{Neumann trace} \quad & \gamma_N \boldsymbol u &=& \dfrac{1}{\kappa} \boldsymbol n \times \mathbf{curl}\, \boldsymbol u\,,\quad \kappa = \omega\, \sqrt{\varepsilon\,\mu}\,. \end{array} $$

#### Notations for relevant trace spaces:

$$ \begin{array}{r rcl l} \textnormal{Dirichlet trace} \quad & \gamma_R \boldsymbol u &\in& H^{-\frac12}\left( \mathrm{curl}_\Gamma,\Gamma\right) \quad &\textnormal{tangential edge projection weakly continuous}\\ \textnormal{rotated Dirichlet trace} \quad & \gamma_D \boldsymbol u &\in& H^{-\frac12}\left( \mathrm{div}_\Gamma, \Gamma\right) \quad & \textnormal{normal edge projection weakly condinuous} \\ \textnormal{Neumann trace} \quad & \gamma_N \boldsymbol u &\in& H^{-\frac12}\left( \mathrm{div}_\Gamma, \Gamma\right) \quad & \textnormal{normal edge projection weakly condinuous}\,. \end{array} $$


## Maxwell equations PEC

Let $\Omega \in \mathbb R^3$ denote a perfect electric conductor and $\gamma_R \boldsymbol E^i = -\boldsymbol m$ the given Dirichlet trace of an incoming time-harmonic electromagnetic signal $\boldsymbol E^i$. The scattered electromagnetic field with components $\boldsymbol E$ and $\boldsymbol H$ solves the following equations in the exterior domain $\Omega^c$:

$$ \left\{ \begin{array}{ccl} \mathbf{curl} \, \boldsymbol H &=& -i\omega\varepsilon \boldsymbol E\,, \\ \mathbf{curl} \, \boldsymbol E &=& i\omega\mu \boldsymbol H\,, \\ \gamma_R\, \boldsymbol E &=& \boldsymbol m \end{array} \right. $$ 

From here we can derive two second order equations: one for the electric field $\boldsymbol E$ and one for the magnetic field $\boldsymbol H$.  


## Electric field boundary integral equations

The electric field $\boldsymbol E$ solves the second order equation with Dirichlet boundary: 

$$\left\{ \begin{array}{rcl l} \mathbf{curl} \, \mathbf{curl}\, \boldsymbol E - \kappa^2 \, \boldsymbol E &=& \boldsymbol 0, \quad &\textnormal{in } \Omega^c \subset \mathbb R^3\,,\\ \gamma_R \,\boldsymbol E &=& \boldsymbol m, \quad & \textnormal{on }\Gamma \\ \left\| \mathbf{curl} \, \boldsymbol E( x) - i\,\omega\,\epsilon \, \boldsymbol E( x)\right\| &=& \mathcal O\left( \displaystyle \frac{1}{\| x\|^2}\right), &\textnormal{for} \; \|x\| \to \infty\,.\end{array} \right. $$ 

From here we can choose an direct or an indirect ansatz. We look first at the direct and then on the indirect ansatz. 

### Direct ansatz 

Representation formula:

$$ \begin{array}{rcl} \boldsymbol E(x) &=& \mathrm{SL}\left( \gamma_N \, \boldsymbol E\right)(x) + \mathrm{DL}\left( \gamma_D\,\boldsymbol E\right)(x) \\[1ex] &=&\underbrace{ \kappa\,\displaystyle \int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \boldsymbol j(y)\, \mathrm{d}\sigma_y + \frac{1}{\kappa} \nabla \displaystyle\int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi}\, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \mathrm{div}_\Gamma \boldsymbol j(y)\, \mathrm{d}\sigma_y }_{\displaystyle{ =\mathrm{SL}(\boldsymbol j) } } + \underbrace{ \nabla \times \displaystyle \int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \boldsymbol n(y) \times \boldsymbol{m}(y)\, \mathrm{d}\sigma_y }_{\displaystyle{ =\mathrm{DL}(\boldsymbol n\times \boldsymbol m)} } \,. \end{array}$$

Variational formulation and discretisation:

$$ \forall \, \boldsymbol v\in H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma): \quad \left\langle \gamma_R \, \mathrm{SL} (\boldsymbol j), \boldsymbol v \right\rangle_{-\frac12} = \left\langle \boldsymbol m, \boldsymbol v\right\rangle_{-\frac12}  - \left\langle \gamma_R\,\mathrm{DL}(\boldsymbol n \times \boldsymbol{m}), \boldsymbol v\right\rangle_{-\frac12} \quad \stackrel{\textnormal{MoM}}{\Longrightarrow} \quad \mathrm{V} \, \mathbf{j} = \left( \frac12 \mathrm{M} - \mathrm{K}\right) \,\mathbf{m} \,, $$ 

### Indirect ansatz

Representation formula:

$$ \begin{array}{rcl} \boldsymbol E(x) &=& \mathrm{SL}\left( \gamma_N \, \boldsymbol E^t\right)(x)  \\[1ex] &=&\underbrace{ \kappa\,\displaystyle \int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \boldsymbol j^t(y)\, \mathrm{d}\sigma_y + \frac{1}{\kappa} \nabla \displaystyle\int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi}\, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \mathrm{div}_\Gamma \boldsymbol j^t(y)\, \mathrm{d}\sigma_y }_{\displaystyle{ =\mathrm{SL}(\boldsymbol j^t) } } \,. \end{array}$$

Variational formulation and discretisation:

$$ \forall \, \boldsymbol v\in H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma): \quad \left\langle \gamma_R \, \mathrm{SL} (\boldsymbol j^t), \boldsymbol v \right\rangle_{-\frac12} = \left\langle \boldsymbol m, \boldsymbol v\right\rangle_{-\frac12}  \quad \stackrel{\textnormal{MoM}}{\Longrightarrow} \quad \mathrm{V} \, \mathbf{j^t} = \mathrm{M} \,\mathbf{m} \,, $$ 


#### Notes

+ The bie from an indirect ansatz is often called **EFIE** (electric field integral equation).  
+ The solution $\boldsymbol j^t$ of the EFIE is the Neumann trace of the total electric field $\boldsymbol E^t = \boldsymbol E + \boldsymbol E^i$. 
+ The density $\boldsymbol j$ from the direct ansatz is not the same as the density $\boldsymbol j^t$ from the indirect ansatz. It holds  
  
$$ \begin{array}{rcl cl l} \boldsymbol j &=& \gamma_N \, \boldsymbol E &=& \dfrac{1}{\kappa} \, \boldsymbol n \times \mathbf{curl}\,\boldsymbol E \quad &\textnormal{Neumann trace of scattered field}  \\ \boldsymbol j^t &=& \gamma_N \, \boldsymbol E^t &=& \dfrac{1}{\kappa} \, \boldsymbol n \times \mathbf{curl}\,\boldsymbol E^t \quad & \textnormal{Neumann trace of total field}\,.  \end{array} $$

+ The densities are related as follows 

$$ \boldsymbol j^t = \boldsymbol j + \boldsymbol j^i$$


## Magnetic field boundary integral equations

The magnetic field solves a second order equation with Neumann boundary conditions:

$$ \left\{ \begin{array}{rcl l} \mathbf{curl} \, \mathbf{curl}\, \boldsymbol H - \kappa^2 \, \boldsymbol H &=& \boldsymbol 0, \quad &\textnormal{in } \Omega^c \subset \mathbb R^3\,,\\ \gamma_N \,\boldsymbol H &=& -\dfrac{i\omega\varepsilon}{\kappa} \, \boldsymbol n\times \boldsymbol m, \quad & \textnormal{on }\Gamma \\[1ex] \left\| \mathbf{curl} \, \boldsymbol H( x) + i\,\omega\,\mu \, \boldsymbol H( x)\right\| &=& \mathcal O\left( \displaystyle \frac{1}{\| x\|^2}\right), &\textnormal{for} \; \|x\| \to \infty\,.\end{array} \right. $$ 


Again, we can choose an direct or an indirect ansatz. We look first at the direct and then on the indirect ansatz. 

### Direct ansatz

Representation formula:

$$\begin{array}{rcl} \boldsymbol H(x) &=& \mathrm{SL}\left( \gamma_N \, \boldsymbol H\right)(x) +\mathrm{DL}\left( \gamma_R\,\boldsymbol H\right)(x) \\[1ex] &=&  -\dfrac{i\omega\varepsilon}{\kappa} \, \Big(\underbrace{ \kappa\, \displaystyle\int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \boldsymbol n(y)\times \boldsymbol m(y)\, \mathrm{d}\sigma_y + \frac{1}{\kappa} \nabla \displaystyle\int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi}\, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \mathrm{curl}_\Gamma \,\boldsymbol m(y)\, \mathrm{d}\sigma_y }_{\displaystyle{\mathrm{SL}(\boldsymbol n \times \boldsymbol m)} } \Big)  + \dfrac{ \kappa }{ i\omega\mu} & \underbrace{ \nabla \times \displaystyle\int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \boldsymbol{j}(y) \, \mathrm{d}\sigma_y }_{\displaystyle{ \mathrm{DL} (\boldsymbol{j}) } }\,.\end{array}$$

Apply rotated Dirichlet trace: 

$$\begin{array}{c rcl} & \gamma_D \,\boldsymbol H &=& -\dfrac{i\omega\varepsilon}{\kappa} \gamma_D \,\mathrm{SL}(\boldsymbol n\times \boldsymbol m) + \dfrac{\kappa}{i\omega\mu} \gamma_D \,\mathrm{DL}( \boldsymbol j) \\[1ex] \Rightarrow & \dfrac{\kappa}{i\omega\mu}\boldsymbol j &=& -\dfrac{i\omega\varepsilon}{\kappa} \gamma_D \, \mathrm{SL}(\boldsymbol n\times \boldsymbol m) + \dfrac{\kappa}{i\omega\mu} \gamma_D \,\mathrm{DL}( \boldsymbol j) \quad \Rightarrow \quad \boldsymbol j = \gamma_D \,\mathrm{SL}( \boldsymbol n\times \boldsymbol m) + \gamma_D\, \mathrm{DL}(\boldsymbol j)  \end{array}$$

Variational formulation and discretisation:

$$ \forall \, \boldsymbol v\in H^{-\frac12}(\mathrm{curl}_\Gamma, \Gamma): \quad \left\langle \boldsymbol v, \boldsymbol j\right\rangle_{-\frac12} - \left\langle \boldsymbol v, \gamma_D \,\mathrm{DL} (\boldsymbol j) \right\rangle_{-\frac12} = \left\langle \boldsymbol v,  \gamma_D\, \mathrm{SL}(\boldsymbol n \times \boldsymbol{m}) \right\rangle_{-\frac12} \quad \stackrel{\textnormal{MoM}}{\Longrightarrow} \quad  \left( \frac12 \mathrm{M}^\intercal + \mathrm{K}^\intercal\right) \,\mathbf{j} = -\mathrm D \, \mathbf m \,, $$ 

### Indirect ansatz

Representation formula:

$$ \boldsymbol H(x) = \mathrm{DL}\left( \gamma_R\,\boldsymbol H^t\right)(x) = \dfrac{ \kappa }{ i\omega\mu} \underbrace{ \nabla \times \displaystyle\int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \boldsymbol{j}^t(y) \, \mathrm{d}\sigma_y }_{\displaystyle{ \mathrm{DL} (\boldsymbol{j}^t) } }\,.$$

Apply rotated Dirichlet trace and use $\boldsymbol j = \boldsymbol j^t - \boldsymbol j^i$: 

$$\begin{array}{ l c rcl}  \gamma_D \,\boldsymbol H =  \dfrac{\kappa}{i\omega\mu} \gamma_D \,\mathrm{DL}( \boldsymbol j^t) \quad &\Rightarrow & \dfrac{\kappa}{i\omega\mu}\boldsymbol j &=&  \dfrac{\kappa}{i\omega\mu} \gamma_D \,\mathrm{DL}( \boldsymbol j^t) \\[2ex] &\Rightarrow &\boldsymbol j^t &=& \gamma_D \,\mathrm{DL}(\boldsymbol j^t) + \boldsymbol j^i  \end{array}$$


Variational formulation and discretisation:

$$ \forall \, \boldsymbol v\in H^{-\frac12}(\mathrm{curl}_\Gamma, \Gamma): \quad \left\langle \boldsymbol v, \boldsymbol j^t\right\rangle_{-\frac12} - \left\langle \boldsymbol v, \gamma_D \,\mathrm{DL} (\boldsymbol j^t) \right\rangle_{-\frac12} = \left\langle \boldsymbol v,  \boldsymbol{j}^i \right\rangle_{-\frac12} \quad \stackrel{\textnormal{MoM}}{\Longrightarrow} \quad  \left( \frac12 \mathrm{M}^\intercal + \mathrm{K}^\intercal\right) \,\mathbf{j} = \mathrm M \, \mathbf j^i \,, $$ 

### Notes

+ The BIE from an indirect ansatz is often called **MFIE** (magnetic field integral equation).  
- Also the magnetic field leads to boundary integral equations  $\boldsymbol j$ and  $\boldsymbol j^t$. The boundary integral equation is are integral equations of second type. 
- consider $H^{-\frac12}(\mathrm{curl}_\Gamma,\Gamma)$ conforming finite elements for test and trial space. Here is the hypersingular operator entry $lk$

$$ \mathrm{D}_{lk} = \int\limits_\Gamma \int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \langle \boldsymbol n(y)\times \boldsymbol v_l(y), \boldsymbol n(x) \times \boldsymbol \varphi_k(x)\rangle\, \mathrm{d}\sigma_y \, \mathrm{d}\sigma_x - \int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi}\, \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \, \mathrm{curl}_\Gamma \,\boldsymbol v_l(y)\, \mathrm{curl}_\Gamma\,\boldsymbol \varphi_k(x) \mathrm{d}\sigma_y \mathrm{d}\sigma_x   $$

- consider a trial function $\boldsymbol \psi_k \in H^{-\frac12}(\mathrm{div}_\Gamma,\Gamma)$ and a test function $\boldsymbol v_l \in H^{-\frac12}(\mathrm{curl}_\Gamma,\Gamma)$. Here is the adjoint double layer potential operator entry $lk$

$$ \mathrm{K}'_{lk} = \int\limits_\Gamma \int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \big\langle \boldsymbol n(y)\times \boldsymbol v_l(y), \nabla_{x} \frac{e^{i\,\kappa\,\|x-y\|}}{\| x-y\|} } \times \boldsymbol \psi_k(y) \big\rangle\, \mathrm{d}\sigma_y \, \mathrm{d}\sigma_x    $$


### Conclusion

| Operator Name | Symbol | trial space | test space | rotation | 
|-|-|-|-|-|
| single layer operator | $\mathrm V $ | $H^{-\frac12}(\mathrm{div}_\Gamma,\Gamma)$ | $H^{-\frac12}(\mathrm{div}_\Gamma,\Gamma)$ | - |
| double layer operator | $\mathrm K $ | $H^{-\frac12}(\mathrm{curl}_\Gamma,\Gamma)$ | $H^{-\frac12}(\mathrm{div}_\Gamma,\Gamma)$ | trial space |
| hypersingular operator | $\mathrm D$ | $H^{-\frac12}(\mathrm{curl}_\Gamma,\Gamma)$ | $H^{-\frac12}(\mathrm{curl}_\Gamma,\Gamma)$ | - |
| adjoint double layer operator | $\mathrm K' $ | $H^{-\frac12}(\mathrm{div}_\Gamma,\Gamma)$ | $H^{-\frac12}(\mathrm{curl}_\Gamma,\Gamma)$ | test space |
| mass matrix | $\mathrm M $ | $H^{-\frac12}(\mathrm{div}_\Gamma,\Gamma)$ | $H^{-\frac12}(\mathrm{curl}_\Gamma,\Gamma)$ | - |

### NG-BEM Python Functions 

| Operator | Python Function  | 
|-|-|
| $\mathrm V $ | MaxwellSingleLayerPotentialOperator | 
| $\mathrm K $ | MaxwellDoubleLayerPotentialOperator | 
| $\mathrm D $ | MaxwellSingleLayerPotentialOperatorCurl |






