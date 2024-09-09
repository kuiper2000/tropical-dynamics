(waves1)=
# Week 6 and 7: Tropical Wave Theory
## Shallow Water Model

In addition to zonal mean circulation, eddy also plays an important role in shaping tropical climatology. Here we adopt the result in {ref}`scale_analysis`. 


```{math}
:label: shallow_water_1
\begin{align}
\frac{\partial u}{\partial t} + [u\frac{\partial u}{\partial x}+ v\frac{\partial u}{\partial y}]-\beta yv & = -\frac{1}{\rho}\frac{\partial p}{\partial x}\\
\frac{\partial v}{\partial t} + [u\frac{\partial v}{\partial x}+ v\frac{\partial v}{\partial y}]+\beta yu & = -\frac{1}{\rho}\frac{\partial p}{\partial y}\\
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z} &= 0 \\
\frac{\partial p}{\partial z} &= -\rho g \\
\end{align}
```

Following {eq}`reduced_gravity2` (reduced gravity) and integrating the third equation from top of the boundary layer to the interface $\eta$, we can rewrite {eq}`shallow_water_1` 

```{math}
:label: shallow_water_2
\begin{align}
\frac{\partial u}{\partial t} + [u\frac{\partial u}{\partial x}+ v\frac{\partial u}{\partial y}]-\beta yv & = -g'\frac{\partial h}{\partial x}\\
\frac{\partial v}{\partial t} + [u\frac{\partial v}{\partial x}+ v\frac{\partial v}{\partial y}]+\beta yu & = -g'\frac{\partial h}{\partial y}\\
w(z=\eta) & = \frac{d \eta}{dt} = -\overline{\eta} (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) + w_0  
\end{align}
```

where $h$ is the location of interface defined as $\eta(x)-\overline{\eta}$ (see {eq}`reduced_gravity2`) and $w(t)$ is the rate change in this displacement i.e., $w(z)=\frac{d h}{dt} = \frac{d \eta}{dt}$. 

Also, according to quasi-equilibrium, we know diabatic heating/cooling is balanced by adiabatic cooling/heating of vertical motion. i.e., {eq}`thermodynamics`. Therefore, the last equation is also subject to a forcing term, which can be written as 

```{math}
:label: shallow_water_3
\begin{align}
w(z) = \frac{d \eta}{dt} = -\overline{\eta} (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) + w_0  & = \frac{g}{N^2}\frac{Q}{c_p T} \\
& = \frac{g}{g\frac{\partial \mathrm{ln}\theta}{\partial z}}\frac{Q}{c_p T} \\
& = \frac{1}{\frac{\partial \mathrm{ln}T}{\partial z}-\frac{R}{c_p}\frac{\partial \mathrm{ln}p}{\partial z}} \frac{Q}{c_p T} \\
& = \frac{1}{\frac{\partial T}{\partial z}+\frac{g}{c_p}}\frac{Q}{c_p} \\
& = \frac{Q}{c_p (\Gamma_d-\Gamma)} \\
\text{using geopotential height to rewrite above equation} & \\
\frac{d \phi}{dt} & =  -g\overline{\eta} (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) = \frac{Q}{ (1-\Gamma/\Gamma_d)}
\end{align}
```

where $\phi$ is $g\eta$. The last equation of {eq}`shallow_water_3` states that the diabatic heating-induced vertical motion differentiation nearly balances the convergence/divergence of horizontal flow. This relation holds for both external (i.e., g is not reduced) and internal modes (i.e., g is reduced). 

:::{note}
One should notice that $ w_0$ in {eq}`shallow_water_2` is not necessarily 0. To ensure non-normal flow at the lower boundary of the domain, it requires that $w(x,y,z=h_B,t) = u\frac{\partial h_B}{\partial x}+v\frac{\partial h_B}{\partial y}$ (at the boundary, $\frac{\partial h_B}{\partial t}=0$). This leads to 
```{math}
:label: boundary_layer
w_0=\underbrace{u\frac{\partial h_B}{\partial x}+v\frac{\partial h_B}{\partial y}}_{\text{momentum advection}}+\underbrace{h_B(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y})}_{mass convergence from lower bound}
```
{eq}`boundary_layer` states that the normal velocity exists through two processes (1) the advection of $h_B$ (inhomogeneous in space) by $u$ and $v$ at $h_B$ or (2) the convergence/divergence of $h_B$ (homogeneous but non-zero)  



:::
