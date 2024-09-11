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

```{figure} ../tropical-dynamics-figures/reduced_gravity2.PNG
---
name: FIG4-1
---
As in {numref}`FIG2-2` except that $h_B\neq0$
```


Following {eq}`reduced_gravity2` (reduced gravity) and integrating the third equation from surface (or top of the boundary) to the interface $\eta$ (i.e., {numref}`FIG4-1`), we can rewrite {eq}`shallow_water_1` 

```{math}
:label: shallow_water_2
\begin{align}
\frac{\partial u}{\partial t} + [u\frac{\partial u}{\partial x}+ v\frac{\partial u}{\partial y}]-\beta yv & = -g'\frac{\partial \eta}{\partial x}\\
\frac{\partial v}{\partial t} + [u\frac{\partial v}{\partial x}+ v\frac{\partial v}{\partial y}]+\beta yu & = -g'\frac{\partial \eta}{\partial y}\\
w(z=\eta) & = \frac{d \eta}{dt} = -\overline{\eta} (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) + w_0  
\end{align}
```

where $\eta$ is the location of interface (see {eq}`reduced_gravity2`) and $w(t)$ is the rate change in this displacement i.e., $w(z)= \frac{d \eta}{dt}$. 

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
w_0=\underbrace{u\frac{\partial h_B}{\partial x}+v\frac{\partial h_B}{\partial y}}_{\text{momentum advection}}+\underbrace{h_B(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y})}_{\text{mass convergence from lower bound}}
```
{eq}`boundary_layer` states that the normal velocity exists through two processes (1) the advection of $h_B$ (inhomogeneous in space) by $u$ and $v$ at $h_B$ or (2) the convergence/divergence of $h_B$ (homogeneous but non-zero).   
:::

### Conservation of Shallow Water PV
{eq}`shallow_water_2` also conserves the potential vorticity. To prove that, we can take $\nabla\times$ of the first equation, we have 


```{math}
:label: shallow_water_4
\begin{align}
\frac{d\zeta}{dt} = \frac{\partial \zeta}{\partial t} + u\frac{\partial \zeta}{\partial x}+v\frac{\partial \zeta}{\partial y} &= -(\zeta+f)(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) \\
& = \frac{\zeta+f}{\overline{\eta}}\frac{d\eta}{dt} \text{  or}\\
& \rightarrow  \frac{1}{\zeta+f}\frac{d\zeta+f}{dt} = \frac{1}{\overline{\eta}}\frac{d\eta}{dt}
\end{align}
```

the last equation of {eq}`shallow_water_4` the rate change in absolute vorticity is proportional to the rate change in height indicating that $\frac{d}{dt}\frac{\zeta+f}{\eta} = 0$. (One should notice that the expression in the last equation of {eq}`shallow_water_2` will yield a more rigorous constraint on PV conservation). One can also prove {eq}`shallow_water_4` from {eq}`shallow_water_2` using the Kelvin circulation theorem. (HW)

### Equatorial Wave Solution 
while {eq}`shallow_water_2` is well simplified, it remains unsolvable due to its nonlinearity. We can start with a linear version or resting basic state assumption ($u=\overline{u}+u'$ where $\overline{u}$ is independent of time or $\overline{u}=0$ for resting basic state) with resting boundary condition $w_0=0$. The equation set is then simplified as 

```{math}
:label: shallow_water_linear
\begin{align}
\frac{\partial u}{\partial t} -\beta yv & = -\frac{\partial \phi'}{\partial x}\\
\frac{\partial v}{\partial t} +\beta yu & = -\frac{\partial \phi'}{\partial y}\\
\frac{\partial \phi'}{\partial t} + g'\overline{\eta} (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) &= 0  
\end{align}
```

where $g'\eta=C_0^2$, $C_0$ is gravity wave speed. The above equation is the famous equation set derived by Matsuno (1966), which is a solvable equation set. 

:::{note}
While the linearization of {eq}`shallow_water_2` significantly simplifies the problem, the trade-off also leads to a few issues. The linear momentum equation doesn't allow the existence of energy cascade into new scales (i.e., filament structure of vorticity won't happen). Such a feature is very important in redistributing the lower troposphere moisture.  
:::


To further simplify the problem, we can non-dimentionalize each variable by representing them with 1 or multiple non-dimentional parameters. 


```{math}
:label: shallow_water_linear
\begin{align}
(x,y) & = \sqrt {\frac{C_0}{\beta}}(x^{\text{non}},y^{\text{non}}) \\
t     & = \frac{t^{\text{non}}}{\sqrt{\beta C_0}} \\
u     & = C_0 u^{\text{non}} \\
\phi^{'}  & = C_0^2 \phi^{^{\text{non}}} 
\end{align}
```
