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

One easy way to make {eq}`shallow_water_3` an internal wave equation is calculating the difference between two vertical layers (such that you have a counteracting force from the fluid on the top). or an alternative but more useful form of the last equation of {eq}`shallow_water_3` is 

```{math}
:label: hydrostatic_thermo
\frac{d \phi_p}{dt} +\sigma \omega = -\frac{Q}{ c_p T} \alpha
```

where $\phi_p =-\alpha$ , $\omega=\overline{\rho}gw$ and $\sigma = \frac{R}{p \rho g}(\Gamma_d-\Gamma) = \frac{N^2}{\rho^2 g^2}=\frac{C_0}{p^2}$. $C_0$ is the external gravity wave speed. {eq}`hydrostatic_thermo` is the hydrostatic thermodynamics equation, which can be derived by rewriting $d \mathrm{ln}\theta$ in {eq}`primitive` with $\frac{c_p}{R}\mathrm{ln}\alpha$. {eq}`hydrostatic_thermo` is an easier form to approach the solution of vertical since the pressure weighting has been considered in the pressure velocity.  

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



### Equatorial Wave Equations  
while {eq}`shallow_water_2` is well simplified, it remains unsolvable due to its nonlinearity. We can start with a linear version or resting basic state assumption ($u=\overline{u}+u'$ where $\overline{u}$ is independent of time or $\overline{u}=0$ for resting basic state) with resting boundary condition $w_0=0$. The equation set is then simplified as 

```{math}
:label: shallow_water_linear
\begin{align}
\frac{\partial u}{\partial t} -\beta yv & = -\frac{\partial \phi'}{\partial x}\\
\frac{\partial v}{\partial t} +\beta yu & = -\frac{\partial \phi'}{\partial y}\\
\frac{\partial u}{\partial x} +\frac{\partial v}{\partial y} + \frac{\partial \omega}{\partial p}& = 0 \\
\frac{\partial \phi_p'}{\partial t} + \frac{C_0^2}{p^2} \omega &= 0  
\end{align}
```

where $g'\eta=C_0^2$, $C_0$ is gravity wave speed. The above equation is the famous equation set derived by Matsuno (1966), which is a solvable equation set. 

One can find {eq}`shallow_water_linear` has four unknowns and four equations, which makes it solvable.  On the other hand, {eq}`shallow_water_linear` has a separable form. Thus one can first assume a solution of $u,v,\phi = \{\widehat{u}, \widehat{v}, \widehat{\phi}\}(x,y,t)\frac{dW(p)}{dp}$ and $w=\widehat{w}(x,y,t)W(p)$ to separate horizontal and vertical structure. Substitute back into {eq}`shallow_water_linear`, we have 

```{math}
:label: shallow_water_linear2
\begin{align}
\frac{\partial \widehat{u}}{\partial t} -\beta yv & = -\frac{\partial \widehat{\phi'}}{\partial x}\\
\frac{\partial \widehat{v}}{\partial t} +\beta yu & = -\frac{\partial \widehat{\phi'}}{\partial y}\\
\frac{\partial \widehat{u}}{\partial x} +\frac{\partial \widehat{v}}{\partial y} + \widehat{\omega}& = 0 \\
\frac{\partial \widehat{\phi'}}{\partial t}\frac{d^2W}{dp^2} + \frac{C_0^2}{p^2} \widehat{\omega}W(p) &= 0  
\end{align}
```

If we rearrange the last equation of {eq}`shallow_water_linear2` as $\frac{d^2 W}{dp^2}/\frac{W}{p^2} = -\frac{C_0^2 \omega}{\widehat{\phi'}}$, left-hand side is solely a function of pressure while the right-hand side is solely a function of $(x,y,t)$. This implies $\frac{d^2 W}{dp^2}/\frac{W}{p^2} = -\frac{C_0^2 \omega}{\widehat{\phi'}}=\lambda$. According to Sturm-Liouville theory, $\lambda$ is a series of non-repeated integers (i.e., $\lambda=1,2,3$). We can further rewrite $\lambda$ as $\frac{C_0^2}{C^2}$ therefore  $C^2=\frac{C_0^2}{\lambda}$.  


:::{note}
While the linearization of {eq}`shallow_water_2` significantly simplifies the problem, the trade-off also leads to a few issues. The linear momentum equation doesn't allow the existence of energy cascade into new scales (i.e., filament structure of vorticity won't happen). Such a feature is very important in redistributing the lower troposphere moisture.  
:::


With the separation of variables, we can rewrite the governing equation

```{math}
:label: shallow_water_linear3
\begin{align}
\frac{\partial \widehat{u}}{\partial t} -\beta y\widehat{v} & = -\frac{\partial \widehat{\phi'}}{\partial x}\\
\frac{\partial \widehat{v}}{\partial t} +\beta y\widehat{u} & = -\frac{\partial \widehat{\phi'}}{\partial y}\\
\frac{\partial \widehat{\phi'}}{\partial t} + C^2 \frac{\partial \widehat{u}}{\partial x} +\frac{\partial \widehat{v}}{\partial y} &= 0 \\
\frac{d^2W(p)}{dp^2}+\frac{C_0^2}{C^2}W(p) &=0
\end{align}
```

where the first three equations only have horizontal structure and the last equation only has vertical structure, which makes both set of equations solvable. We will start with the horizontal structure and followed by the discussion of vertical structure. 


### Horizontal Wave Solutions 
To further simplify the problem, we can non-dimentionalize each variable by subtracting non-dimentional parameters. 


```{math}
:label: shallow_water_linear4
\begin{align}
(x,y) & = \sqrt {\frac{C}{\beta}}(x^{\text{non}},y^{\text{non}}) \\
t     & = \frac{t^{\text{non}}}{\sqrt{\beta C}} \\
u     & = C u^{\text{non}} \\
\phi^{'}  & = C^2 \phi^{\text{non}} 
\end{align}
```

where the superscript $\text{non}$ represents the non-dimensional variable. With such scaling, we can isolate the independent parameters and the yielded equations can be written as follows (we drop the superscript for simplification)

```{math}
:label: shallow_water_linear_no_dimension1
\begin{align}
\frac{\partial u}{\partial t} -yv & = -\frac{\partial \phi}{\partial x}\\
\frac{\partial v}{\partial t} +yu & = -\frac{\partial \phi}{\partial y}\\
\frac{\partial \phi}{\partial t} + (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) &= 0   \\
\end{align}
```

One advantage of using such non-dimensional equation is that we can apply it to a wide range of parameter space. For example, as long as the shallow water assumption holds, regardless of Rossby number, {eq}`shallow_water_linear_no_dimension1`. We can further employ the following assumption to separate the meridional and zonal structure of waves. 

```{math}
:label: shallow_water_linear_no_dimension2
\begin{align}
u &= U(y)e^{ikx-i\omega t} \\
v &= V(y)e^{ikx-i\omega t} \\
\phi &= \Phi(y)e^{ikx-i\omega t} \\
\end{align}
```

Substitute {eq}`shallow_water_linear_no_dimension2` back into {eq}`shallow_water_linear_no_dimension1` and elimating the terms of $U(y)$ and $\Phi(y)$. One can have the another differential equation which is used to describe the meridional structure of tropical waves. 

```{math}
:label: y_differential_equation
\frac{d^2}{dy^2} V^2 + (\omega^2-k^2+\frac{k}{\omega}-y^2)V = 0
```

Here we implement a boundary condition, where $V(y)$ tapers toward 0 at $y= \pm \infty$, to solve the above equations. Another interesting condition happens when we implement the long-wave approximation ($v$ vanish, because the wave is so long that its zonal boundaries (meridional wind happens) is barely observed in a limited domain), which leads to 

```{math}
:label: shallow_water_linear_no_dimension3
\begin{align}
\frac{\partial u}{\partial t} & = -\frac{\partial \phi}{\partial x}\\
yu & = -\frac{\partial \phi}{\partial y}\\
\frac{\partial \phi}{\partial t} + (\frac{\partial u}{\partial x}) &= 0   \\
\end{align}
```

and the corresponding differential equation can be written as 
```{math}
:label: Kelvin_wave_equation
\frac{\partial ^2 u}{\partial t^2}-\frac{\partial ^2 u }{\partial x^2} = 0
```

{eq}`y_differential_equation` and {eq}`Kelvin_wave_equation` leads to two different dispersion relationships, where wave's amplitude is conserved over the characteristic lines. 

```{math}
:label: dispersion
\begin{align}
\frac{\omega}{k} &= \mp 1 \\
\omega^2-k^2+\frac{k}{\omega} &= 2m+1 \text{   where } m\in[0,1,2,\cdots]\\
\end{align}
```

One should notice that the minus sign in the first equation doesn't really exist due to physical constrain. (Readers will try to figure out this point in the homework assignment). $2m+1$ in the second equation comes from the Sturm-Liouville theorem. 

Observing {eq}`dispersion`, it's not hard to find that the first equation is non-dispersive (i.e.,  wave length won't influence the propagation speed). In addition, it represents a pure gravity wave, which can be proved by calculating the PV in {eq}`shallow_water_linear_no_dimension3` (I will leave the practice to the reviewer). For the second equation, we can category the terms into three groups. (1) Rossby wave dominated regimes (2) Gravity wave dominated regimes, and (3) in-between. 

#### Case 1: Rossby wave regime
The main balance happens between $-k^2+\frac{k}{\omega}\approx 2m+1$ (i.e., $\omega$ is really small). In this case, $\omega \approx \frac{k}{k^2+2m+1}$, which is similar to what we derived in barotropic vorticity equation. Indeed, when $\omega$ is small, it represents the low-frequency limit. 

#### Case 2: Inertia gravity wave 
The main balance happens between $\omega^2+k^2\approx 2m+1$ (i.e., $\omega$ is big). In such case, $\omega = \pm \sqrt{k^2+2m+1}$. It's not hard to find that when $k$ is really big, $\omega \approx \pm k$ indicating it's dominated by gravity wave propagation in both directions. (Readers will complete the discussion of small $k$ case). 

#### Case 3: Mixed Rossby gravity wave (Yanai wave)
There is a special case where $m=0$, then {eq}`dispersion` becomes $(\omega-k)(\omega^2+k\omega-1)=0$. The three roots of this equation are $\omega=k$, $\omega = -\frac{k}{2}+\sqrt{(\frac{k}{2})^2+1}$ and $\omega = -\frac{k}{2}-\sqrt{(\frac{k}{2})^2+1}$. $\omega=k$ apparently has a gravity wave-like behavior. $\omega = -\frac{k}{2}+\sqrt{(\frac{k}{2})^2+1}$ and $\omega = -\frac{k}{2}-\sqrt{(\frac{k}{2})^2+1}$ can be further categorized into three cases for discussion (eastward propagation, westward propagation with small k and westward propagation with large k). 