(waves1)=
# Week 6 and 7: Tropical Wave Theory
## Shallow Water Model

Aside from the zonal mean circulation, eddies significantly influence tropical climatology. We refer to the results from {ref}`scale_analysis` in this section for the discussion of this week. In a hydrostatic, and low-latitude atmosphere (weak temperature gradient), the governing equation can be written as. 


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


By integrating the third equation from the surface (or the boundary top) to the interface $\eta$ (refer to {numref}`FIG4-1`) and using the reduced gravity concept from {eq}`reduced_gravity2`, we can rewrite {eq}`shallow_water_1`.

```{math}
:label: shallow_water_2
\begin{align}
\frac{\partial u}{\partial t} + [u\frac{\partial u}{\partial x}+ v\frac{\partial u}{\partial y}]-\beta yv & = -g'\frac{\partial \eta}{\partial x}\\
\frac{\partial v}{\partial t} + [u\frac{\partial v}{\partial x}+ v\frac{\partial v}{\partial y}]+\beta yu & = -g'\frac{\partial \eta}{\partial y}\\
w(z=\eta) & = \frac{d \eta}{dt} = -\overline{\eta} (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) + w_0  
\end{align}
```

Here, $\eta$ represents the interface location (see {eq}`reduced_gravity2`), and $w(z)$ denotes the rate of displacement change, expressed as $\frac{d \eta}{dt}$.

Furthermore, based on quasi-equilibrium principles, we understand that diabatic heating and cooling are counterbalanced by the adiabatic cooling and heating from vertical motion, as indicated by {eq}`thermodynamics`. Consequently, a forcing term can be introduced into the final equation, expressed as:

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

This uses the geopotential height to rewrite the previous equation. The last line of {eq}`shallow_water_3` shows that diabatic heating-driven vertical motion almost balances horizontal flow convergence/divergence. This applies to both external (with unreduced gravity) and internal modes (with reduced gravity).

:::{note} It's important to note that $w_0$ in {eq}`shallow_water_2` isn't necessarily zero. For non-normal flow at the lower boundary, we need $w(x, y, z = h_B, t) = u\frac{\partial h_B}{\partial x} + v\frac{\partial h_B}{\partial y}$, assuming $\frac{\partial h_B}{\partial t} = 0$ at the boundary. This results in:

```{math}
:label: boundary_layer
w_0=\underbrace{u\frac{\partial h_B}{\partial x}+v\frac{\partial h_B}{\partial y}}_{\text{momentum advection}}+\underbrace{h_B(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y})}_{\text{mass convergence from lower bound}}
```
{eq}`boundary_layer` shows that normal velocity arises from two processes: (1) momentum advection of $h_B$ by $u$ and $v$, and (2) mass convergence/divergence at $h_B$. 
:::

To simplify {eq}`shallow_water_3` into an internal wave equation, we can calculate the difference between two vertical layers, introducing a counteracting force from the upper fluid. Alternatively, a more useful form of the final equation in {eq}`shallow_water_3` is:


```{math}
:label: hydrostatic_thermo
\frac{d \phi_p}{dt} +\sigma \omega = -\frac{Q}{ c_p T} \alpha
```

Here, $\phi_p = -\alpha$, $\omega = \overline{\rho} g w$, and $\sigma = \frac{R}{p \rho g} (\Gamma_d - \Gamma) = \frac{N^2}{\rho^2 g^2} = \frac{C_0}{p^2}$, where $C_0$ is the external gravity wave speed. {eq}`hydrostatic_thermo` represents the hydrostatic thermodynamic equation, derived by rewriting $d \mathrm{ln}\theta$ from {eq}`primitive` with $\frac{c_p}{R} \mathrm{ln}\alpha$. This form makes it easier to solve vertical motion since pressure velocity weighting is already factored in.

### Conservation of Shallow Water PV
{eq}`shallow_water_2` also conserves potential vorticity. To demonstrate this, we can take the curl of the first equation, resulting in:


```{math}
:label: shallow_water_4
\begin{align}
\frac{d\zeta}{dt} = \frac{\partial \zeta}{\partial t} + u\frac{\partial \zeta}{\partial x}+v\frac{\partial \zeta}{\partial y} &= -(\zeta+f)(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) \\
& = \frac{\zeta+f}{\overline{\eta}}\frac{d\eta}{dt} \text{  or}\\
& \rightarrow  \frac{1}{\zeta+f}\frac{d\zeta+f}{dt} = \frac{1}{\overline{\eta}}\frac{d\eta}{dt}
\end{align}
```

The last equation of {eq}`shallow_water_4` indicates that the rate of change in absolute vorticity is proportional to the rate of change in height, implying that $\frac{d}{dt}\left( \frac{\zeta + f}{\eta} \right) = 0$. (Note that the expression derived in the final equation of {eq}`shallow_water_2` imposes a stricter constraint on potential vorticity conservation.) Furthermore, {eq}shallow_water_4 can also be derived from {eq}`shallow_water_2` using the Kelvin circulation theorem (homework).



### Equatorial Wave Equations  
Although {eq}`shallow_water_2` is significantly simplified, it remains unsolvable due to its nonlinear nature. To address this, we can begin with a linear approximation or assume a resting basic state ($u = \overline{u} + u'$ where $\overline{u}$ is independent of time, or $\overline{u} = 0$ for a resting basic state) with a boundary condition of $w_0 = 0$. This simplifies the equations as follows:

```{math}
:label: shallow_water_linear
\begin{align}
\frac{\partial u}{\partial t} -\beta yv & = -\frac{\partial \phi'}{\partial x}\\
\frac{\partial v}{\partial t} +\beta yu & = -\frac{\partial \phi'}{\partial y}\\
\frac{\partial u}{\partial x} +\frac{\partial v}{\partial y} + \frac{\partial \omega}{\partial p}& = 0 \\
\frac{\partial \phi_p'}{\partial t} + \frac{C_0^2}{p^2} \omega &= 0  
\end{align}
```

where $g'\eta = C_0^2$, with $C_0$ representing the gravity wave speed. This set of equations is the well-known formulation derived by Matsuno (1966) and is a solvable system.

As {eq}`shallow_water_linear` contains four unknowns and four equations, it is inherently solvable. Additionally, the equation has a separable form, allowing us to assume a solution of the form $u, v, \phi = {\widehat{u}, \widehat{v}, \widehat{\phi}}(x, y, t) \frac{dW(p)}{dp}$ and $w = \widehat{w}(x, y, t) W(p)$ to separate the horizontal and vertical structures. Substituting this back into {eq}`shallow_water_linear`, we obtain:



```{math}
:label: shallow_water_linear2
\begin{align}
\frac{\partial \widehat{u}}{\partial t} -\beta yv & = -\frac{\partial \widehat{\phi'}}{\partial x}\\
\frac{\partial \widehat{v}}{\partial t} +\beta yu & = -\frac{\partial \widehat{\phi'}}{\partial y}\\
\frac{\partial \widehat{u}}{\partial x} +\frac{\partial \widehat{v}}{\partial y} + \widehat{\omega}& = 0 \\
\frac{\partial \widehat{\phi'}}{\partial t}\frac{d^2W}{dp^2} + \frac{C_0^2}{p^2} \widehat{\omega}W(p) &= 0  
\end{align}
```

Rearranging the last equation of {eq}`shallow_water_linear2` as $\frac{d^2 W}{dp^2} / \frac{W}{p^2} = -\frac{C_0^2 \omega}{\widehat{\phi'}}$, the left-hand side is solely a function of pressure, while the right-hand side depends only on $(x, y, t)$. This leads to $\frac{d^2 W}{dp^2} / \frac{W}{p^2} = -\frac{C_0^2 \omega}{\widehat{\phi'}} = -\lambda$. According to Sturm-Liouville theory, $\lambda$ represents a series of distinct integers (i.e., $\lambda = 1, 2, 3$). We can rewrite $\lambda$ as $\frac{C_0^2}{C^2}$, so that $C^2 = \frac{C_0^2}{\lambda}$.

:::{note} While linearizing {eq}`shallow_water_2` significantly simplifies the problem, it introduces limitations. The linear momentum equation prevents the energy cascade into new scales (i.e., it does not capture the filamentary structure of vorticity). This feature is crucial for redistributing moisture in the lower troposphere. 
:::

With the separation of variables, we can rewrite the governing equations as:


```{math}
:label: shallow_water_linear3
\begin{align}
\frac{\partial \widehat{u}}{\partial t} -\beta y\widehat{v} & = -\frac{\partial \widehat{\phi'}}{\partial x}\\
\frac{\partial \widehat{v}}{\partial t} +\beta y\widehat{u} & = -\frac{\partial \widehat{\phi'}}{\partial y}\\
\frac{\partial \widehat{\phi'}}{\partial t} + C^2 \frac{\partial \widehat{u}}{\partial x} +\frac{\partial \widehat{v}}{\partial y} &= 0 \\
\frac{d^2W(p)}{dp^2}+\frac{C_0^2}{C^2}W(p) &=0
\end{align}
```

Here, the first three equations describe the horizontal structure, while the last equation describes the vertical structure, making both sets of equations solvable. We will begin by solving for the horizontal structure, followed by a discussion of the vertical structure.



### Horizontal Wave Solutions 
To further simplify the problem, we can non-dimensionalize each variable by introducing non-dimensional parameters:


```{math}
:label: shallow_water_linear4
\begin{align}
(x,y) & = \sqrt {\frac{C}{\beta}}(x^{\text{non}},y^{\text{non}}) \\
t     & = \frac{t^{\text{non}}}{\sqrt{\beta C}} \\
u     & = C u^{\text{non}} \\
\phi^{'}  & = C^2 \phi^{\text{non}} 
\end{align}
```

Here, the superscript $\text{non}$ indicates the non-dimensional variable. With this scaling, we can isolate the independent parameters, and the resulting equations can be written as follows (for simplicity, the superscripts are dropped):

```{math}
:label: shallow_water_linear_no_dimension1
\begin{align}
\frac{\partial u}{\partial t} -yv & = -\frac{\partial \phi}{\partial x}\\
\frac{\partial v}{\partial t} +yu & = -\frac{\partial \phi}{\partial y}\\
\frac{\partial \phi}{\partial t} + (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}) &= 0   \\
\end{align}
```

One advantage of using this non-dimensional form is its broad applicability to different parameter spaces. For instance, as long as the shallow water assumption holds, this set of equations remains valid regardless of the Rossby number. We can further simplify by assuming a separable structure for the meridional and zonal components of the waves:


```{math}
:label: shallow_water_linear_no_dimension2
\begin{align}
u &= U(y)e^{ikx-i\omega t} \\
v &= V(y)e^{ikx-i\omega t} \\
\phi &= \Phi(y)e^{ikx-i\omega t} \\
\end{align}
```

Substituting {eq}`shallow_water_linear_no_dimension2` back into {eq}`shallow_water_linear_no_dimension1` and eliminating terms involving $U(y)$ and $\Phi(y)$, we derive a differential equation that describes the meridional structure of tropical waves:



```{math}
:label: y_differential_equation
\frac{d^2}{dy^2} V + (\omega^2-k^2-\frac{k}{\omega}-y^2)V = 0
```

We can solve this equation by applying the boundary condition that $V(y)$ tapers off to zero as $y \to \pm \infty$ and assume $V(y)=\widehat{V(y)}e^{-y^2}$. This assumption can convert {eq}`y_differential_equation` into a typical Hermite differential equation, where $\widehat{V(y)}=H_m(y)$ is so-called Hermite polynomial. 

Another special case arises when we use the long-wave approximation (where $v$ vanishes due to the wave's length being so large that the zonal boundaries, where meridional winds occur, are negligible in a limited domain). This leads to:

```{math}
:label: shallow_water_linear_no_dimension3
\begin{align}
\frac{\partial u}{\partial t} & = -\frac{\partial \phi}{\partial x}\\
yu & = -\frac{\partial \phi}{\partial y}\\
\frac{\partial \phi}{\partial t} + (\frac{\partial u}{\partial x}) &= 0   \\
\end{align}
```

The corresponding differential equation is:

```{math}
:label: Kelvin_wave_equation
\frac{\partial ^2 u}{\partial t^2}-\frac{\partial ^2 u }{\partial x^2} = 0
```

Equations {eq}`y_differential_equation` and {eq}`Kelvin_wave_equation` yield two distinct dispersion relationships where the wave amplitude is conserved along characteristic lines:

```{math}
:label: dispersion
\begin{align}
\frac{\omega}{k} &= \mp 1 \\
\omega^2-k^2-\frac{k}{\omega} &= 2m+1 \text{   where } m\in[0,1,2,\cdots]\\
\end{align}
```

Note that the negative sign in the first equation is unphysical due to certain constraints (this is an exercise left for the reader). The term $2m+1$ in the second equation is derived from Sturm-Liouville theory.

Looking at {eq}`dispersion`, we see that the first equation is non-dispersive (i.e., the wavelength does not affect the wave propagation speed). This represents a pure gravity wave, which can be proven by calculating the potential vorticity in {eq}`shallow_water_linear_no_dimension3` (another exercise for the reader). The second equation can be broken down into three distinct regimes: (1) Rossby wave-dominated, (2) gravity wave-dominated, and (3) mixed regimes.

#### Case 1: Rossby wave regime
In this case, the dominant balance is $-k^2-\frac{k}{\omega} \approx 2m+1$, meaning $\omega$ is small. Therefore, $\omega \approx -\frac{k}{k^2+2m+1}$, similar to the result from the barotropic vorticity equation. This regime represents the low-frequency limit (i.e., the time scales dominated by the Earth's rotation)

#### Case 2: Inertia-gravity wave regime
Here, $\omega^2-k^2 \approx 2m+1$, meaning $\omega$ is large. In this case, $\omega = \pm \sqrt{k^2+2m+1}$. For large $k$ (larger than $2m+1$ but about the same scale of $\omega$), $\omega \approx \pm k$, indicating the wave is dominated by gravity wave propagation. The reader can continue the analysis for small $k$ and will find it slightly tilting westward due to the existence of $2m+1$, indicating this wave can still feel the rotation of the earth given its planetary scale circulation. 

#### Case 3: Mixed Rossby-gravity wave (Yanai wave)
A special case occurs when $m=0$, leading to the simplified dispersion relation $(\omega+k)(\omega^2+k\omega-1)=0$. The three roots are $\omega=-k$, $\omega = \frac{k}{2}+\sqrt{\left(\frac{k}{2}\right)^2+1}$, and $\omega = \frac{k}{2}-\sqrt{\left(\frac{k}{2}\right)^2+1}$. The root $\omega=-k$ behaves like a gravity wave, while the other two can be categorized into three cases (eastward propagation, westward propagation with small $k$, and westward propagation with large $k$). In the eastward-propagating case, it tapers toward the inertia gravity wave's dispersion relationship. In low wave number to westward propagating cases, it behaves like a Rossby wave. Therefore, this special case is called mixed Rossby-gravity wave, which has the characters of both waves. In the real world, MRG is an important embryo for cyclone genesis. 

These three cases are summarized and shown in the figure below. 

```{figure} ../tropical-dynamics-figures/Dispersion_Relationship.png
---
name: FIG4-2
---
The dispersion relationship of tropical waves. 
```

An easy way to understand wave dispersion is to plot the Hovmoller diagram of variables of interest in a longitude-time plane. (i.e., the figure below {numref}`FIG4-2`). If we track along a specific feature (let's say ), the slope of that feature on the Hovmoller diagram represents one single dot on the space-time plot i.e., $(\omega,k)$ is a specific pair. One can also take Fourier transform twice (one along time and the other along longitude coordinate) of the Hovmoller diagram, the yielded result is also the space-time diagram.     

```{figure} ../tropical-dynamics-figures/TRMM_DYNAMO.png
---
name: FIG4-3
width: 400px
---
The meridional averaged precipitation from TRMM product during the DYNAMO field campaign. 
```

### Vertical Normal Mode
An important thing is {eq}`shallow_water_linear4` is scaled by $C$ rather than $C_0$ indicating that we still need to determine $C$ to close the problem. According to the separation of variable, we know $\frac{C_0^2}{C^2}=\lambda$ or $C = \sqrt{\frac{C_0^2}{\lambda}}$. $\lambda$ on the other hand represents the "ID" of a wave's vertical structure. To see what $\lambda$ stands for, we will complete the derivation of vertical normal mode. 


```{math}
:label: vertical_normal_mode1
\begin{align}
\frac{d^2 W(p)}{dp^2}+\lambda \frac{W^2}{p^2} = 0 
\end{align}
```

Here we implement two boundary conditions. At $p=p_s$ and $p=p_t$, $W(p)=0$, which represents the ridge-lid boundary condition. (one should notice that the lower boundary is not necessarily 0 given the frictional convergence plays a certain role in wave dynamics).  

```{math}
:label: vertical_normal_mode2
\begin{align}
& \frac{d^2 W(p)}{dp^2}+\lambda \frac{W^2}{p^2} = 0 \\
\end{align}
```
{eq}`vertical_normal_mode2` is an Euler differential equation, which can be solved by assuming a solution form of $W(p)=\sum_{n=0}^{\infty} C_n p^{n+r}$. Using this series solution, we can have 

```{math}
:label: vertical_normal_mode3
\begin{align}
W(p) = Ap^{\frac{1}{2}+\sqrt{\frac{1}{4}-\lambda}}+Bp^{\frac{1}{2}-\sqrt{\frac{1}{4}-\lambda}}
\end{align}
```

Substitute two boundary conditions into {eq}`vertical_normal_mode3` we have 

```{math}
:label: vertical_normal_mode4
\begin{cases} 
W(p_s) &= Ap_s^{\frac{1}{2}+\sqrt{\frac{1}{4}-\lambda}}+Bp_s^{\frac{1}{2}-\sqrt{\frac{1}{4}-\lambda}} \\
W(p_t) &= Ap_t^{\frac{1}{2}+\sqrt{\frac{1}{4}-\lambda}}+Bp_t^{\frac{1}{2}-\sqrt{\frac{1}{4}-\lambda}} 
\end{cases}
```

which implies $\sqrt{\frac{1}{4}-\lambda}=\frac{i \pi m}{\mathrm{ln}(\frac{p_s}{p_t})}$ and $B=-A p_t^{2\sqrt{\frac{1}{4}-\lambda}}$ (readers will complete the steps in the HW). The vertical structure solution can therefore be written as 

```{math}
:label: vertical_normal_mode5
W(p) = Ap^{\frac{1}{2}+\sqrt{\frac{1}{4}-\lambda}}-A p_t^{2\sqrt{\frac{1}{4}-\lambda}}p^{\frac{1}{2}-\sqrt{\frac{1}{4}-\lambda}}
```

It is obvious that two basis functions constitute the vertical structure solution i.e., $p^{\frac{1}{2}}p^{\frac{i\pi m}{\mathrm{ln}(\frac{p_s}{p_t})}}$ and $p^{\frac{1}{2}}p^{\frac{-i\pi m}{\mathrm{ln}(\frac{p_s}{p_t})}}$. Results are plotted in the figure below 

```{figure} ../tropical-dynamics-figures/Vertical_Normal_Mode.png
---
name: FIG4-4
width: 700px
---
The vertical normal mode for m=1,2,3,4
```

Recall that $\lambda = \frac{C_0^2}{C^2}$, also $\frac{1}{4}-\lambda=-\frac{\pi^2 m^2}{(\mathrm{ln}(\frac{p_s}{p_t}))^2}$. Two equations lead to $C^2=\frac{C_0^2}{\frac{1}{4}+\frac{\pi^2m^2}{(\mathrm{ln}(\frac{p_s}{p_t}))^2}}$. 


It's not hard to find that bigger $m$ corresponds to a more complicated vertical structure and slower propagating speed. Indeed, higher $m$ usually represents are more stratified fluid and therefore the density difference between two adjacent layers is smaller as well. Small density difference leads to more reduction in gravitational acceleration. 


If one substitutes the non-dimension parameters with typical tropical conditions, we will have $C_0\sim 76\frac{m}{s}$, $C_1\sim 50\frac{m}{s}$, $C_2\sim 26.2\frac{m}{s}$,... The first baroclinic structure has a wave speed around 50$\frac{m}{s}$ which is much faster than the observed tropical wave phase speed. Madden-Julian oscillation (one of the dominant tropical variability), has a phase speed around 5-10$\frac{m}{s}$ closer to the fourth baroclinic mode. MJO is, however, dominated by the first baroclinic mode structure. This inconsistency drives the community to study MJO in two different directions (1) the MJO is governed by wave dynamics with various vertical structures or (2) the MJO's characteristic time sccales are determined by something else other than the internal timescales of free-wave solution. We will have a more in-depth discussion next week.  


## Convective Coupled Wave 

While {eq}`dispersion` seems over-simplified (linear, resting basic states...), the observation is however closely tied to the analytical solution of tropical wave solution. The figure below shows the symmetric components of tropical OLR (averaged over 15$^{\circ}$N-15$^{\circ}$S). One can find two peaks stand out, one follows the dispersion relationship of Kelvin wave and the other one follows the equatorial Rossby waves. Given that OLR is a good proxy of convection, {numref}`FIG4-5` suggests that the convectively coupled tropical waves to some degree can be considered as a forced response of free-wave solution, where the forcing time scales are either much shorter or equivalent to the internal time scales of free-wave solutions.   


```{figure} ../tropical-dynamics-figures/Dispersion_Relationship_OLR.png
---
name: FIG4-5
width: 700px
---
The OLR power spectrum superimposed with the analytical dispersion relationship. The three curves (from top to bottom) for each type of wave correspond to different equivalent depths (50,25 and 12 meter). 
```

Another interesting but less addressed feature is the maximum around $k=1-4$ and $\text{CPD}\leq 0.05$ (i.e., period lower than 20 days). This is the famous Madden-Julian oscillation. Its distinct spectrum suggests the linear wave theory may not be able to explain the corresponding dynamics. Finding the governing equation for modeling MJO has been the Holy Grail in the scientific community of tropical large-scale dynamics. We will go through more debates in between different communities in the next week. 