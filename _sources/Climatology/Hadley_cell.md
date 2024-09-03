(Hadley_Cell)=
# Week 4 and 5: Hadley Cell 
## Held-Hou Model

Follow the discussion in previous two weeks, this week, we will look at how energy is residtributed from tropics to other latitudes. The early work on general circulation is motivated by the observation from Halley (1686) and Hadley (1735). They concluded that the tropical trade wind is part of large-scale ciruclation which redistributes the latitudinal difference in solar energy. Even without the upper-troposphere observation during that period, they also predicted the upper troposphere westerly, which is quite impressive. 

The zonal averaged Hadley circulation is shown in {numref}`FIG3-1` (mass stream function). The existence of Hadley cell has two main purposes (1) redistributing the meridional potential energy, and (2) conserving the angular momentum. We will see this in a minute. According to these two purposes, a few interesting features shown in {numref}`FIG3-1` lead to intriguing conclusions. First, the maximum ascending motion centers around $5^o$N suggesting the inter-hemisphere imbalance of incoming energy. Therefore, the Hadley cell will transport the energy from one hemisphere to the other hemisphere. Second, the descending branch is around 30$^o$N, which is relatively asymmetric about the hemisphere. (Indeed, if one flips {numref}`FIG1-1` with respect to the equator, a nearly identical incoming shortwave radiation can be found in both hemispheres). This implies the energy balance within this region is quite symmetric.   


```{figure} ../tropical-dynamics-figures/mass_stream_function.png
---
name: FIG3-1
---
The long-term averaged mass stream function. 
``` 

To model how the energy and angular momentum are redistributed via Hadley cell, there were a few attempts traced back to 1980s. Held and Hou (1980) is one of the pioneers in this field. Held and Hou (1980) model can be briefly summarized by {numref}`FIG3-2`. 


```{figure} ../tropical-dynamics-figures/Held-and-Hou.png
---
name: FIG3-2
---
Held and Hou (1980) Hadley Cell model
```

Held and Hou model is a two-layer model mimicking the observed Hadley cell circulation which has poleward flow in the upper layer transporting higher potential energy to other latitudes and equatorward return flow at the lower layer. In the middle, an interface of radiative-equilibrium potential temperature represents the difference in received energy over latitudes, which can be formulated as 

```{math}
:label: HH_model_theta_e
\begin{align}
\frac{d \theta}{d t} & = \frac{\theta_E-\theta}{\tau} \\
\theta_E (\phi) & = \theta_0 -\frac{1}{3}\Delta \theta(3 \mathrm{sin}^2\phi-1)  = \theta_0 -\Delta \theta \frac{y^2}{a^2} \\
& \text{given that $\phi \approx \mathrm{sin}\phi \approx \frac{y}{a}$}\\
\end{align}
```
where the heating rate of a given column is determined by the energy imbalance between $\theta_E$ and $\theta$ as well as the relaxation time (first equation of {eq}`HH_model_theta_e`). The former represents the radiative equilibrium profile of received energy and the latter is the atmospheric temperature. $\theta_E$ follows the second Legendre mode which has a maximum at the equator and minimum at the pole. 


To understand how the angular momentum is balanced and energy is redistributed, we start with a zonal mean zonal momentum equation. 

The Earth's angular momentum is the combination of two components: (1) the tangential velocity from solid body rotation and (2) the wind relative to the Earth's motion. According to the definition of angular momentum of unit mass: 

```{math}
:label: high_school_angular_momentum
m = r \overrightarrow{u} = a \mathrm{cos}\phi (\Omega a \mathrm{cos}\phi + u)  
```

where $m$ is the Earth angular momentum, $a$ is the Earth's radius and $u$ is the zonal wind. The conservation of angular momentum across latitudes implies that $\frac{\partial m}{\partial y}=0$. Readers can derive the above expression from zonal momentum equation. 

```{math}
:label: zonal_angular_momentum_equation
\frac{\partial \overline{u}}{\partial t} - (f+\overline{\zeta})\overline{v} + \overline{w}\frac{\partial \overline{u}}{\partial z} = -\frac{1}{\mathrm{cos}^2\phi} \frac{\partial }{\partial \phi} (\mathrm{cos}^2\phi \overline{u'v'})-\frac{\partial \overline{u'w'}}{\partial z}
```

One can derive {eq}`zonal_angular_momentum_equation` by writing non-linear momentum advection into the rotational acceleration and the convergence/divergence of kinetic. (see Vallis 2017). In {eq}`zonal_angular_momentum_equation`, prime and bar represent zonal mean and eddy components respectively. $(f+\overline{\zeta})\overline{v}$ is the acceleration of zonal wind by the torque of Coriolis force and relative vorticity. $\overline{w}\frac{\partial \overline{u}}{\partial z}$ is the vertical advection zonal wind. The two double prime terms on the right are eddy covariance terms (it will be discussed in a moment).


If we assume a dynamical equilibrium (i.e., Earth's rotation is not spinning up or down), zonally symmetric (axisymmetric) atmosphere, where the vertical advection of zonal wind is negligible. Then we have 

```{math}
:label: high_school_angular_momentum2
(f+\overline{\zeta})\overline{v} = 0 
```

and further assume $v$ is not zero, then it leads to $(f+\overline{\zeta})=0$. An alternative way to write such relationship is

```{math}
:label: high_school_angular_momentum3
2\Omega \mathrm{sin}\phi = \frac{1}{a}\frac{\partial \overline{u}}{\partial \phi}-\frac{\overline{u}\mathrm{tan}\phi}{a}
```

By assuming the zonal wind at $\phi=0$, we can solve the above ODE by integrating over $\phi$ i.e., 

```{math}
:label: high_school_angular_momentum4
\overline{u}(\phi) = \frac{\Omega a \mathrm{sin}^2\phi}{\mathrm{cos}\phi} \rightarrow \overline{u}_M=\frac{\Omega}{a}y^2  \\ 
\text{given that $\phi \approx \mathrm{sin}\phi \approx \frac{y}{a}$ and 0 surface wind}
```
{eq}`high_school_angular_momentum4` shows how upper-level wind varies with latitude. Given that both hydrostatic and geostrophic balances hold in the tropics (as long as the aspect ratio is small), we can further link {eq}`high_school_angular_momentum4` to {eq}`HH_model_theta_e`. i.e., ($f\frac{\partial \overline{u}}{\partial z}= -\frac{g}{\theta_0}\frac{\partial \theta}{\partial y}$)

```{math}
:label: thermal_wind
\begin{align}
f\frac{\partial \overline{u}}{\partial z} &= \frac{u_M}{H} \approx 2\Omega \frac{y}{a}= \frac{\Omega}{aH}y^2 \\
\frac{\partial \theta_M}{\partial      y} &= -\frac{2\Omega^2\theta_0}{a^2g H}y^3 
\end{align}
```

The second equation represents the necessary meridional temperature gradient to maintain an angular-momentum conserved thermal wind. Integrating over y, we have 

```{math}
:label: angular_momentum_theta
\begin{align}
\theta_M = \theta_{M0}-\frac{\Omega^2\theta_0}{2a^2gH}y^4
\end{align}
```

{eq}`angular_momentum_theta` can be summarized as the following figure. 

```{figure} ../tropical-dynamics-figures/angular_momentum_energy_balance.png
---
name: FIG3-3
---
The angular-momentum conserved $\theta$ (black) versus radiative-equilibrium $\theta$ (orange). 
``` 

In {numref}`FIG3-3`, if we move along the black curves, the total angular momentum is conserved. On the other hand, the orange curve represents the conservation of radiative energy. Therefore, if the observed temperature follows the black curve, the Hadley should redistribute energy from the tropics (the regions where orange curve surpasses black curve) to the energy sink of subtropics (black curve surpasses orange curve). 



