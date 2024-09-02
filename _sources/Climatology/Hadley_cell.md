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


To understand how the angular momentum is balanced and energy is redistributed, we start with a zonal mean zonal momentum equation. The Earth's angular momentum is the combination of two components: (1) the tangential velocity from solid body rotation and (2) the wind relative to the Earth's motion. 

According to the definition of angular momentum of unit mass: 

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


If we assume a dynamical equilibrium (i.e., Earth's rotation is spinning up or down and so is the atmosphere), zonally symmetric (axisymmetric) atmosphere, where the vertical advection of zonal wind is negligible. Then we have 

```{math}
:label: high_school_angular_momentum2
(f+\overline{\zeta})\overline{v} = 0 
```

If we assume $v$ is not zero, then it leads to $(f+\overline{\zeta})=0$. An alternative way to write such relationship 

```{math}
:label: high_school_angular_momentum3
2\Omega \mathrm{sin}\phi = \frac{1}{a}\frac{\partial \overline{u}}{\partial \phi}-\frac{\overline{u}\mathrm{tan}\phi}{a}
```

By assuming the zonal wind at $\phi=0$, we can solve the above ODE by integrating over $\phi$ i.e., 

```{math}
:label: high_school_angular_momentum4
\overline{u}(\phi) = \frac{\Omega a \mathrm{sin}^2\phi}{\mathrm{cos}\phi}
```
{eq}`high_school_angular_momentum4`


