(Hadley_Cell)=
# Week 4 and 5: Hadley Cell 
## Held-Hou Model

Follow the discussion in previous two weeks, this week, we will look at how energy is residtributed from tropics to other latitudes aka Hadley Cell. The early work on general circulation is motivated by the observation from Halley (1686) and Hadley (1735). They concluded that the tropical trade wind is part of large-scale ciruclation which redistributes the latitudinal difference in solar energy. Even without the upper-troposphere observation during that period, they also predicted the upper troposphere westerly, which is quite impressive. 

The zonal averaged Hadley circulation is shown in {numref}`FIG3-1` (mass stream function). The existence of Hadley cell has two main purposes (1) redistributing the meridional potential energy, and (2) conserving the angular momentum. We will see this in a minute.  


```{figure} ../tropical-dynamics-figures/mass_stream_function.png
---
name: FIG3-1
---
The long-term averaged mass stream function. 
``` 

There are a few attempts to model the observed Hadley circulation (i.e., the strength and width). Held and Hou (1980) is one of the most famous. To understand how the angular momentum is balanced and energy is redistributed, we start with a zonal mean zonal momentum equation. The Earth's angular momentum is the combination of two components: (1) the tangential velocity from solid body rotation and (2) the wind relative to the Earth's motion. 

According to the definition of angular momentum of unit mass: 

```{math}
:label: high_school_angular_momentum
m = r \overrightarrow{u} = a \mathrm{cos}\phi (\Omega a \mathrm{cos}\phi + u)  
```

where $m$ is the Earth angular momentum, $a$ is the Earth's radius and $u$ is the zonal wind. The conservation of angular momentum across latitudes implies that $\frac{\partial m}{\partial y}=0$. Readers can derive the above expression from zonal momentum equation. 

```{math}
:label: zonal_angular_momentum_equation
\frac{\partial u}{\partial t} - (f+\overline{\zeta})\overline{v} + \overline{w}\frac{\partial \overline{u}}{\partial z} 
```

