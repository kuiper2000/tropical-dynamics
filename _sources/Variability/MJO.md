(MJO)=
# Week 8: Madden-Julian oscillation


The Madden-Julian Oscillation (MJO) is a distinctive phenomenon that occurs in the tropical Indian Ocean and the Western Pacific. It is characterized by large-scale atmospheric circulation patterns (with zonal wave numbers ranging from 1 to 4) and operates on intraseasonal timescales, lasting between 20 to 90 days. Due to its intraseasonal nature, it is often referred to as tropical intraseasonal variability. The MJO's signal is so strong that it can even be detected in raw-sounding data and water vapor channels. {numref}`FIG5-1` below illustrate the power spectrum of tropical sea-level pressure from various sounding stations , as presented in the {cite}`Madden1971-ao` paper. A notable peak around 40-50 days is a clear indicator of the MJO. Also, {numref}`FIG5-2` shows an observed MJO event located at the western Pacific. 


```{figure} ../tropical-dynamics-figures/MJO_power_spectrum.png
---
name: FIG5-1
---
:width: 700px
The power spectrum of tropical intraseasonal variability
```

```{figure} ../tropical-dynamics-figures/latest72hrs.gif
:width: 700px
---
name: FIG5-2
---
The total precipitable water of an MJO event in wester Pacific from MIMIC-TPW. 
```



The MJO has significant global effects. It plays a key role in modulating the variability of the summer monsoon in South and Southeast Asia, influencing related extreme events such as tropical cyclones and droughts. Additionally, its interaction with the subtropical jet stream can generate quasi-stationary Rossby waves, allowing wave energy to propagate toward North America, creating the Pacific-North America (PNA) pattern. This, in turn, impacts winter storms and the polar vortex. Despite the considerable global influence of the MJO, modeling its dynamics has posed challenges for the scientific community since its discovery. Two critical aspects in understanding the MJO are its (1) frequency (or period) and (2) the evolution of its vertical structure from shallow to deep convection. We will explore these various characteristics and examine how they affect the ability to accurately model the MJO.


## Theories for MJO
### wave perspective
Two primary theories have been developed to explain the MJO: (1) wave dynamics and (2) the moisture mode framework. Each theory addresses different physical thresholds. Regarding wave dynamics, notable contributions were made in the 1990s by Dr. Bin Wang, Dr. Tim Li, and C.-P. Chang (from our department). The wave dynamics theory, often referred to as the wave school, originates from the shallow water equations, where the influence of gravity and Rossby waves is significant. For example:

```{math}
:label: MJO_wave
\begin{align}
\frac{\partial \text{Div}}{\partial t} &\sim O( \nabla \cdot (f+\zeta)\times \vec{u})  \\
\frac{\partial \zeta}{\partial t} &\sim O(\nabla \cdot (f+\zeta)\times \vec{u}) \\
\frac{\partial \phi_p}{\partial t} &\sim O(w\frac{N^2}{g}-\frac{Q}{c_p T}) \\
Q &\sim -Lw\frac{\partial \tilde{q}}{\partial z} \\
\end{align}
```

in {eq}`MJO_wave` $\text{Div}$, $\zeta$ and $\phi_p$ represent the divergence, vorticity, and temperature, respectively. The tendency terms in these equations have magnitudes similar to other terms, such as advection, indicating that the inertial effects are still significant. This suggests the system has not yet reached a full dynamical equilibrium, and wave processes continue to play a role. Additionally, diabatic heating is governed by the vertical redistribution of moisture, similar to the convective instability of the second kind (CISK).


### moisture mode framework
In the moisture mode framework, the assumption is that the internal adjustment of gravity waves has already reached equilibrium, leading to the removal of $\frac{\partial \text{Div}}{\partial t}$ from the system. Furthermore, under the weak temperature gradient approximation, the third equation simplifies as follows:

```{math}
:label: MJO_moisture_mode
\begin{align}
\frac{\partial \text{Div}}{\partial t} &\sim 0  \\
\frac{\partial \zeta}{\partial t} &\sim O(\nabla \cdot (f+\zeta)\times \vec{u}) \\
w\frac{N^2}{g}&=\frac{Q}{c_p T} \\
Q &\sim L\frac{dq}{dt}+Q_R(T,p,q,\text{Cloud})\\
\end{align}
```

A key distinction in the moisture mode framework is that the forcing is considered a prognostic variable rather than a diagnostic one. This means that the forcing evolves over time, rather than being predetermined. From the comparison between {eq}`MJO_wave` and {eq}`MJO_moisture_mode`, we can see that the treatment of both the internal wave processes and the forcing mechanism is crucial.

It's important to note that these two frameworks represent opposite ends of the parameter spectrum, with one being more wave-like and the other more akin to moist static energy (MSE) variability. A substantial amount of tropical variability lies somewhere between these two perspectives.


## Some Evidence from Wave Dynamics Perspective

From the perspective of wave dynamics, the MJO has long been regarded as a type of convectively coupled wave for several reasons. First, the circulation pattern associated with the MJO closely resembles that of a forced Kelvin wave, where the horizontal winds exhibit a Gill-like response (which we will explore in more detail later). Specifically, there are easterly winds to the east of the convection and westerly winds to the west. The eastern portion of the MJO shows similarities to observed convectively coupled Kelvin waves, while the western part is often characterized by twin cyclones. Moreover, the easterly phase is associated with shallower circulation (and lower cloud heights), while the westerly phase is linked to stratiform clouds, featuring a top-heavy heating profile.

The shallow circulation aids in recharging moist entropy into the free troposphere through shallow divergent flow, which exports lower dry static energy. Conversely, deep circulation efficiently discharges moist entropy from the column to other regions. This recharge-discharge process of column energy is crucial to MJO dynamics, as it directly influences the intensity and distribution of diabatic heating. This aspect has been largely absent from traditional wave dynamics theories, and we will delve deeper into this topic when discussing the Gill solution.

The vertical structure of the MJO, which tilts from east to west, strongly resembles the structure observed in equatorial Kelvin waves. In tropical wave theory, the presence of shallow circulation can be explained by momentum balance within a frictional planetary boundary layer. Upward motion at the top of this boundary layer induces shallow heating (around the 4th baroclinic mode), which in turn slows the propagation speed. During the convective phase, the maximum upward motion near 300 hPa, coupled with divergence in the upper troposphere, creates a first baroclinic mode structure. This process allows the column to export energy, which accelerates the propagation speed. The coupling between shallow and deep heating profiles results in a wave speed that falls between these two modes. However, even with shallow heating factored in, the phase speed remains faster than what is observed.

To better understand how these different components affect the characteristic timescales of the MJO, the following sections will explore these concepts in more detail.

```{figure} ../tropical-dynamics-figures/MJO_vertical_structure.png 
:width: 700px
---
name: FIG5-2
---
The vertical structure of MJO and the cross-section at westerly/easterly phase. From {cite}`Adames2015-rm`
```

## Some Evidence from Moisture Mode Perspective 
On the other hand, there are also evidences supporting the moisture mode framework. First, the weak temperature gradient approximation is more solid on intraseasonal time scales (i.e., $w\frac{N^2}{g}=\frac{Q}{c_p T}$ is a better approximation on intraseasonal timescales than wave time scales). Second, 




```{bibliography}
```

