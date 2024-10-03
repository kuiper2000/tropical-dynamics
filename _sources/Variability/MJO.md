(MJO)=
# Week 8: Madden-Julian oscillation


The Madden-Julian Oscillation (MJO) is a distinctive phenomenon that occurs in the tropical Indian Ocean and the Western Pacific. It is characterized by large-scale atmospheric circulation patterns (with zonal wave numbers ranging from 1 to 4) and operates on intraseasonal timescales, lasting between 20 to 90 days. Due to its intraseasonal nature, it is often referred to as tropical intraseasonal variability. The MJO's signal is so strong that it can even be detected in raw sounding data. The figure below illustrates the power spectrum of tropical sea-level pressure from various sounding stations, as presented in the {cite}`Madden1971-ao` paper. A notable peak around 40-50 days is a clear indicator of the MJO.


```{figure} ../tropical-dynamics-figures/MJO_power_spectrum.png
:width: 700px
---
name: FIG5-1
---
The power spectrum of tropical intraseasonal variability
```

The MJO has significant global effects. It plays a key role in modulating the variability of the summer monsoon in South and Southeast Asia, influencing related extreme events such as tropical cyclones and droughts. Additionally, its interaction with the subtropical jet stream can generate quasi-stationary Rossby waves, allowing wave energy to propagate toward North America, creating the Pacific-North America (PNA) pattern. This, in turn, impacts winter storms and the polar vortex. Despite the considerable global influence of the MJO, modeling its dynamics has posed challenges for the scientific community since its discovery. Two critical aspects in understanding the MJO are its (1) frequency (or period) and (2) the evolution of its vertical structure from shallow to deep convection. We will explore these various characteristics and examine how they affect the ability to accurately model the MJO.


## Theories for MJO
### wave perspective
There are two main theories for MJO: (1) wave dynamics and (2) moisture mode framework. Each of them satisfies different thresholds. For wave dynamics, some iconic work has been done by Dr. Bin Wang, Dr. Tim Li, and C.-P. Chang (in our department) back in 90s. The origins of the theory in _wave-school_ start from shallow water equation, where the tendency term associated with gravity/Rossby wave is non-negligible i.e., 

```{math}
:label: MJO_wave
\begin{align}
\frac{\partial \text{Div}}{\partial t} &\sim O( \nabla \cdot (f+\zeta)\times \vec{u})  \\
\frac{\partial \zeta}{\partial t} &\sim O(\nabla \cdot (f+\zeta)\times \vec{u}) \\
\frac{\partial \phi_p}{\partial t} &\sim O(w\frac{N^2}{g}-\frac{Q}{c_p T}) \\
Q &\sim -Lw\frac{\partial \tilde{q}}{\partial z} \\
\end{align}
```

in {eq}`MJO_wave` $\text{Div}$, $\zeta$ and $\phi_p$ are the divergence, vorticity and the temperature respectively. One can find that the tendency terms in {eq}`MJO_wave` have equivalent (i.e., the difference falls within an order) amplitude to other terms such as advection. The existence of an inertial term indicates the system itself _hasn't_ reached a dynamical equilibrium and therefore the internal timescales of waves still play a role. The diabatic heating is directly determined by the vertical redistribution of moisture (similar to the convective instability of the second kind, i.e., CISK processes).

### moisture mode framework
In the moisture mode framework, it is assumed the internal adjustment of the gravity wave has reached an equilibrium state (the first equation is discarded). Also, according to the weak temperature gradient approximation, the third equation is written as:    

```{math}
:label: MJO_moisture_mode
\begin{align}
\frac{\partial \text{Div}}{\partial t} &\sim 0  \\
\frac{\partial \zeta}{\partial t} &\sim O(\nabla \cdot (f+\zeta)\times \vec{u}) \\
w\frac{N^2}{g}&=\frac{Q}{c_p T} \\
Q &\sim L\frac{dq}{dt}+Q_R(T,p,q,\text{Cloud})\\
\end{align}
```

In the moisture mode framework, there is one additional component that makes it distinct from the traditional wave framework. The forcing itself is a prognostic variable rather than a diagnostic variable. From {eq}`MJO_wave` and {eq}`MJO_moisture_mode`, we can conclude that how we treated (1) the internal processes of wave and (2) forcing matters. 

One should notice that these two are sitting on either side of parameter spectrum (wave-like vs MSE variability-like) and there is a lot of tropical variability sitting in the middle ground.  


## Some Evidence from Wave Dynamics Perspective

From wave dynamics perspective, the MJO has long been considered as a type of convectively coupled wave for reasons. First, MJO's circulation pattern is very similar to those found in the forced Kelvin wave, where the horizontal wind is characterized by a Gill-like response (we will talk a bit more about Gill solution later) with easterly/westerly to the east/west of convection. The eastern part is similar to the observed convectively coupled Kelvin wave and the western part is usually characterized by twin cyclones. In addition, the easterly region is usually characterized by shallower circulation (as well as cloud height) while the westerly phase is characterized by stratiform (top-heavy heating profile). The shallow circulation can recharge moist entropy into the free troposphere due to shallow divergent flow (export lower dry static energy). On the contrary, the deep circulation can efficiently discharge the column's moist entropy to other regions. Such recharge-discharge of column energy plays a key role in MJO dynamics, which in turn determines the diabatic heating intensity/distribution. (This has been a missing part in the traditional wave dynamics theory and more discussion will be provided in the section of a Gill solution). 

The tilting vertical structure from east to west closely mimics the observed structure in equatorial Kelvin waves. In tropical wave theory, such shallow circulation can be explained through the momentum balance in a frictional planetary boundary. The upward motion at the top of the boundary layer can induce shallow heating (close to 4th baroclinic mode) which further slows the propagating speed. In the convective phase, the  maximum ascending motion around 300-hPa and upper troposphere divergence is similar to the first baroclinic mode structure, which exports the column energy but also makes the propagating speed faster. Thus, the coupling between shallow and deep heating profile can create a wave speed in between two modes. However, even with the inclusion of shallow heating, the phase speed is still much faster than the observed feature. To understand how each of these components influences the characteristic timescales of MJO, the following few sections will provide more in-depth discussion. 


```{figure} ../tropical-dynamics-figures/MJO_vertical_structure.png 
:width: 700px
---
name: FIG5-2
---
The vertical structure of MJO and the cross-section at westerly/easterly phase. From {cite}`Adames2015-rm`
```






```{bibliography}
```

