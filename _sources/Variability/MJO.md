(MJO)=
# Week 8: Madden-Julian oscillation


The Madden-Julian Oscillation (MJO) is a distinctive phenomenon that occurs in the tropical Indian Ocean and the Western Pacific. It is characterized by large-scale atmospheric circulation patterns (with zonal wave numbers ranging from 1 to 4) and operates on intraseasonal timescales, lasting between 20 to 90 days. Due to its intraseasonal nature, it is often referred to as tropical intraseasonal variability. The MJO's signal is so strong that it can even be detected in raw-sounding data and water vapor channels. {numref}`FIG5-1` below illustrate the power spectrum of tropical sea-level pressure from various sounding stations , as presented in the {cite}`Madden1971-ao` paper. A notable peak around 40-50 days is a clear indicator of the MJO. Also, {numref}`FIG5-2` shows an observed MJO event located at the western Pacific from unfiltered column water vapor.  


```{figure} ../tropical-dynamics-figures/MJO_power_spectrum.png
---
name: FIG5-1
width: 700px
---

The power spectrum of tropical intraseasonal variability
```

```{figure} ../tropical-dynamics-figures/latest72hrs.gif
---
name: FIG5-2
width: 700px
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
Q &\sim -Lw\frac{\partial \tilde{q}}{\partial z}(\text{SST}) \\
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

The shallow circulation aids in recharging moist entropy into the free troposphere through shallow divergent flow, which exports lower dry static energy but imports moist-laden air from the boundary layer. Conversely, deep circulation efficiently discharges moist entropy from the column to other regions. This recharge-discharge process of column energy is crucial to MJO dynamics, as it directly influences the intensity and distribution of diabatic heating. This aspect has been largely absent from traditional wave dynamics theories (or more focus on the recharge of DSE from wave dynamics perspective). We will delve deeper into this topic when discussing the Gill solution.

The vertical structure of the MJO, which tilts from east to west, strongly resembles the structure observed in equatorial Kelvin waves. In tropical wave theory, the presence of shallow circulation can be explained by momentum balance within a frictional planetary boundary layer. Upward motion at the top of this boundary layer induces shallow heating (around the 4th baroclinic mode), which in turn slows the propagation speed. During the convective phase, the maximum upward motion near 300 hPa, coupled with divergence in the upper troposphere, creates a first baroclinic mode structure. This process allows the column to export energy, which accelerates the propagation speed. The coupling between shallow and deep heating profiles results in a wave speed that falls between these two modes. However, even with shallow heating factored in, the phase speed remains faster than what is observed.


```{figure} ../tropical-dynamics-figures/MJO_vertical_structure.png 
---
name: FIG5-3
width: 700px
---
The vertical structure of MJO and the cross-section at westerly/easterly phase. From {cite}`Adames2015-rm`
```

### The boundary layer momentum balance
In the boundary layer, the momentum equation can be written as 

```{math}
:label: wave_boundary_layer
\begin{align}
\frac{\partial u_B}{\partial t} -\beta y v_B &= -\frac{\partial \phi}{\partial x}-E u_B \\
\frac{\partial v_B}{\partial t} +\beta y v_B &= -\frac{\partial \phi}{\partial y}-E v_B \\
\end{align}
```

where the subscript _B_ represents boundary layer. One should notice that we don't include the thermodynamics equation in the boundary layer equation given that the diabatic processes associated with shallow cumulus is weaker those those found in the free troposphere (not a perfect but workable assumption). Another important feature is that the boundary and lower free troposphere share the same pressure gradient, which can be attributed to the fast adjustment in the planetary boundary. Leveraging the same assumption (fast adjustment) for scale analysis, we can also drop the tendency term due to the fast adjustment, leading to   

```{math}
:label: wave_boundary_layer2
\begin{align}
-\beta y v_B &= -\frac{\partial \phi}{\partial x}-E u_B \\
\beta y v_B &= -\frac{\partial \phi}{\partial y}-E v_B \\
\end{align}
```

{eq}`wave_boundary_layer2` is a diagnostic equation, i.e., as long as we know the pressure gradient force in the free troposphere, we know $u_B$ and $v_B$ in the boundary layer. The boundary layer in turn feedback to the free troposphere through {eq}`boundary_layer` ($\tilde{w}$). The shallow vertical motion can trigger shallow mode of heating and recharge the column MSE.  

From {eq}`wave_boundary_layer2` we also know the steady motion of boundary layer convergence on an equatorial beta plane can be written as 

```{math}
:label: wave_boundary_layer3
\begin{align}
D=-E(\underbrace{\nabla^2 \phi}_{\text{pressure gradient force}}+\underbrace{\beta u}_{\text{Ekman Pumping}} +\underbrace{\beta^2 yv/E}_{\text{Sverdrup Balance}})/(E^2+\beta^2y^2)
\end{align}
```

In large-scale motion (where the pressure gradient is small), the second term is the dominant forcing where the lower-troposphere westerly/esterly is characterized by boundary layer top convergence/divergence. Such boundary-top convergence/divergence will further lead to the downward/upward motion in the lower troposphere. In the off-equatorial region, Severdrup balance (negative planetary PV advection balanced by stretching) is the main cause of the lower-troposphere vertical motion. 

One key take-home message is that the boundary-layer top easterly collocated with boundary ascending motion, which is a quarter phase leading the mid-troposphere ascending motion. Thus, as long as the surface friction exists, we can always see such phase difference between the boundary layer ascending motion and the mid-troposhere ascending motion.   

:::{note}
Severdrup balance has also been applied to explain the stationary Rossby wave response such as the monsoon-desert mechanism proposed by Dr. Brian Hoskins back in 1980s. See {cite}`Rodwell1996-bz`
:::


## Some Evidence from Moisture Mode Perspective 
Three necessary conditions for a moisture mode variability happens, (1) the diabatic heating is dominated by the moisture-related processes (2) weak temperature gradient (3) Precipitation anomaly is highly correlated with moisture anomaly. 

According to these criteria, there are evidences supporting the moisture mode framework, especially over the Indo-Pacific regions. First, the weak temperature gradient approximation (i.e., $w\frac{N^2}{g}=\frac{Q}{c_p T}$ is a good approximation on intraseasonal timescales (even a more solid approximation than applying it to wave time scales). This indicates that the dynamics are nearly in a state of equilibrium. Second, the moistening/drying in the lower troposphere ($\sim$700-hPa) has an equivalent contribution to the residual of the column process (i.e., vertical moisture advection balanced by the physical processes), suggesting what happens in horizontal directions is equally important. Precipitation anomaly is in generally highly correlated to column moisture if it's linearized in intraseasonal timescales. 

Overall, MJO has satisfied all three criteria in the region of the Indian Ocean, and it meets either one or two criteria over other basins. Details can be found in XXX. 


### Gill Solution (the forced wave solution)
The Gill solution is a steady state, forced response of the wave solution. To acquire the steady-state wave solution, we adopt {eq}`shallow_water_linear3` 


```{math}
:label: Gill_solution
\begin{align}
\epsilon u -\beta yv & = -\frac{\partial \phi}{\partial x}\\
\epsilon v +\beta yu & = -\frac{\partial \phi}{\partial y}\\
\epsilon \phi + C_m^2 (\frac{\partial u}{\partial x} +\frac{\partial v}{\partial y}) &= Q_{m} \\
\frac{Q}{c_p T} &=\sum_{m} Q_{m} W_m(p) \\ 
\end{align}
```

where $W_m(p)$ is the eigen basis of vertical normal mode decomposition and $Q_{m}$ is the coefficient of the corresponding vertical mode. (see {eq}`vertical_normal_mode2`). Thus, the first step to solve Gill solution is decomposing the diabatic heating into various normal modes, where each mode corresponds to different equivalent depths (gravity wave speed). Then, the derived coefficients (i.e., $Q_m$ in the third equation of {eq}`Gill_solution`) of each mode are functions of longitude and latitude, which will be used to solve the horizontal wave solution (i.e., first 3 equations in {eq}`Gill_solution`) 

Gill's solution has been applied to explain the observed circulation on a beta plane with a resting basic state. Some real-world cases include the monsoon-desert relationship, MJO circulation, and North American low-level jet. In addition, it's the backbone of the moisture mode framework. In moisture mode, heating is dominated by the variation of moisture, which is redistributed by circulation. The circulation can further be diagnosed using Gill solution. 

Following the discussion in the first week, one can find that the moisture mode starts from an energetic perspective and emphasizes on the role of moisture in determining the intraseasonal variation of moist static energy. Thus, the intrinsic timescales of moisture mode variability is determined by the redistribution of MSE rather than the adjustment of DSE. 

More details can be found in {cite}`Adames2016-cm`. 

### Gross-moist stability 




```{bibliography}
```

