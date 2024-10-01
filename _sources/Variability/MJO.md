(MJO)=
# Week 8: Madden-Julian oscillation
## Weak Temperature and Steady State Response 


MJO is a unique feature existing in the tropical Indian Ocean and the Western Pacific, which characterized by planetary scale circulation (zonal wave number around $1-4$ and ) intraseasonal timescales (20-90 days of period). That's why it's also called _tropical intraseasonal variability_. Its signal is so strong that we can even observe it from the raw sounding data. Figure below shows the power spectrum of tropical sea-level pressure collecting at various sounding sites from {cite}`Madden1971-ao` paper. The local maximum around periods of 40-50 days shows the finger-print of MJO. 


```{figure} ../tropical-dynamics-figures/MJO_power_spectrum.png
:width: 700px
---
name: FIG5-1
---
The power spectrum of tropical intraseasonal variability
```

The MJO has pronounced global impacts. It modulates South/Southeast Asia summer monsoon variability and the associated extremes such as tropical cyclone and drought. Its interation with subtropical jet can creates quasi-stationary Rossby wave, where the wave energy propagates into the North America (so-called Pacific-North America pattern, PNA) and influences the winter storm activity as well as the polar vortex. Although the global influenes from MJO are significant, the communities have struggled to model its underlying dynamics since its first discovery. Two key variables are especially important (1) frequency (or period) and (2) shallow-to-deep vertical structure evolution. We will begin with various aspects about the MJO. 


## Theory for MJO
There are two main theories for MJO: (1) wave dynamics and (2) moisture mode framework. Each of them satisfy different thresholds. For wave dynamics, the tendency associated with gravity/Rossby wave is non-neglegible and quasi-equilibrium i.e., 

```{math}
:label: MJO_wave
\begin{align}
\frac{\partial \text{Div}}{\partial t} \neq 0 \\
\frac{\partial \zeta}{\partial t} \neq 0 \\
\frac{\partial \phi_p}{\partial t} \neq 0 \\
\end{align}
```

where $\text{Div}$, $\zeta$ and $\phi_p$ are the divergence, vorticity and the temperature respectively. In wave dynamics   

## Wave Dynamics

From wave dynamics perspective, the MJO has long been considered as a type of convectively coupled wave for two reanson. First, MJO's circulation pattern is very similar to those found in the forced Kelvin wave, where the horizontal wind is characterized by a Gill-like response (we will talk a bit more about Gill solution later) with easterly/westerly to the east/west of convection. In addition, the easterly region is usually characterized by shallower circulation (as well as cloud height) while the westerly phase is characterized by stratiform (top-heavy heating profile). The shallow circulation can recharge moist entropy into the free troposphere due to shallow divergent flow (export lower dry static energy). On the contrary, the deep circulation can efficiently discharge the column moist entropy to other regions. Such recharge-discharge of column energy plays a key role in MJO dynamics, which in turn determines the diabatic heating intensity/distribution. This has been a missing part in the traditional wave dynamics theory. (more discussion will be provided in the section of a Gill solution). 

The tilting vertical structure from east to west closely mimics the observed structure in equatorial Kelvin waves. In tropical wave theory, such shallow circulation can be explained through the momentum balance in a frictional planetary boundary. The upward motion at the top of boundary layer can induce shallow heating (close to 4th baroclinic mode) which further slows the propagating speed. In the convective phase, the  maximum ascending motion around 300-hPa and upper troposphere divergence is similar to the first baroclinic mode structure, which exports the column energy but also makes the propagating speed faster. Thus, the coupling between shallow and deep heating profile can create a wave speed in between two modes. However, even with the inclusion of shallow heating, the phase speed is still much faster than the observed feature. To understand how each of these components influences the characteristic timescales of MJO, the following few sections will provide more in-depth discussion. 


```{figure} ../tropical-dynamics-figures/MJO_vertical_structure.png 
:width: 700px
---
name: FIG5-2
---
The vertical structure of MJO and the cross-section at westerly/easterly phase. From {cite}`Adames2015-rm`
```






```{bibliography}
```

