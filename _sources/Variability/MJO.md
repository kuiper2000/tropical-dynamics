(MJO)=
# Week 8: Madden-Julian oscillation


MJO is a unique feature existing in the tropical Indian Ocean and the Western Pacific, which is characterized by planetary scale circulation (zonal wave number around $1-4$ and ) intraseasonal timescales (20-90 days of period). That's why it's also called _tropical intraseasonal variability_. Its signal is so strong that we can even observe it from the raw-sounding data. The figure below shows the power spectrum of tropical sea-level pressure collecting at various sounding sites from {cite}`Madden1971-ao` paper. The local maximum of around 40-50 days shows the fingerprint of MJO. 


```{figure} ../tropical-dynamics-figures/MJO_power_spectrum.png
:width: 700px
---
name: FIG5-1
---
The power spectrum of tropical intraseasonal variability
```

The MJO has pronounced global impacts. It modulates South/Southeast Asia summer monsoon variability and the associated extremes such as tropical cyclones and drought. Its interaction with the subtropical jet can create quasi-stationary Rossby waves, where the wave energy propagates into North America (so-called Pacific-North America pattern, PNA) and influences the winter storm activity and the polar vortex. Although the global influences from MJO are significant, the communities have struggled to model its underlying dynamics since its first discovery. Two key characteristics are essential in the MJO problem (1) frequency (or period) and (2) shallow-to-deep vertical structure evolution. We will begin with various aspects of the MJO and look into how the details in these aspects influence the modeled behavior of MJO.  


## Theories for MJO
### wave perspective
There are two main theories for MJO: (1) wave dynamics and (2) moisture mode framework. Each of them satisfies different thresholds. For wave dynamics, the tendency associated with gravity/Rossby wave is non-negligible i.e., 

```{math}
:label: MJO_wave
\begin{align}
\frac{\partial \text{Div}}{\partial t} &\sim O( \nabla \cdot (f+\zeta)\times \vec{u})  \\
\frac{\partial \zeta}{\partial t} &\sim O(\nabla \cdot (f+\zeta)\times \vec{u}) \\
\frac{\partial \phi_p}{\partial t} &\sim O(w\frac{N^2}{g}-\frac{Q}{c_p T}) \\
Q &\sim -Lw\frac{\partial \tilde{q}}{\partial z} \\
\end{align}
```

where $\text{Div}$, $\zeta$ and $\phi_p$ are the divergence, vorticity and the temperature respectively. One can find that the tendency terms in {eq}`MJO_wave` have equivalent (i.e., the difference falls within an order) amplitude to other terms such as advection. The existence of inertial term  

### moisture mode framework
In the moisture mode framework, it is assumed the internal adjustment of gravity has reached a dynamical equilibrium (the first equation is discarded). Also, according to the weak temperature gradient approximation, the third equation is written as:   

```{math}
:label: MJO_moisture_mode
\begin{align}
\frac{\partial \text{Div}}{\partial t} &\sim 0  \\
\frac{\partial \zeta}{\partial t} &\sim O(\nabla \cdot (f+\zeta)\times \vec{u}) \\
w\frac{N^2}{g}&=\frac{Q}{c_p T} \\
Q &\sim L\frac{dq}{dt}+Q_R(T,p,q,\text{Cloud})\\
\end{align}
```

In the moisture mode framework, there is one additional component that makes it distinct from the traditional wave framework. The forcing itself is a prognostic variable rather than a diagnostic variable. From these two differences, we can conclude that (1) wave dynamics focus on the timescales where gravity waves haven't reached a dynamical equilibrium state and the moisture is a diagnostic variable (dominated by tge vertical moisture advection) (2) moisture mode theory focuses on the timescales where gravity waves ``have'' reached a dynamical equilibrium while the moisture is prognostic. 

One should notice that these two are sitting on either side of parameter spectrum (wave-like vs MSE variability-like) and there is a lot of tropical variability sitting in the middle ground.  


## Some Evidence from Wave Dynamics Perspective

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

