(QBO)=
# Week 9: Quasi-biennial Oscillation

Quasi-biennial Oscillation (QBO) is a quasi-periodic oscillation in the stratospheric, zonal-mean zonal wind, characterized by a period of around 28 months. It was first discovered by Richard and Reed in 1960. The signal is so regular that one can observe it without implementing any filtering (see figure below)

```{figure} ../tropical-dynamics-figures/QBO.png
---
name: FIG6-1
width: 700px
---

The zonal mean zonal wind in the stratospheric tropics. 
```

The QBO can efficiently change the static stability around the tropopause, which has strong modulation on the convective activity. For example, during the easterly phase of QBO, the convection of MJO is much stronger than the westerly phase of QBO. Along with the MJO global impacts, one can expect that the fingerprint of QBO can be observed in most places around the globe. 

Given the QBO is high up in the stratosphere, the corresponding dynamics is relatively simple. A theory of QBO was first proposed by Charney and Drazin (1968), {cite}`Lindzen1968-zv` and {cite}`Holton1972-ff`, which incorporates the well-known wave-mean flow interaction dynamics. According to their theory, the oscillation is due to the interaction of upward-propagating Kelvin and Yanai waves with the mean flow. More recently, it has been proposed that eastward and westward-propagating gravity waves play an important role and the low-frequency waves were de-emphasized. 

The various equatorially trapped waves discovered by {cite}`Matsuno1966-tv` (discussed {ref}`waves1`) are observed to produce wave energy which further propagates upward and downward. It is believed that the energy source comes from diabatic processes such as convection (Nitta 1972). While in {ref}`waves1`, a resting basic state is considered, it is very important to consider mean flow due to the _Doppler Shift_. 

## Doppler-Shift of wave and the deposition of wave energy 

Like other waves, the atmospheric waves also experience Doppler Shift as long as the relative motion happens between the observer and the source of waves. For example, the westward propagating wave looks more stationary when riding on a westerly mean flow (i.e., frequency $\sim \infty$) and the eastward propagating wave looks more transient in the same mean flow (i.e., complete multiple oscillations in a short period). The figure below demonstrates how Doppler Shift changes the frequency of a wave. The siren frequency from an ambulance sounds higher when the ambulance approaches the person listening to the sound and vice versa. In such cases, the observed wave number (or frequency) may exceed the critical value where the mean flow can sustain the corresponding wave propagation (i.e., you can still tell it's a wave!). At the end, the wave breaks and deposits momentum. We will walk through more details in the later of this week.  


```{figure} ../tropical-dynamics-figures/Doppler_Shift.jpeg
---
name: FIG6-2
width: 700px
---

An example of the Doppler Shift in a moving object. 
```


To understand how these waves influence the phase transition of QBO (zonal mean flow), we gonna start with something classic, the Ellassen Palm theory for mountain wave (one should notice that the upward motion from the boundary is similar to inhomogeneous topography at the lower boundary.) There are two important ingredients in this theory (1) the direction of momentum transport and (2) when the wave momentum is deposited.  

```{figure} ../tropical-dynamics-figures/Mountain_wave.png
---
name: FIG6-3
width: 700px
---

An example of how the mountain wave propagates Eastward($+C$)/Westward($-C$) and upward at the same time.  
```


### The direction of momentum transport 
In the panel of {numref}`FIG6-3`, we can find that when the inhomogeneity of the lower boundary exists (such as mountain/mass flux from the lower boundary), it will trigger gravity propagating eastward and westward. The eastward propagating waves (westerly to the mean state) generally transport westerly momentum upward due to their zonal height tilting (tilting eastward with height) and the same concept can be applied to the westward propagating waves. 

However, when the mean westerly exists, the zonal height tilting due to the westward propagating wave vanishes, which is not the case for the eastward propagating wave. Therefore, not only do the types of gravity waves matter but whether they are filtered by mean flow also matters. 

### When momentum is deposited 
One should notice that the presence of vertical momentum transport does not necessarily indicate the change in mean flow. Because when the input and output have an equivalent amount, then there is no acceleration/deceleration of mean flow. The necessary condition of the presence of vertical momentum flux convergence can be derived through the angular momentum and eddy kinematic energy equations. 

To have both equations, we will begin with the zonal momentum equation, 

```{math}
:label: QBO_zonal_momentum
\begin{align}
u_t + U u_x + w U_z + \frac{1}{\rho_0} p_x = 0 \\ 
w_t + U w_x + \frac{1}{\rho_0} p_z + g = 0 \\
u_x + w_z = 0
\end{align}
```

Assume a wave solution of 


```{math}
:label: QBO_wave_solution
\begin{align}
u = \hat{u}(z)e^{ikx-ikct} \\ 
w = \hat{w}(z)e^{ikx-ikct} \\ 
p = \hat{p}(z)e^{ikx-ikct} + p_0(z)\\ 
\end{align}
```
(where $\frac{1}{\rho_0}p_0(z)+g=0$)

substitute into {eq}`QBO_zonal_momentum` 
```{math}
:label: QBO_spectral_form
\begin{align}
ik (U-c) u + w U_z + ik\frac{1}{\rho_0} p = 0 \\ 
ik (U-c) w + \frac{1}{\rho_0} p_z = 0 \\
ik u + w_z = 0
\end{align}
```

Multiply each equation in {eq}`QBO_spectral_form` by (1) $u$ (2) $w$ and (3) $p$ respectively.  We have...

```{math}
:label: QBO_mechanical_energy
\begin{align}
& \begin{cases}
[ik \rho_0 (U-c) u^2 + ik up] & = -\rho_0 uw U_z \\ 
[ik \rho_0 (U-c) w^2 ]+ wp_z  & = 0 \\
[ik up + w_z p] & = 0
\end{cases} \\
\rightarrow & \underbrace{ik [E(U-c)+pu]}_{\text{the horizontal flux of mechanical energy }} + \underbrace{(pw)_z}_{\text{the vertical flux of mechanical energy}} = \underbrace{-\rho_0 uw U_z}_{\text{Shear production}} 
\end{align}
```

The above formula is similar to the _turbulence kinetic energy_ used in boundary layer dynamics. Integrate the last equation over zonal direction, we have the zonal mean mechanical energy 


```{math}
:label: QBO_mechanical_energy_zonal_mean
\underbrace{\int_{0}^{2\pi}(pw)_z dx}_{\text{the vertical flux of mechanical energy}} = \underbrace{\int_{0}^{2\pi} -\rho_0 uw (U-c)_z dx}_{\text{Shear production}} 
```


We can also link the above equation to the angular momentum equation. To achieve this, we multiply the first equation of {eq}`QBO_spectral_form` by $\rho_0$ (U-c) u + p and integrate it over zonal direction (which eliminates the terms associated with zonal gradient, i.e., "$ik$" terms ) We have 


```{math}
:label: QBO_mechanical_energy_zonal_mean2
\begin{align}
\int_0^{2\pi} pw dx = -(U-c) \rho_0 \int_0^{2\pi} uw dx 
\end{align}
```

{eq}`QBO_mechanical_energy_zonal_mean` and {eq}`QBO_mechanical_energy_zonal_mean2` suggests that $\rho_0 \int_0^{2\pi} uw dx =\text{const}$ when $U-c\neq 0$. i.e., no momentum will be deposited until the wave reaches the critical level (the level where $U=c$). At a critical level, we can find some analogs in our ambulance example. It corresponds to where ambulance and sound travel with the same speed and direction. Therefore, in a limited traveling length of the wave, we can observe a nearly infinite number of waves making the finite assumption of wave dynamics no longer hold and momentum is deposited. 

In addition, when $\int_0^{2\pi} uw$ is positive through a critical level (i.e., keeps transporting mechanical energy upward), it implies that $uw$ term must change signs above and below the critical level, which breaks the assumption of $\rho_0 \int_0^{2\pi} uw dx =\text{const}$. 

To solve the problem, all of the momentum _must_ be absorbed at the critical level. Booker and Bretherton (1967) provide a useful formula...  

```{math}
:label: QBO_mechanical_energy_zonal_mean3
\begin{cases}
\rho_0 \int_0^{2\pi} uw dx = A \\
\rho_0 \int_0^{2\pi} uw dx = -A e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}} 
\end{cases}
\begin{align}
& \text{ for z below critical level} \\
&  \text{ for z above critical level} 
\end{align}

```

{eq}`QBO_mechanical_energy_zonal_mean3` suggests that $\rho_0 \int_0^{2\pi} uw dx$ is constant below critical level. Right above the critical level, it needs to taper toward 0 in a short range of traveling distance. Given the Richardson number is always greater than 1, the second equation of {eq}`QBO_mechanical_energy_zonal_mean3` indeed satisfies the condition we need. {eq}`QBO_mechanical_energy_zonal_mean3` also suggests the momentum absorbed at the critical level is $A[1+e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}}]$. 

However, the absorption of momentum at a single level will lead to a shock-like signal and modeling difficulty. Also, we need to determine the sign of A to make all necessary conditions consistent. 

### A spectral solution of momentum deposition 
To circumvent the problem in the previous section, we can approach it with a spectral perspective of wave propagation. 
```{math}
:label: QBO_mechanical_energy_zonal_mean4
\begin{align}
\rho_0 \int_0^{2\pi} uw dx                       & = \int^{\infty}_{-\infty} f(c) dc = f(U) |dU| \\
(\rho_0 \int_0^{2\pi} uw dx)_k|^{z_0+}_{z_0-}/dz & = -A[1+e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}}]/dz = f(U_0) [1+e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}}]  |dU/dz| \\
& \text{ where $U_0$ is the critical mean flow for wave $k$}\\
\end{align}
```

It is assumed that the disturbance that transports momentum consists of a spectrum of waves with a continuous distribution of phase speed. For a limited range of phase speed $c\pm dc$, we can find a limited range of critical level ($U\pm dU$) which absorbs the momentum. Such momentum absorption only applies to a finite range of wave and keep the rest unaffected by the mean flow. (i.e., the second equation of {eq}`QBO_mechanical_energy_zonal_mean4`). The second equation of {eq}`QBO_mechanical_energy_zonal_mean4` comes from the fact that $f(U) |dU|$ has a unit of momentum change per unit time. Thus, the eddy momentum flux convergence (i.e., the rate change of mean flow due to eddy can be written as)

```{math}
:label: QBO_mechanical_energy_zonal_mean5
-\frac{d F_{WM}}{dz} = - f(U_0)  [1+e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}}]  |\frac{dU}{dz}|
```
where $\Delta t=|\frac{dU}{dz}|$. $F_{WM}$ is the eddy momentum flux (i.e., the amount of absorbed momentum) by wave with $U_0=c$ and $-\frac{d F_{WM}}{dz}$ is the corresponding momentum flux convergence at the critical level. 

Here, we will use two cases to analyze the sign of $f(U_0)$. The first case represents the westerly shear (westerly increases with height) and the second case represents the easterly shear. 

```{figure} ../tropical-dynamics-figures/QBO_dynamics.png
---
name: FIG6-4
width: 700px
---

Momentum fluxes in (a) westerly shear and (b) easterly shear
```

In the westerly shear, if $\int_0^{2\pi} pw dx$ remains positive, for the regions below the critical level (where $c>U$), $\int_0^{2\pi} uw dx $ must be positive and so does $f(U)^{z-}$ according to {eq}`QBO_mechanical_energy_zonal_mean4` and {eq}`QBO_mechanical_energy_zonal_mean2`. For the regions above where (where $U>c$), $\int_0^{2\pi} uw dx $ must be negative. This also makes $(f(U)^{z+}-f(U)^{z-})/dz = f(U_0)  [1+e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}}]  |\frac{dU}{dz}|$ negative according to {eq}`QBO_mechanical_energy_zonal_mean4`. Since $|\frac{dU}{dz}|$ is positive definite, $f(U)=f(U_0)  [1+e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}}]$ will be positive definite. 

Readers will also find $- f(U_0)  [1+e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}}]  |\frac{dU}{dz}|$ is equivalent to $- |f(U_0)|  [1+e^{-2\pi\sqrt{\mathbf{Ri}-\frac{1}{4}}}]  \frac{dU}{dz}$ by examining both westerly and easterly shear case. (HW)

## Lidzen and Holton Model 
Substitute the above conclusion back into zonal momentum equation, the yielded {eq}`LH68_01` looks identical to a simple advection model, where the _downward_ propagation of zonal wind is influenced by the flux convergence of eddy momentum. The downward propagation speed is determined by the intensity of momentum flux convergence. 

```{math}
:label: LH68_01
\rho\frac{\partial \overline{u}}{\partial t} + \overline{w}\frac{\partial \overline{u}}{\partial z} = -f(U)\frac{\partial \overline{u}}{\partial z}
```

In this model, some simple setup is given. First, we assume there is no mean vertical motion $\overline{w}$. Second, the downward motion due to the momentum flux convergence only happens over a limited range of wind speed. i.e., 

```{math}
:label: LH68_02
\begin{cases}
f(u) &= \text{constant} \\  
f(u) &= 0
\end{cases}
\begin{align}
\text{for $-c<u<c$, where $c=5m/s$}\\
\text{for $U>=|c|$, where $c=5m/s$} 
\end{align}
```

Using the example below, at the initial time, only regions between 37.5 and 22.5 km can experience westerly acceleration while other regions are left out. Therefore, for wind speed fall within this ranges will keep accelerating until it hit the upper bound of westerly. The same idea can be applied to mean easterly. 

At the end, one can observe the downward propagation of zonal wind anomaly. Currently, the scientific community doesn't have a good parameterization of bulk eddy. (given the challenges of higher-order closure problem) Therefore, the downward propagation is determined through the observation.  

```{figure} ../tropical-dynamics-figures/QBO_model.png
---
name: FIG6-5
width: 700px
---
The simulated QBO evolution in the westerly shear using {eq}`LH68_01` 
```

The study of QBO also shows the importance of properly parameterized gravity wave drag given its global scale of change in zonal wind. In earlier years, most GCM experienced severe cold bias over the polar regions and stratosphere. One reason is that most gravity waves don't deposit enough momentum to change the mean flow. 



```{bibliography}
```

