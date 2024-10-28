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

Given the QBO is high up in the stratosphere, the corresponding dynamics is relatively simple. A theory of QBO was first proposed by Charney and Drazin (1968), Lidzen and Holton (1968) and Holton and Lidzen (1972), which incorporates the well-known wave-mean flow interaction dynamics. According to their theory, the oscillation is due to the interaction of upward-propagating Kelvin and Yanai waves with the mean flow. More recently, it has been proposed that eastward and westward-propagating gravity waves play an important role and the low-frequency waves were de-emphasized. 

The various equatorially trapped waves discovered by Matsuno (1966) (discussed {ref}`waves1`) are observed to produce wave energy which further propagates upward and downward. It is believed that the energy source comes from diabatic processes such as convection (Nitta 1972). While in {ref}`waves1`, a resting basic state is considered, it is very important to consider mean flow due to the _Doppler Shift_. 

## Doppler-Shift of wave and the deposition of wave energy 

Like other waves, the atmospheric waves also experience Doppler Shift as long as the relative motion happens between the observer and the source of waves. For example, the westward propagating wave looks more stationary when riding on a westerly mean flow (i.e., frequency $\sim \infty$) and the eastward propagating wave looks more transient in the same mean flow (i.e., complete multiple oscillations in a short period). The figure below demonstrates how Doppler Shift changes the frequency of a wave. The siren frequency from an ambulance sounds higher when the ambulance approaches the person listening to the sound and vice versa. In such cases, the observed wave number (or frequency) may exceed the critical value where the mean flow can sustain the corresponding wave propagation (i.e., you can still tell it's a wave!). At the end, the wave breaks and deposits momentum. We will walk through more details in the later of this week.  


```{figure} ../tropical-dynamics-figures/Doppler_Shift.jpeg
---
name: FIG5-2
width: 700px
---

An example of the Doppler Shift in a moving object. 
```


To understand how these waves influence the phase transition of QBO (zonal mean flow), we gonna start with something classic, the Ellassen Palm theory for mountain wave (one should notice that the upward motion from the boundary is similar to inhomogeneous topography at the lower boundary.) There are two important ingredients in this theory (1) the direction of momentum transport and (2) when the wave momentum is deposited.  

```{figure} ../tropical-dynamics-figures/Mountain_wave.png
---
name: FIG5-3
width: 700px
---

An example of how the mountain wave propagates Eastward($+C$)/Westward($-C$) and upward at the same time.  
```


### The direction of momentum transport 
In the panel of {numref}`FIG5-3`, we can find that when the inhomogeneity of lower boundary exists (such as mountain/mass flux from the lower boundary), it will trigger gravity propagating eastward and westward. The eastward propagating waves (westerly to the mean state) generally transport westerly momentum upward due to its zonal height tilting (tilting eastward with height) and the same concept can be applied to the westward propagating waves. 

However, when the mean westerly exists, the zonal height tilting due to the westward propagating wave vanishes, which is not the case for the eastward propagating wave. Therefore, not only the types of gravity wave matters but whether they are filtered by mean flow also matters. 

### When momentum is deposited 
One should notice that the presence of vertical momentum transport does not necessarily indicate the change in mean flow. Because when the input and output have an equivalent amount, then there is no acceleration/deceleration of mean flow. The necessary condition of the presence of vertical momentum flux convergence can be derived through the angular momentum and eddy kinematic energy equations. 

To have both equation, we will begin with zonal momentum equation, 

```{math}
:label: QBO_zonal_momentum
\begin{align}
u_t + U u_x + w U_z + \frac{1}{rho_0} p_x = 0 \\ 
\end{align}
```





```{bibliography}
```

