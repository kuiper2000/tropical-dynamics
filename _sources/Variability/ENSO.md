(ENSO)=
# Week 10: El Ni\~no Southern Oscillation

The El Nino Southern Oscillation (ENSO) is one of the most important variability on seasonal to interannual time scales. It has profound global influence and even influences across timescales. For example, the spring tornado frequency over the continental US, the North Western subtropical high variability, tropical cyclone frequency, and global warming patterns are all shaped by ENSO. These are only a few. ENSO is characterized by 2-7 years, making it highly predictable on seasonal timescales. Given its global impact and long predictability, successfully modeling ENSO has been one of the Holy Grails of the climate community. 

```{figure} ../tropical-dynamics-figures/El_Nino_Image.png
---
name: FIG10-1
width: 700px
---
The 2016 El Nino events from. Credits: NOAA_PSL. 
```

## Background 
The earliest record of El Nino Southern oscillation was by a sailor in Peru back in 1892, long before its theory was developed. It wasn't until 1969, that Jacob Bjerknes proposed a conceptual model for positive feedback in tropical air-sea interaction. The ENSO theory is more completed when Drs. Mark Cane and Stephan Zebiak proposed a series of modeling frameworks, which incorporates (1) atmosphere, (2) ocean mixed layer, and (3) thermocline. Building on their prototype ENSO model, two schools of theory were proposed to explain the observed variability: (1) Wave theory and (2) Recharge-discharge of warm water volume. Among the wave theories, the Delayed oscillator by Drs. Max Suarez, Paul S. Schopf, David Battisti, and Anthony Hirts is the most famous. The recharge-discharge oscillator is proposed by Dr. Fei-Fei Jin (also Tim Li in the following year).  


## Warm Water Volume (Recharge-Discharge Mechanism)
We will go through some fundamental processes consisting of ENSO. One should notice that both mechanisms (recharge-discharge and delayed oscillator) are simplified versions of the Zebiak and Cane and McCreary models for better physical interpretation given the more complex structure of these two intermediate models. 


Recharge-Discharge Mechanism is summarized in the diagram below. 


```{figure} ../tropical-dynamics-figures/ENSO_diagram.png
---
name: FIG10-2
width: 700px
---
The (known) feedback processes in ENSO.  
```

(1)+(2) is the well-known Bjerknes feedback. When the central or eastern Pacific is characterized by warm SST (where is climatologically cold), it will favor the development of convection over that region, i.e., weakened Walker circulation. The enhanced westerly will advect warm water from the western Pacific to the eastern Pacific. 

(3)+(4)+(5) is the negative feedback through subtropical cells. When the enhanced westerly happens over the equatorial surface, it increases cyclonic wind stress curl at the subtropical regions. To balance the reduced relative vorticity, a poleward geostrophic current brings negative planetary vorticity to this region. Such geostrophic currents also bring warm water to high latitudes and shallow the equatorial thermocline. An opposite process happens during the La Nina year. 

:::{note}
Sverdrup balance was first used to explain the existence of Western Boundary current. It was then modified by Dr. Biran Hoskins to explain the extension of subtropical high and monsoon gyre. (See Hoskins and Rodwell). Dr. Fei-Fei Jin was working with Dr. Brian Hoskins  
:::


(6)+(7) The enhanced westerly during the El Nino year also suppresses the equatorial upwelling by Ekman pumping. (i.e., Ekman feedback)

(8)+(9) On the other hand, the deepened thermocline also makes it harder to bring cold water below the thermocline (i.e., thermocline feedback). 

Considering all of these processes as a whole is the well-known recharge-discharge oscillator by Jin (1997). One should notice, that the existence of _oscillation_ should involve at least one negative feedback. The diagram above also indicates that the entire process can be reduced to a two-variable system, where SST and thermocline depth are only predictors. The atmospheric-related processes are in a steady state due to their transient timescales compared to ocean processes. Thus, the entire system can be formulated as follows: 


```{math}
:label: Recharge-discharge
\begin{align}
\frac{dT}{dt}=CT+Dh-\epsilon T^3 \\
\frac{dh}{dt}=-ET-R_h
\end{align}
```

Details of each term will be provided in the final project. 


## Wave School (the oscillators)
### Delayed oscillator
The earliest ENSO oscillator is proposed by McCreary (1983) (Fig. {numref}`McCreary`), which coupled an oceanic shallow water model with a diagnostic atmospheric component. 

```{figure} ../tropical-dynamics-figures/ENSO.gif
---
name: McCreary
width: 700px
---
The simulated ENSO evolution is based on McCreary model. 
```

McCreary hypothesized that the reflection of oceanic Rossby waves can help generate the interannual variability of SST. By emphasizing on the western boundary, Suarez and Schopf (1988) proposed the prototype of the delayed oscillator. Zebiak and Cane then coupled an atmospheric Gill model with a 1.5-layer SST model (mixed layer + thermocline dynamics) becoming the first who successfully predicted the ENSO. Battisti and Hirst (1989) used such an intermediate model to propose the well-know delayed oscillator mechanism. 

Different from the Recharge-discharge oscillator, the delayed oscillator focuses on the transient dynamics and explains the change of warm water volume through a _wave_ lens. 

The delayed oscillator can be formulated as follows: 

```{math}
:label: western_oscillator
\begin{align}
\frac{dT}{dt} = AT-BT(t-\eta)-\epsilon T^3
\end{align}
```
 where T represents the SST anomaly in the _equatorial eastern Pacific_. The first term on the right-hand side is the Bjerknes feedback between the ocean and the atmosphere. The second term represents the negative feedback due to the wave reflection on the western boundary. The cubic term is a higher-order damping (doesn't rule out the periodicity). Overall, the Delayed oscillator focuses on the western boundary processes but can capture most of the variability found in the Zebiak and Cane model. 


### Western-Pacific Oscillator
Supported by observational evidence and other modeling studies, Weisberg and Wang (1997) proposed the Western-Pacific oscillator. One can consider it a more complicated version of the delayed oscillator but more focused on the role of the western Pacific. The Western-Pacific oscillator can be summarized in the following figure.  

```{figure} ../tropical-dynamics-figures/delayed_oscillator.jpeg
---
name: FIG10-4
width: 700px
---
The western Pacific oscillator in the Delayed oscillator. From Wang (2018): A review of ENSO theories. 
```

and the corresponding equation 

```{math}
:label: western_oscillator
\begin{align}
\frac{dT}{dt} &= a \tau_1 +b_2\tau_2(t-\delta)-\epsilon T^3 \\ 
\frac{dh}{dt} &= -c \tau_1(t-\lambda)-R_hh \\ 
\frac{d \tau_1}{dt} &= dT -R_{\tau_1}\tau_1 \\
\frac{d \tau_2}{dt} &= eh -R_{\tau_2}\tau_2 \\
\end{align}
```


Where, T, h, $\tau_1$ and $\tau_2$ are illustrated in {numref}`FIG10-4` 

We first start with equatorial westerly. The equatorial westerly ($\tau_1$) at Nino 4 driven by tropical convection (like MJO) will drive a downwelling Kelvin wave (warm anomaly) which further propagates eastward to increase Nino 3 temperature. This is represented in the first term of equation 1. At the same time, the atmospheric Gill response also drives twin cyclones off the equator. The twin cyclone can induce oceanic upwelling due to Ekman pumping, which further expands westward to the Nino 6 regions (as indicated in the time-delayed term of $\tau_1(t-\lambda)$). The cold SST then triggers surface high, which induces equatorial easterly at Nino 5. Unlike the equatorial westerly, the equatorial easterly will induce upwelling Kelvin wave, reversing the warm anomaly pattern. Such equatorial easterly has been the main focus over the past few years (i.e., Dr. Wayne Lee's PhD work). 

### Advective-reflective oscillator
The advective-reflective oscillator is very similar to the previous two oscillators except looking into the wave reflection on both sides of the basin.

## The Primitive Equation-based ENSO models
While these oscillators/theories focus on different aspects of ENSO, they all got some success in simulating ENSO due to the strong air-sea coupling, i.e., omitting one or a few variables won't hurt the main variability as long as the omitted variables can be represented by other processes. 

For a more completed picture, the McCreary (1983) and Zebiak and Cane (1987) models have incorporated all of the ingredients above using the primitive equations. Both models are a combination of (1) one oceanic shallow water model and (2) a diagnostic atmosphere component. Here we use McCreary Model as an example. 

For the oceanic component, 

```{math}
:label: McCreary_ocean
\begin{align}
& u_t -\beta y v + p_x = F + \nu_h \nabla^2 u \\
& v_t +\beta y u + p_y = G + \nu_h \nabla^2 v \\
& p_t + c^2 (u_x+v_x) = 0 
\end{align}
```

The ocean model {eq}`McCreary_ocean` is nearly identical to the atmospheric Gill model. The only difference is the existence of wind stress forcing (i.e., $F$ and $G$). The wind stress forcing is the main process in triggering the oceanic waves.   


```{math}
:label: McCreary_atmosphere
\begin{cases}
\tau_h &= \tau_{0h}X(x-x_h)Y_h(y) \\
\tau_w &= \tau_{0w}X(x-x_h)Y_w(y) \\
\tau_b &= \tau_{0b}X(x-x_h)Y_b(y) \\
\end{cases}
```


```{math}
:label: McCreary_atmosphere2
\begin{align}
X(x)   & = cos (2\pi x/D) \\ 
Y_h(y) & = 1/2[1-cos (2\pi y/\lambda_h)] \\
Y_w(y) & = 1/2[1+cos (2\pi y/\lambda)] \\
Y_b(y) & = 1/2[1-cos (2\pi y/\lambda)] 
\end{align}
\begin{cases}
|x|\leq D/4 \\
|y|\leq \lambda_h \\ 
|y|\leq \lambda/2 \\
|y|\leq \lambda   \\
X = Y_h = Y_w = Y_b = 0, \text{otherwise} 
\end{cases}
```

