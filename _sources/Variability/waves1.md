(waves1)=
# Week 6 and 7: Tropical Wave Theory
## Shallow Water Model

In addition to zonal mean circulation, eddy also plays an important role in shaping tropical climatology. Here we repeat the same scale analysis in  {ref}`scale_analysis`. 


```{math}
:label: shallow_water_1
\begin{align}
\frac{\partial u}{\partial t} + [u\frac{\partial u}{\partial x}+ v\frac{\partial u}{\partial y}+w\frac{\partial u}{\partial z}]-\beta yv & = -\frac{1}{\rho}\frac{\partial p}{\partial x}\\
\frac{\partial v}{\partial t} + [u\frac{\partial v}{\partial x}+ v\frac{\partial v}{\partial y}+w\frac{\partial v}{\partial z}]+\beta yu & = -\frac{1}{\rho}\frac{\partial p}{\partial y}\\
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z} &= 0 \\
\end{align}
```

Following {eq}`reduced_gravity2`, we can rewrite {eq}`shallow_water_1` 

```{math}
:label: shallow_water_1
\begin{align}
\frac{\partial u}{\partial t} + [u\frac{\partial u}{\partial x}+ v\frac{\partial u}{\partial y}]-\beta yv & = -g\frac{\partial H}{\partial x}\\
\frac{\partial v}{\partial t} + [u\frac{\partial v}{\partial x}+ v\frac{\partial v}{\partial y}]+\beta yu & = -g\frac{\partial H}{\partial x}\\
w(z) & = -z (\frac{\partial u}{\partial x}+\frac{\partial v}{\partial x}) + w_0 
\end{align}
```

where $w(z)$ is the 