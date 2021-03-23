

## Organelle transport models - Theory
### Cargo motion (no molecular motors)
The initial model here considered is the one-dimensional generalized Langevin equation (GLE), which describes
subdiffusive particle motion when in a viscoelastic medium in equilibrium conditions:


```math
mx(t) = - \int_0^t \gamma(t − t')x(t')dt'+\xi(t) ,
```

where, $x(t)$ is the position of the particle, $m$ is its mass, and $\gamma(t)$ is the frictional kernel given by

\gamma(t) = \frac{\gamma_0}{(1 − \alpha)^t}t^{\alpha},

where $\xi(t)$ is the time-correlated Gaussian thermal noise that satisfies the fluctuation-dissipation relation

- <img src="\xi(t)\xi(t') = k_{BT}\gamma(t − t')." />

GLE coupled to a stochastic stepping dynamics is used to model the motion of particles driven by molecular motors in ithe cytoplasm, here represented by a viscoelastic medium.

In the absence of motors, the model produces subdiffusive motion of particles characterized by a power-law scaling of the mean square displacement versus the lag time as $t^\alpha$, with $0<\alpha<1$, similar to the one observed in cells.

#### Cargo transport meditated by molecular motors
In active transport, it is usual to study the cargo response when external forces are applied. 
In experimental settings, due to either hydrodynamic friction of the media, or either the effect of optical traps, external forces are applied to the organelle.For the latter, the cargo displacement is calculated when it is confined in an armonic potential. Under this consideration, a fixed value of the external force is added to the initial equation of movement. Hence the equation for organelle transport mediated by molecular motors is given by:


\eta_0\dot{x}(t)=-\sum_{i=0}^{N-1}k_i(x(t)-x_i(t))+\xi_0(t)+F_m(x(t),x_m(t))-L_{ext},


where $L_{ext}$ is the external force due to the optical trap, 0$\leq L_{ext}\leq7$pN.




 
