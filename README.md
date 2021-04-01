# Langevin Equation in Organelle Transport models 

This repository provides a numerical approach to study the effect of a molecular motor in organelle transport. 

## Models
We consider two different sets of models. The first set of models (1_no_motors/)contemplates the most fundamental form of transport, driven by thermal fluctuations and cellular crowding, but in the absence of molecular motors:

- In viscous media only ([Difussion](1_no_motors/1-1_viscous/particle-difussion.f)).
- In viscoelastic media ([Difussion with virtual particles](1_no_motors/1-2_viscoelastic/particle_difussion+virtual_particles.f)).


For the second set, we consider organelle transport mediated by molecular motors. This set includes:

- Standard Model Langevin + Monte Carlo ([SLMC](2_motors/SLMC/SLMC-pv.f)),
- Standard Model Langevin + Monte Carlo with external force ([SLMC-Fext](2_motors/SLMC_Fext/SLMC_Fext.f)),

As a third additional set, we include the model of a molecular motor with no cargo .


Each model includes a jupyter notebook or tex file with the respective theoretical model explained.

## Folders

Each folder specifies the type  regime of the particle, and includes the initial conditions (.ini file). As an example, the directory tree for the first set of models with no motors is shown below

```
├── 1_no_motors
│   ├── 1-1_viscous
│   │   ├── cond_diff.ini
│   │   └── particle-difussion.f
│   ├── 1-2_viscoelastic
│   │   ├── cond_vp.ini
│   │   └── particle_difussion+virtual_particles.f
│   ├── no_motors_theory.pdf
│   └── no_motors_theory.tex
```

## References

I. Goychuck. Coexistence and eficiency of normal and anomalous transport by molecular motors in living cells. _Phys. Rev._, 80:1-74, 2009.

Bouzat S. Infuence of molecular motors on the motion of particle in viscoelastic media. _Phys. Rev. E_, 89:062707, 2014. 8, 10