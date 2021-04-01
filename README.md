# Langevin Equation in Organelle Transport models 

This repository provides a numerical approach to study the effect of a molecular motor in organelle transport. 

## Models
We consider two different models in organelle transport. The first set of models contemplates the most fundamental models of organelle transport:

- In viscous media only (Langevin equation),
- In viscoelastic media (Modified Langevin equation).

For the second set, we consider organelle transport mediated by molecular motors. This set includes:

- Standard Model Langevin + Monte Carlo ([SLMC](SLMC/SLMC.ipynb)),
- Standard Model Langevin + Monte Carlo with continous step ([SLMC-CS](SLMC_CS/SLMC_CS.ipynb)),
- Time-averaged Force (TAF).
- Two-dimensional model (2DM)

The last model here included is the extension in two dimensions of the SLMC-CS model. We include an additional folder with the model of a molecular motor with no cargo.


Each model includes a jupyter notebook or tex file with the respective theoretical model explained.
