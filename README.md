The SIR model is a compartmental model of infectious disease progression (e.g., https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model). This repository contains some very simple SIR models I tried to implement. This was done mostly out of curiosity during a COVID lockdown.

_sir.c_

This program just solves the very simple numerical differential equations for a given starting point.

_sird.c_

This implements a diffusion version of the SIR model. A 2D grid is created, and certain points 'infected'. A Lorentzian map is used to enable the infection to diffuse across the surface.

I've also added some quite nice (in my opinion!) videos generated from these.
