# Electrostatic-Optics
Compute the isopotentials for a planar Einzel lens and trace charged particles through the lens as a simulated beam. Extract current density profile from image plane.
## Purpose
Simulate the performance of a simple electrostatic charged particle optical element.
Provides an alternative method to geometric optical formulations to characterize the focusing power of an Einzel lens.
## Applications
The mass and the charge of the particle can be changed to account for electrons and multiple ion species.
The focusing power can be adjusted by changing the applied potential.
The lens can function in accel- or decel- mode by controlling the sign of the potential on the center element.
### Sample Application
Find the voltages which cause crossovers in the beam at certain distances from the lens element and find the closest crossover in accel- mode without exceeding breakdown voltage of the elements.   
One way to do this is to script the Matlab code with the electrode voltage as a parameter and find crossover points and density functions to determine the current through an 'x' size virtual aperture located at the right edge of the domain.
