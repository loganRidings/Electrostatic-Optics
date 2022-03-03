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
<p align="center">
<img src=https://github.com/loganRidings/Electrostatic-Optics/blob/main/sampleBeam.png alt="Beam Traces" width="400"/>  
<img src=https://github.com/loganRidings/Electrostatic-Optics/blob/main/currentDensityPlot.png alt="Beam Traces" width="400"/>
</p>

####  Required Matlab Add-ons
Parallel and Distributed Computing Toolbox (optional)   
Parial Differential Equation Toolbox (required)

####  Running without Parallel Toolbox
Just change the `parfor` to a `for` loop.
I built it to easily switch back and forth, however, the runtime is at least 4x for sequential ODE solving...
