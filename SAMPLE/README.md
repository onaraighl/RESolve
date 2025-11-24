# RESolve/SAMPLE

**RESolve** is a suite of Matlab codes to solve the Richards Equation in 1D.  `RESolve/SAMPLE` is a sample implementation equipped with synthetic meteorological inputs.

# Code Structure:

The model is written in a single Matlab function richards pdepe.m, which includes a main function setting up the spatial grid, time vector, soil hydraulic functions, rainfall input,
and calls to pdepe. The main function pdepe.m contains subfunctions:

* `richards_pde` - a function to define the actual PDE being solved, including definitions of the capacity c, the flux f, and the source s;
* bc_fun - a function to set up the boundary conditions
* initial_h - a function to define the initial conditions.






