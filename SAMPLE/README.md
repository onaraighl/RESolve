# RESolve/SAMPLE

**RESolve** is a suite of Matlab codes to solve the Richards Equation in 1D.  `RESolve/SAMPLE` is a sample implementation equipped with synthetic meteorological inputs.

# Code Structure:

The model is written in a single Matlab function, `richards_pdepe` which is contained in a `.m` file of the same name.    The main function `richards_pdepe` includes instructions for setting up the spatial grid, the time vector, and the physical parameters.  The main function also handles calls to Matlab's built-in PDE solver, `pdepe`.  The main function also contains important subfunctions:

* `richards_pde` - a function to define the actual PDE being solved, including definitions of the capacity c, the flux f, and the source s;
* `bc_fun` - a function to set up the boundary conditions
* `initial_h` - a function to define the initial conditions.

This structure is deliberate, as the subfunctions do utilize some global variables defined in `richards_pdepe`.

Additional important `.m` files are provided in this directory.  These are subfunctions which are called by `richards_pdepe`, but in a modular fashion.  These include:

* `get_EP.m' - a function to compute the potential evapotranspiration at a given time t.  All units are SI.
* `get_rainfall.m' - a function to compute the rainfall at a given time t.  All units are SI.
* `mathias_model.m` - a function to compute the root distribution function and the plant stress function.
* `vg_model.m` - a function to compute the Van Genuchten model parameters as a function of pressure head h.  Outputs are volumetric water content θ(h), capacity C(h), and relative permeability K(h).

# Parameters:

Van Genuchten parameters are defined in `richards_pdepe`, along with the soil column depth.  Synthetic meteorological data are also provided in `richards_pdepe`.  Parameters relating to the root distribution function and the plant stress function are hard-coded in `mathiasa_model.m`.




