# RESolve/SAMPLE

**RESolve** is a suite of Matlab codes to solve the Richards Equation in 1D.  `RESolve/SAMPLE` is a sample implementation equipped with synthetic meteorological inputs.

# Code Structure:

The model is written in a single Matlab function, `richards_pdepe` which is contained in a `.m` file of the same name.    The main function `richards_pdepe` includes instructions for setting up the spatial grid, the time vector, and the physical parameters.  The main function also handles calls to Matlab's built-in PDE solver, `pdepe`.  The main function also contains important subfunctions:

* `richards_pde` - a function to define the actual PDE being solved, including definitions of the capacity c, the flux f, and the source s;
* `bc_fun` - a function to set up the boundary conditions
* `initial_h` - a function to define the initial conditions.

This structure is deliberate, as the subfunctions do utilize some global variables defined in `richards_pdepe`.

Additional important `.m` files are provided in this directory.  These are subfunctions which are called by `richards_pdepe`, but in a modular fashion.  These include:

* `get_EP.m` - a function to compute the potential evapotranspiration at a given time t.  All units are SI.
* `get_rainfall.m` - a function to compute the rainfall at a given time t.  All units are SI.
* `mathias_model.m` - a function to compute the root distribution function and the plant stress function.
* `vg_model.m` - a function to compute the Van Genuchten model parameters as a function of pressure head h.  Outputs are volumetric water content θ(h), capacity C(h), and relative permeability K(h).
* `vg_inv_Se.m` - a function to invert the Van Genuchten model.  So, for an input value θ, an output value of h is computed.

# Parameters:

Van Genuchten parameters are defined in `richards_pdepe`, along with the soil column depth.  Synthetic meteorological data are also provided in `richards_pdepe`.  Parameters relating to the root distribution function and the plant stress function are hard-coded in `mathiasa_model.m`.

# Running the code:

Currently, one types at the Matlab command prompt: 

`[z,t,h_tz,theta_tz,sol]=richards_pdepe();`

Here, the inputs are NULL and the outputs are:

*  z        - z-coordinate
*  t        - vector of time values
*  h_tz     - spacetime array of values of pressure head (negative for unstaturated).
*  theta_tz - spacetime array of values of theta, in units of mm^3/mm^3.
*  sol      - matlab output (included here as diagnostic).

Note that **SI units are used throughout**.  In `SAMPLE/` it is recommended to save the results as a structure in a `.m` file:

`temp=struct('t',t,'z',z,'h_tz',h_tz,'theta_tz',theta_tz)`
`save('data.mat','temp')`

Then, some postprocessing steps can be run, using `data.mat' as an input.

# Postprocessing

Two postprocessing functions are currently provided:

*  `postprocessTheta.m'
*  `postprocessPlot.m'

The function `postprocessTheta.m' takes **inputs**

*  t        - vector of time values
*  z        - z-coordinate
*  theta_tz - spacetime array of values of theta, in units of mm^3/mm^3.

and produces as an **output** a vector `THETA`, which is the storage water level, as a function of time, in units of meters.  The water storage level is typically denoted in the literature by (Θ, big-theta), and is computed from θ via the trapezoidal rule.

In addition, `postprocessPlot.m' is a function taking NULL inputs and producing NULL outputs which operates on the aforementioned `data.mat`.  The function loads `data.mat` and calls `postprocessTheta.m`, producing a pre-formatted plot of the storage water level (Θ, big-theta) over time.




