# RESolve/OPTIMIZE

Suite of Matlab codes codes to fit a Richards Equation model to data from Johnstown Castle Co. Wexford

# Project overview:

* Data digitized from the paper of [Diamond and Sills](https://t-stor.teagasc.ie/entities/publication/c5d35858-c39a-44c5-84ba-942b8b9f2b3f)
* Nonlinear least-squres fitting with K_sat and alpha as parameters.
* Procedure descriibed in [notes2](https://github.com/onaraighl/RESolve/blob/main/OPTIMIZE/notes2.pdf).

# Detailed code structure:

Optimization procedure is called from function `myOptimization.m`.  Input arguments are null and the output arguments are `p_opt` and `fval`:

* `p_opt` contains the optimum values of K_sat and alpha;
* `fval` is the corresponding value of the cost function.

The cost function is computed by comparing the model and the observations across two years, from January 1st 1998 to Decemmber 31st 1999.  The cost function is the sum of the squares of differences, divided by the number of days in this interval.  The observations are loaded from a `.csv' file for the purpose of computing the cost function.

The code uses surrogate optimization, which is a good technique when each evaulation of the cost function requires the solution of a PDE.  Furthermore, the code is run in parallel, meaning that multiple evaluations of the cost function can be done simultaneously.

The cost function calls `RESolve`, through the Matlab function `richards_pdepe.m`.  The input arguments are stored in a vector `inputParams`, whose entries are Ksat and alpha.  The solution of the PDE requires the calling of various subfunctions, as documented elsewhere, including:

* `get_EP.m`
* `get_rainfall.m`
* `mathias_model.m`
* `vg_inv_Se.m`
* `vg_model.m`





