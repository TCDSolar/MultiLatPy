# MultiLatPy
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This package contains the Python tools/algorithms to perform and/or aid multilateration techniques of their choice. These methods are purely algerbraic and only require the spacecraft positions used for analysis (as long as each position is in the same corrdinate system) and the Time Differnce of Arrival (TODA) for each reciever. The code is set up to process and compute the estimated position of a source using the Bancroft (bcroft) Method - one of many different techniques and provide useful aids (such as the Geometric Dilution of Precision or GDOP) to elimate regions of high uncertainty which will likely provide erroneous results.

## Using the data

Both the GDOP and bancroft algorithms require satellite positions (x,y,z) in the same coordinate system. Default system is Heliocentric Earth Eclitpic (HEE). The positions should be in units of solar radii (Rsun).

TDOA for each spacecraft can be determined from timestamped datasets and have units of seconds.

## GDOP calculator
Uses the satellite positions to create a normalised psuedorange matrix for each reciever. From this, the user selects the ```scale``` of their GDOP grid of length ```nx``` x ```ny``` (defaulted to 1000 points) to evaluate the GDOP value at each position. 

```shell
python gdop.py
```

```gdop.py``` contains the ```gdop2D``` function which calualtes the GDOP in a 2D system. This algorithm uses matrix algerba to determine the covariance matrix of the pseudoranges. The GDOP is then the trace of this covariance matrix. Other values can be computed such as the Time Dilution of Precision (the final diagonal entry).


**Input:** ```scale``` (size of grid in units of AU); ```stations_pos``` (positions of N-satellites used in calculation, where N >= 2); ```nx``` (size of grid in x-direction; ```ny``` size of grid in y-direction.

**Output:** ```sat_err``` (array of GDOP values for each position on grid); ```tim_err``` (array of TDOP calues for each position on grid); ```cov_Q``` (covariance matrix used for each posiiton, may be of use for further analysis)

Note: GDOP values > 20 are considered to be regions of high uncertainty and should be eliminated from further analysis

The GDOP values can be plotted onto the grid using:

```shell
python gdop_plot.py
```

This script contains two functions for plotting:

<ol>
  <li> <code>gdop2D_plotter</code>: plot GDOP unitless values from given satellite positions used in <code>gdop2D</code> function </li>
  <li> <code>gdop2D_unc</code>: plot errors in location in terms of total radial uncertianty from GDOP (assuming spherical symmetry) and time cadence of the observations given statellite positions used in <code>gdop2D</code> function </li>   
</ol>

<sub>**Authors:** Dale Weigt, Shane Maloney, Alberto Canizares, Sophie Murray, Peter Gallagher</sub>

Both require the same scaling used for the initial grid to compute the GDOP values, the resulting array of GDOPs and the satellite positions used. **NOTE: the plotting tool assumes that the grd is square, i.e. <code>nx = ny</code>**

<code>gdop2D_unc</code> also require inputs <code>prop_v</code> (propagation velocity of the source) and <code>delta_t</code> (time cadence of observations) to compute the total uncertainty at each grid point.

Both functions out the resulting figure and axes objects for further editing if required. ```gdop2D_unc``` also provides an array of the total uncertainty calulated at each point: ```zvals_r```.

## Requirements

### References
<ul>
  <li></li>
</ul>

