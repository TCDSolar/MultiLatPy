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

<sub>**Authors:** Dale Weigt, Shane Maloney, Alberto Canizares, Sophie Murray, Peter Gallagher</sub>



## Requirements


