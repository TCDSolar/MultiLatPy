# MultiLatPy
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/629036694.svg)](https://zenodo.org/badge/latestdoi/629036694)

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

Both require the same scaling used for the initial grid to compute the GDOP values, the resulting array of GDOPs and the satellite positions used. **NOTE: the plotting tool assumes that the grd is square, i.e. <code>nx = ny</code>**

<code>gdop2D_unc</code> also require inputs <code>prop_v</code> (propagation velocity of the source) and <code>delta_t</code> (time cadence of observations) to compute the total uncertainty at each grid point.

Both functions out the resulting figure and axes objects for further editing if required. ```gdop2D_unc``` also provides an array of the total uncertainty calulated at each point: ```zvals_r```.

<sub>**Authors:** Dale Weigt, Shane Maloney, Alberto Canizares, Sophie Murray, Peter Gallagher</sub>

## Bancroft multilateration
The Bancroft method uses linear algerba to obtain a direct solution of the source/reciever position (and clock offset), without the need for 
any priori knowledge of source/reciever location [1]. A quick and efficient algorithm to find a souce location from 4+ satellites using psuedorange matrix of satellites. The code can be adapted to account for fewer recievers. Mathematical derivation found here: https://gssc.esa.int/navipedia/index.php/ 

**NOTE: calculation here does not include ionospheric and tropospheric terms and assumes the emission is emitted isotropically. The location calulated can be made frequency dependent if you specify times at a certain frequency**

The Bancroft algorithm here is adaped from code provided by 10GeV: https://codereview.stackexchange.com/questions/253982/bancrofts-method-implementation and can be accesed using:

```shell
python bcroft.py
```

The function ```tdoa_ban``` is inisde this script and utilises the 4x4 Minkowski matrix <code>diag(1,1,1,-1)</code> and <a href="https://mathworld.wolfram.com/LorentzianInnerProduct.html">Lorentz Inner product</a> to calculate the positions.

**Inputs:**
<ul>
  <li> <code>stations_pos</code> : <b>array of N elements of shape (1, 3)</b>; where N >= 4 is the number of recievers (x,y,z) coordinates of satellite/reciever positions. Units of R_sun </li>.
  <li> <code>alpha</code>: <b>scalar</b>; scaling used to adapt speed of light <code>(2.99792458e8/sol_2_m)*alpha</code>. Speed in units of R_sun/s </li>
  <li> <code>tdoa</code>: <b>array of length (1,n) </b>; Array of determiend TDOAs to be used for Bancroft algorithm in units of seconds.
</ul>

**Outputs:**

<code>Vertexer</code> class used to calculate Bancroft location of source. Creates a class once the matrix equation has been solved. 'Vertexer' class first used to find:
<ol> 
  <li>Solve Bancroft formulation using N satelite geometry of the form  (for N >= 4):
                
   <code> [variable] = Vertexer(np.array([[x_1, y_1, z_1],[x_2, y_2, z_2],[x_3, y_3,z_3], [x_4, y_4,z_4],...,[x_N, y_N, z_N]]))</code>
                  </li>
  <li>Inset given TDOA to establish source location:
    
    [variable].find(np.array([[t_1], [t_2], [t_3], [t_4],...,[t_N]]))
  </li>
</ol>
Using these, the final output is the location of the source. Algorithm produces 2 solutions: correct/most optimum solution ia taken to be one of the smallest Residiual sum of squares (RSS) error.

The output is therefore the position of the source in (x,y,z) coordinates in the system used for the satellite positions. In units of Rsun.

The remaining 'bcroft' files consisit of unit tests to evaluate the accuracy of ```tdoa_ban```:

<ol>
  <li><code>bcroft_static_test.py</code>: evaluates accuracy of algorithm using various source/satellite positions. No noise is added to TDOAs </li>
  <li><code>bcroft_sc_noise_test.py</code> and <code>bcroft_source_noise_test.py</code>: identical to previous test must adds Gaussian noise of varying degrees to the calculated TDOAs
 </ol>

All tests use <code>pytest</code> module to run.

<sub>**Authors:** Dale Weigt, Shane Maloney, Alberto Canizares, Sophie Murray, Peter Gallagher</sub>


## Requirements
<ul>
  <li><b>python 3.9.13</b></li>
  <li>dataclasses</li>
  <li>matplotlib 3.5.2</li>
  <li>numpy 1.23.1</li>
  <li>pandas 1.4.3</li> 
  <li>pytest 7.1.2</li>
</ul>

### References
<ul>
  <li> [1] Bancroft, S. (1985), <i> An Algebraic Solution of the GPS Equations </i>, IEEE Transactions on Aerospace Electronic Systems, <b>21</b>, 1, 56-59, doi: <a href="https://ieeexplore.ieee.org/document/4104017">10.1109/TAES.1985.310538</a> </li>
</ul>

