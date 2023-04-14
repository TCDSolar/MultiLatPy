#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 11:19:30 2023

@author: daleweigt
"""
from bcroft import tdoa_ban 
import numpy as np
from astropy.constants import c, m_e, R_sun, e, eps0, au

# content of test_expectation.py
import pytest


@pytest.mark.parametrize(
    "source",
    [([0,0,1E-3]),([10,0,1E-3]),([0,10,1E-3]),([10,10,1E-3])]
)
@pytest.mark.parametrize(
    "NoiseLevel",
    [0,0.5,1,2]
)
def test_source_pos(source,NoiseLevel):
    
    """
    

    Parameters
    ----------
    source : array of source locations in (x,y,z) HEE coordinates. units of Rsun
        position vectors of selected spacecraft positions to run tests over

    NoiseLevel : float
        float to adjust noise level added to data (in this case, the times used)

    Returns
    -------
    - Evaluation of if simulated bancroft values are within the threshold tolerance of true source
    values.
    - relative tolerance (rtol) = 1E-7 (default)
    - absolute tolerance (atol) = 2 Rsun (default) -> both can be changed
    - Tolerance equation for bancroft positions (x_bc) is:
        
        (x_bc * rtol) + atol 
        
        where x_bc passes unit test if it is within this threshold of the true value (e.g., here ~2 Rsun,
        for atol = 2)
    - Allows to assess how noisy data affects the bancroft algorithm. 

    """
    
    stations_rsun = np.array([[-26, -165, 10], [110, -174,-10], [-110, 174,10], [215, 0,-10]])
    scale = 1/6.957E8 # m/s to R_sun/s
    # TEST SOURCE
    xx, yy, zz = source  # n = 2
    x_true = np.array([xx, yy, zz])  # true source position (units: R_sun)
    v_true = c.value * scale # speed of light (m/s)
    d_true = np.linalg.norm(stations_rsun - x_true, axis=1)
    t1_true = d_true / v_true #- 10*60  # true time of flight values

    t_obs = t1_true + np.abs(t1_true - np.random.normal(t1_true, NoiseLevel, *t1_true.shape))# noisy observations
    
    tdoa = t_obs - np.min(t_obs) # TDOAs calculated from true source
    x_bc = tdoa_ban(stations_pos=stations_rsun, scale=scale, tdoa=tdoa) # calculated bancroft location
    np.testing.assert_allclose(x_true, x_bc,atol=2) # compares the values to absolute tolerance
