#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 11:21:21 2023

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

def test_source_pos(source):
    """
    
    Parameters
    ----------
    source : array of source positions in (x,y,z) HEE coordinates. units of Rsun
        position vectors of selected source positions to run tests over, without applying 
        Gaussian noise

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

    """
    stations_rsun = np.array([[-26, -165, 10], [110, -174,-10], [-110, 174,10], [215, 0,-10]])
    scale = 1/6.957E8 # m/s to R_sun/s
    # TEST SOURCE
    xx, yy, zz = source  
    x_true = np.array([xx, yy, zz])  # true source position (units: R_sun)
    v_true = c.value * scale # speed of light (m/s)
    d_true = np.linalg.norm(stations_rsun - x_true, axis=1)
    t_obs = d_true / v_true #- 10*60  # true time of flight values

    tdoa = t_obs - np.min(t_obs) # TDOAs calculated from true source

    x_bc = tdoa_ban(stations_pos=stations_rsun, scale=scale, tdoa=tdoa) # calculated bancroft location

    np.testing.assert_allclose(x_true, x_bc,atol=2) # compares the values to absolute tolerance

@pytest.mark.parametrize(
    "sc",
    [(np.array([[-200, 200,0.001], [200, 200,-0.001], [200, -200, 0.001], [-200,-200, -0.001]])),
     (np.array([[-200, 0, 0.001], [0, 200, 0.001], [200, 0, 0.001], [0,-200,0.001]])),
     (np.array([[-26, -165, 10], [110, -174,-10], [-110, 174,10], [215, 0,-10]])),
     (np.array([[300, 0, 5], [256, 300,-5], [215, -230, 5], [216, 100,-5]]))]
)

def test_sc_pos(sc):
    
    """
    
    Parameters
    ----------
    sc : array of spacecraft in (x,y,z) HEE coordinates. units of Rsun
        position vectors of selected spacecraft positions to run tests over, without applying 
        Gaussian noise

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

    """
    
    stations = sc
    scale = 1/6.957E8 # m/s to R_sun/s
    # TEST SOURCE
    xx, yy, zz = [25,31,0]  # n = 2
    # xx, yy = -7.3, 8.1        #  n = 2
    x_true = np.array([xx , yy , zz])  # true source position (m)
    v_true = c.value * scale  # speed of light (m/s)
    d_true = np.linalg.norm(stations - x_true, axis=1)
    t_obs = d_true / v_true #- 10*60  # true time of flight values
    t_obs = d_true / v_true #- 10*60  # true time of flight values

    tdoa = t_obs - np.min(t_obs) # TDOAs calculated from true source

    x_bc = tdoa_ban(stations_pos=stations, scale=scale, tdoa=tdoa) # calculated bancroft location

    np.testing.assert_allclose(x_true, x_bc,atol=2) # compares the values to absolute tolerance
