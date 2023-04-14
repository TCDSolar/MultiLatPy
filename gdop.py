#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 11:04:35 2023

@author: daleweigt
"""

def gdop2D(scale, stations_pos, nx=1000, ny=1000):
    
    """

    Parameters
    ----------
    scale: sets up the size of grid needed in AU (e.g. 1.5 = 1.5 AU in x- and y-dir)
    
    stations_pos : array of Ne elements of shape (1, 2), where N >= 2 is the number of recievers
        (x,y,z) coordinates of satellite/reciever positions.
        
    nx: default is 1000 - number of grid points along x-axis 
    
    ny: default is 1000 - number of grid points along y-axis


    Returns
    -------
    sat_err: array of calculated GDOPs found from the covariance matrix of normalised psuedoranges of each satellite
            (see below for more details)
            
    tim_err: array of calculated TDOPs found from the covariance matrix - e.g., works out uncertainty of how the clocks of each
            satellite are synchronised.
            
    cov_Q: list of covariance matrices for each position. Provides info on the variances in each coordinate and how they are
            related.
                  
    
    Motivation
    ----------
    - Geometric Dilation of Precision (GDOP) is a useful measure of how accurate multilateration (or source tracking) is at a 
        given position on a grid.
    - Widely used in GPS measurements.
    - Values > 20 are considered non-accurate and should be discarded.
    - Here, we apply this to a constellation spacecraft setting in a more simplified manner.
    - The aim is to determine which positions in the inner heliosphere (e.g., within Earth's orbit) would provide the best coverage
        to track various phenomena (e.g., solar radio burts).
    - Useful to be used in tandem with your favourite multilateration technique.
    - Purely algerbraic - uses the covariance matrix of the normalised psudeoranges of th spacecraft (Q = (A.T @ A)^-1)
    - GDOP is the trace of the covariance matrix 
    - TDOP is the final element along the diagonal
    - See slides for further details...
    """
    
    
    import numpy as np
    from numpy.linalg import inv, norm

    
    AU_2_sol = 215.032 # AU -> solar radii
    
    #setup grid for GDOP calculation

    x = np.linspace(-1*scale*AU_2_sol, scale*AU_2_sol, nx+1) # in AU
    y = np.linspace(-1*scale*AU_2_sol, scale*AU_2_sol, ny+1)
    xgrid, ygrid = np.meshgrid(x, y)
    # take each point on the grid as a source
    xdata, ydata = np.ravel(xgrid), np.ravel(ygrid)
    pnts = np.vstack((xdata,ydata)).T
    
    """Satellite parameters and GDOP calc"""
    
    # (x,y) coordinate of each spacecraft
    satx, saty = stations_pos[:,0], stations_pos[:,1]
    
    sat_err = []
    tim_err = []
    cov_Q = []

    for i  in range(len(xdata)):
        # caluclate distance between each simulated source and spacecraft...
        R_dist = np.sqrt((satx - xdata[i])**2 + (saty - ydata[i])**2)
        
        matA = []
        for num in range(len(R_dist)):
            matA.append([(satx[num] - xdata[i])/R_dist[num],
                        (saty[num] - ydata[i])/ R_dist[num], 1])
        # create  psuedorange matrix using differend in distances (and normalised)  
        matA = np.asarray(matA)
        
        covA = inv(np.dot(matA.T, matA)) # determine covariance matrix
        gdop = np.sqrt(np.trace(covA))# calculate GDOP...
        tdop = np.sqrt(np.diag(covA)[len(covA)-1]) # calculate TDOP
        sat_err.append(gdop)
        cov_Q.append(covA)
        tim_err.append(tdop)

    sat_err = np.array(sat_err)
    tim_err = np.array(tim_err)
    
    return sat_err, tim_err, cov_Q 
    
    
