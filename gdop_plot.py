#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 11:15:09 2023

@author: daleweigt
"""

def gdop2D_plotter(scale, gdop, stations_pos):
    
    """

    Parameters
    ----------
    scale: sets up the size of grid needed in AU (e.g. 1.5 = 1.5 AU in x- and y-dir)
           NOTE: scale should be the same used for GDOP calculation
           
    gdop: array of calculated GDOPs from stations_pos array
        
    stations_pos : array of satellites position used in GDOP calculation. Units of solar
            radius


    Returns
    -------
    gdop_fig, ax: figure and axes used for GDOP plot.
    
    Motivation
    ----------
    - To plot GDOP from given statellite positions used in gdop2D function
    """
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt 
    from matplotlib.colors import LogNorm
    import matplotlib.patheffects as pe

    AU_2_sol = 215.032 # AU -> solar radii
    
    # Setting up same grid as used in GDOP calculation
    
    x = np.linspace(-1*scale*AU_2_sol, scale*AU_2_sol, np.sqrt(len(gdop)).astype(int)) # in AU
    y = np.linspace(-1*scale*AU_2_sol, scale*AU_2_sol, np.sqrt(len(gdop)).astype(int))
    xgrid, ygrid = np.meshgrid(x, y)
    xdata, ydata = np.ravel(xgrid), np.ravel(ygrid)

    # Coordinates of Sun
    sunx, suny, rsun = 0.0, 0.0, 1.0 

    # Prepping source coordinates, GDOP and TDOP data
    df = pd.DataFrame(dict(x=xdata, y=ydata, err=gdop))
    xcol, ycol, zcol = "x", "y", "err"
    df = df.sort_values(by=[xcol, ycol])
    xvals = df[xcol].unique()
    yvals = df[ycol].unique()
    zvals = df[zcol].values.reshape(len(xvals), len(yvals)).T
    
    
    # Generate plot of earth's orbit and Sun's position
    earth_orb = plt.Circle((sunx, suny), 1.0*AU_2_sol, color='k', linestyle='--', fill=False,
                           label='Earth orbit')
    sun = plt.Circle((sunx, suny), rsun, color='orange')#, label='Sun (0,0)')
    
    # Set up figure for plotting GDOP    
    gdop_fig, ax = plt.subplots(figsize=(10,10)) 
    ax.axhline(0, alpha=.1, color='white')
    ax.axvline(0, alpha=.1, color='white')
    levels = np.logspace(0,3,101)
    
    # Define color bar limits 
    vmax = levels[-1]
    gplot = plt.contourf(xvals, yvals, zvals, levels=levels, cmap='gnuplot_r',norm=LogNorm(vmax=vmax)) 
    
    # Plot Earth orbit and Sun's position
    ax.add_artist(earth_orb)
    ax.add_artist(sun)
    
    # Plot position and labels for each satellite used in calculation...
    for jj in range(0, len(stations_pos)):
        ax.scatter(stations_pos[jj][0], stations_pos[jj][1], marker='^', color='white',
                   edgecolor='k', s=100, label='Receiver' if jj == 0 else "")
    
    
        ax.text(stations_pos[jj][0]-0.15*AU_2_sol, stations_pos[jj][1],'R{}'.format(jj+1), size=14,
                color='white', path_effects=[pe.withStroke(linewidth=2, foreground="black")])
    
    
    # and Sun label                   
    ax.text(sunx+0.02*AU_2_sol, suny-0.06*AU_2_sol,'Sun', size=14, color='white',
            path_effects=[pe.withStroke(linewidth=2, foreground="black")])
    
    ax.set_xlabel('x (HEE: R$_\odot$)',size=14)
    ax.set_ylabel('y (HEE: R$_\odot$)', size=14)
    ax.tick_params(labelsize=14)
    plt.legend(fontsize = 10, framealpha=1.0)
    cbar = plt.colorbar(gplot, fraction = 0.03, pad=0.04)
    
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('GDOP', size=16)

    return gdop_fig, ax



def gdop2D_unc(scale, gdop, stations_pos, delta_t, prop_v):
    
    """

    Parameters
    ----------
    scale: sets up the size of grid needed in AU (e.g. 1.5 = 1.5 AU in x- and y-dir)
           NOTE: scale should be the same used for GDOP calculation
           
    gdop: array of calculated GDOPs from stations_pos array
    
    stations_pos : array of satellites position used in GDOP calculation. Units of solar
          radius.
    
    delta_t: time cadence of observations - generate propagtion speeds. units of seconds
    
    prop_v: propagtion velocity of source. units of m/s

    Returns
    -------
    unc_fig, ax: figure and axes used for GDOP plot.
    
    zvals_r: uncertainty in position due to GDOP measurements. units of solar radius
    Motivation
    ----------
    - To plot GDOP in terms of radial uncertianty from given statellite positions used in gdop2D function
    """
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt 
    from matplotlib.colors import LogNorm
    import matplotlib.patheffects as pe
    
    
    AU_2_sol = 215.032 # AU -> solar radii
    sol_2_m = 6.957E8 # Solar radius -> meteres

    c = prop_v/sol_2_m # convert propagation speed from m/s -> solar radii/s
    # Setting up same grid as used in GDOP calculation
    
    x = np.linspace(-1*scale*AU_2_sol, scale*AU_2_sol, np.sqrt(len(gdop)).astype(int)) # in AU
    y = np.linspace(-1*scale*AU_2_sol, scale*AU_2_sol, np.sqrt(len(gdop)).astype(int))
    xgrid, ygrid = np.meshgrid(x, y)
    xdata, ydata = np.ravel(xgrid), np.ravel(ygrid)

    # Coordinates of Sun
    sunx, suny, rsun = 0.0, 0.0, 1.0 

    # Prepping source coordinates, GDOP and TDOP data
    df = pd.DataFrame(dict(x=xdata, y=ydata, err=gdop))
    xcol, ycol, zcol = "x", "y", "err"
    df = df.sort_values(by=[xcol, ycol])
    xvals = df[xcol].unique()
    yvals = df[ycol].unique()
    zvals = df[zcol].values.reshape(len(xvals), len(yvals)).T
    
    
    # Generate plot of earth's orbit and Sun's position
    earth_orb = plt.Circle((sunx, suny), 1.0*AU_2_sol, color='k', linestyle='--', fill=False,
                           label='Earth orbit')
    sun = plt.Circle((sunx, suny), rsun, color='orange')#, label='Sun (0,0)')
    
    # Plotting GDOP -> positional error
    zvals_r = np.sqrt(((((3*zvals)/(4*np.pi))**1/3)*(c*delta_t/2))**2 + (c*delta_t/2)**2) # assuming error is radial
    sun = plt.Circle((sunx, suny), rsun, color='orange')#, label='Sun (0,0)')
    earth_orb = plt.Circle((sunx, suny), 1.0*AU_2_sol, color='k', linestyle='--', fill=False)

    unc_fig, ax = plt.subplots(figsize=(10,10)) 
    ax.axhline(0, alpha=.1, color='white')
    ax.axvline(0, alpha=.1, color='white')
    levels = np.logspace(np.log10(np.amin(zvals_r)),2,101)
    
    # Define color bar limits 
    vmax = levels[-1]
    gplot = plt.contourf(xvals, yvals, zvals_r, levels=levels, cmap='gnuplot_r',norm=LogNorm(vmax=vmax)) 
    ax.add_artist(sun)
    ax.add_artist(earth_orb)

    # Plot position and labels for each satellite used in calculation...
    for jj in range(0, len(stations_pos)):
        ax.scatter(stations_pos[jj][0], stations_pos[jj][1], marker='^', color='white',
                    edgecolor='k', s=100, 
                    label='Receiver' if jj == 0 else "")

        ax.text(stations_pos[jj][0]-0.15*AU_2_sol, stations_pos[jj][1],'R{}'.format(jj+1), size=14,
                color='white', path_effects=[pe.withStroke(linewidth=2, foreground="black")],)


    ax.text(sunx+0.02*AU_2_sol, suny-0.06*AU_2_sol,'Sun', size=14, color='white',
            path_effects=[pe.withStroke(linewidth=2, foreground="black")])

    ax.set_xlabel('x (HEE: R$_\odot$)',size=14)
    ax.set_ylabel('y (HEE: R$_\odot$)', size=14)
    ax.tick_params(labelsize=14)
    plt.legend(fontsize = 10, framealpha=1.0)#, loc="upper left")
    cbar = plt.colorbar(gplot, fraction = 0.03, pad=0.04)

    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('Radial uncertainty (GDOP + c$\delta$t: R$_\odot$)', size=16)

    return unc_fig, ax, zvals_r
