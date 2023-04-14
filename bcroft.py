#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 15:08:15 2023

@author: daleweigt
"""

#M = np.diag([1,1,-1])

def tdoa_ban(stations_pos, alpha, tdoa):
    """
    

    Parameters
    ----------
    stations_pos : array of N elements of shape (1, 3), where N >= 4 is the number of recievers
        (x,y,z) coordinates of satellite/reciever positions. Units of R_sun
    alpha: scalar
        scaling used to adapt speed of light (2.99792458e8/sol_2_m)*alpha. Speed in units of 
        R_sun/s
    M: diagonal matrix (4 x 4), of form diag([1,1,1,-1])
        Minkowski matrix required for Lorentz inner product. [Already in function]
    tdoa: array of length (1,n)
        Array of determiend TDOAs to be used for Bancroft algorithm

    Returns
    -------
    'Vertexer' class used to calculate Bancroft location of source
        Creates a class once the matrix equation has been solved. 'Vertexer' class 
        first used to find:
            1) Solve Bancroft formulation using N satelite geometry of the form 
              (for N >= 4):
                  <variable> = Vertexer(np.array([[x_1, y_1, z_1],[x_2, y_2, z_2],[x_3, y_3,z_3], [x_4, y_4,z_4],...,[x_N, y_N, z_N]])) 
        
            2) Inset given TDOA to establish source location:
                  <variable>.find(np.array([[t_1], [t_2], [t_3], [t_4],...,[t_N]]))
                  
    Using these, the final output is the location of the source.
    Algorithm produces 2 solutions - most optimum solution ia taken to be one of the smallest RSS error
    
    
    Motivation
    ----------
    - Obtain a direct solution of the source/reciever position (and clock offset), without the need for 
    any priori knowledge of source/reciever location
    - Quick look way to find souce location from 4+ satellites using psuedorange matrix of satellites
    - Adpated from code by 10GeV: https://codereview.stackexchange.com/questions/253982/bancrofts-method-implementation
    - Mathematical derivation found here: https://gssc.esa.int/navipedia/index.php/Bancroft_Method#:~:text=The%20Bancroft%20method%20allows%20obtaining,knowledge%20for%20the%20receiver%20location.
    - NOTE: ionospheirc and tropospheric terms not included
    """
    
    
    from dataclasses import dataclass
    import numpy as np
    from numpy.linalg import inv, norm

    sol_2_m = 6.957E8
    M = np.diag([1, 1, 1, -1]) # Minkowski matrix

    @dataclass
    class Vertexer:
    
        nodes: np.ndarray
        
        v = (2.99792458e8/sol_2_m)*alpha
# speed of light
        # Defaults
    
    
        def errFunc(self, point, times):
            
            """
            Parameters
            ----------
            point: array of shape (n, 3), where n >= 4 is the number of recievers
                (x,y,z) coordinates of satellite/reciever positions.
            times:
                Array of determiend TDOAs to be used for Bancroft algorithm

            Returns
            -------
            Error function giving the sum of squared residuals (RSS) of the TOA of predicted point at each satellite 
            - algorithm takes sum of squared residuals (RSS error) of the time of arrival of the predicted point at 
            each node, minus the given time of arrival.
            """
            # Return RSS error
            error = 0
    
            for n, t in zip(self.nodes, times):
                error += ((np.linalg.norm(n - point) / self.v) - t)**2
    
            return error
    
        def find(self, times):
            def lorentzInner(v, w):
                """
                Parameters
                ----------
                v, w (vector/matrix): v must have shape (m x 4), n must have shape (4 x n)
                M: diagonal matrix (4 x 4), of form diag([1,1,1,-1])
                
                Returns
                -------
                Output of Lorentz inner product deifned as:
                    
                    <v, w> = v.T @ M @ w
                                                                        [w1]                              
                                                                        |w2|
                          = [v1, v2, v3, v4,..., vN] @ diag(1,1,1,-1) @ |w3|
                                                                        [w4]
                
                """
                # Return Lorentzian Inner-Product
                return np.sum(v * (w @ M), axis = -1)
    
            A = np.append(self.nodes, times * self.v, axis = 1) # generates psudeorange matrix for N satellites
            
            # Express linear equation in a form of a quadratic...
            b = 0.5 * lorentzInner(A, A) # ... column vector of inner product of satellite positions and psudeoranges... 
            oneA = np.linalg.solve(np.dot(A.T, A), np.dot(A.T,np.ones(len(stations_pos)))) # setting up matrix algebra
            invA = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, b))
    
            # Solve for lamba in bancroft equation...
            solution = []
            for Lambda in np.roots([ lorentzInner(oneA, oneA),          
                                    (lorentzInner(oneA, invA) - 1) * 2, 
                                      lorentzInner(invA, invA),          
                                    ]):
                
    
                X, Y, Z, T = M @ np.linalg.solve((A.T @ A), (A.T @ (Lambda * np.ones(len(stations_pos)) + b)))#, rcond=None)
                solution.append(np.array([X,Y,Z]))
                #... provides two possible positional solutions...
            
            return min(solution, key = lambda err: self.errFunc(err, times))
            #return solution - takes argument with the smallest RMS error
            
    # Genertes arrays of:
        #(1) matrix positional coordinates in (x,y,z) heliocentric earth ecliptic (HEE) coordinates
        #(2) corresponding TDOA foudn for each satellite
    
    sat_matrix = []
    tdoa_matrix = []
    for num in range(len(stations_pos)):
        sat_matrix.append([stations_pos[num][0], stations_pos[num][1], 
                          stations_pos[num][2]])
        
    for tt in tdoa:
        tdoa_matrix.append([tt])    
            
    myVertexer = Vertexer(np.array(sat_matrix)) # Generates positions in the bcroft class...
    bcroft_loc = np.real(myVertexer.find(np.array(tdoa_matrix))) #... and finds positions based on TDOAs.

    return bcroft_loc