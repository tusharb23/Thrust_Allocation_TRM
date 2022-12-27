#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 11:41:33 2022

@author: saumya
"""

import numpy as np
import math
from numpy.linalg import inv

# The following function forms the longitudinal state space with the following state vectors
# X = [V, alpha, q, theta]; U = [delta_e delta_i]' -- 9 variable output 
# alpha-q and V-gamma
def longitudinalss(dx):
    # Now to get the A matrix for the system
    # Since for longitudinal flight gamma + alpha = theta  -using theta is enough
    A = np.zeros((4,4))
    A[0,0] = dx[0,0];
    A[0,1] = dx[0,8];
    A[0,2] = dx[0,3];
    A[0,3] = dx[0,5];
    
    A[1,0] = dx[7,0]; #---> theta_derivative instead of gamma_derivative for now
    A[1,1] = dx[7,8];
    A[1,2] = dx[7,3];
    A[1,3] = dx[7,5];
    
    A[2,0] = dx[2,0];
    A[2,1] = dx[2,8];
    A[2,2] = dx[2,3];
    A[2,3] = dx[2,5];
    
    A[3,0] = dx[4,0]; 
    A[3,1] = dx[4,8];
    A[3,2] = dx[4,3];
    A[3,3] = dx[4,5];

# To assess effect of control variables obtain the B matrix but for now A is enough for stability matrix

    return A
# The following function forms the lateral state space with the following state vectors
# X = [beta, r, p, phi]; U = [delta_a, delta_r, delta_i]' -- 9 variable output 
# First mode is the pure roll mode
# Second mode is the side-slip oscillation related to beta and r
# Thirds mode is the spiral mode related to phi
def lateralss(dx):
    
    # Now to get the A matrix for the system
    A = np.zeros((4,4))
    A[0,0] = dx[1,1];
    A[0,1] = dx[1,6];
    A[0,2] = dx[1,4];
    A[0,3] = dx[1,7];
    
    A[1,0] = dx[5,1];
    A[1,1] = dx[5,6];
    A[1,2] = dx[5,4];
    A[1,3] = dx[5,7];
    
    A[2,0] = dx[3,1];
    A[2,1] = dx[3,6];
    A[2,2] = dx[3,4];
    A[2,3] = dx[3,7];

    A[3,0] = dx[6,1];
    A[3,1] = dx[6,6];
    A[3,2] = dx[6,4];
    A[3,3] = dx[6,7];

# To assess effect of control variables obtain the B matrix but for now A is enough for stability matrix

    return A