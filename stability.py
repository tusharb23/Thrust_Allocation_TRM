#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 15:55:24 2022
<<<<<<< Updated upstream

@author: saumya
"""

"""=======
This module calls the equations.py and computes the total Jacobian of the aircraft.
From this jacobian, the eigenvalues are calculated and the stability is analysed for the entire aircraft.
Then inference is made regarding the dynamic stability of the system and different modes.
@author: saumya
"""

import equation as e
import numpy as np
import AeroForcesDECOL

def Jacobian(x, fix, CoefMatrix, atmo, g, PropWing):
    
    """The number of equations from equation.py are 10 but derivatives are obtained from the first eight
    equations. The number of columns is same as the number of fixed variables -(omega since only in last equation)
    plus the number of variables in x. This gives the total uncoupled Jacobian of the aircraft for analysis."""
    
    step_vector = 0.001; # Let
    Velocity=(10,15,20,25,30,35)
    V = fix[0] # ---> To account for change in co-efficient matrix due to velocity change
    V1 = V + 0.001
    V2 = V - 0.001
    
    length_x = len(x);
    dx = np.zeros((8,length_x+3));
    Coef=AeroForcesDECOL.CoefInterpol(V, CoefMatrix, Velocity)
    tuple_1 = (CoefMatrix, atmo, g, PropWing)
    Coef=AeroForcesDECOL.CoefInterpol(V1, CoefMatrix, Velocity)
    tuple_V1 = (Coef, atmo, g, PropWing)
    Coef=AeroForcesDECOL.CoefInterpol(V2, CoefMatrix, Velocity)
    tuple_V2 = (Coef, atmo, g, PropWing)
    
    # Centered difference
    fix_high = np.zeros(np.size(fix));
    fix_high[0] = step_vector; # --V
    fix_low = np.zeros(np.size(fix));
    fix_low[0] = - step_vector;
    k = (e.Dynamics_DEP(x, fix+fix_high, *tuple_V1) - e.Dynamics_DEP(x, fix+fix_low, *tuple_V2) )/(2*step_vector)
    k = np.reshape(k,(10,1))
    dx[:,0] = k[:8,0]
    
    
    fix_high = np.zeros(np.size(fix));
    fix_high[1] = step_vector; # --beta
    fix_low = np.zeros(np.size(fix));
    fix_low[1] = - step_vector;
    k = (e.Dynamics_DEP(x, fix+fix_high, *tuple_1) - e.Dynamics_DEP(x, fix+fix_low, *tuple_1) )/(2*step_vector)
    k = np.reshape(k,(10,1))
    dx[:,1] = k[:8,0]
    
    fix_high = np.zeros(np.size(fix));
    fix_high[2] = step_vector; # --gamma
    fix_low = np.zeros(np.size(fix));
    fix_low[2] = - step_vector;
    k = (e.Dynamics_DEP(x, fix+fix_high, *tuple_1) - e.Dynamics_DEP(x, fix+fix_low, *tuple_1) )/(2*step_vector)
    k = np.reshape(k,(10,1))
    dx[:,2] = k[:8,0]
    
    for i in range(length_x):
        
        x_high = np.zeros(np.size(x));
        x_high[i] = step_vector; 
        x_low = np.zeros(np.size(x));
        x_low[i] = - step_vector;
        k = (e.Dynamics_DEP(x+x_high, fix, *tuple_1) - e.Dynamics_DEP(x+x_low, fix, *tuple_1) )/(2*step_vector)
        k_trans = np.reshape(k,(10,1))
        dx[:,i+3] = k_trans[:8,0]
        
    return dx

# V, beta, gamma, alpha, p, q, r, phi, theta, delta_a, delta_e. delta_r, delta_x
# V, beta, alpha, p, q, r, phi, theta
# From the above dx get A and B and use it to obtain the step response of the system or something for a small step
# For control fixed response set B as zero


