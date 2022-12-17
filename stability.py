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

def Jacobian(fix, x, CoefMatrix, atmo, g, PropWing):
    
    """The number of equations from equation.py are 10 but derivatives are obtained from the first eight
    equations. The number of columns is same as the number of fixed variables -(omega since only in last equation)
    plus the number of variables in x. This gives the total uncoupled Jacobian of the aircraft for analysis."""
    
    step_vector = 0.001; # Let
    
    length_x = len(x);
    dx = np.zeros(8,(length_x+3));
    tuple_1 = (CoefMatrix, atmo, g, PropWing)
    
    # Centered difference
    fix_high = fix;
    fix_high[0] = fix[0] + step_vector; # --V
    fix_low = fix;
    fix_low[0] = fix[0] - step_vector;
    dx[:,0] = (e.Constraints_DEP(x, fix_high, *tuple_1) - e.Constraints_DEP(x, fix_low, *tuple_1) )/(2*step_vector)
    
    fix_high = fix;
    fix_high[1] = fix[1] + step_vector; # --beta
    fix_low = fix;
    fix_low[1] = fix[1] - step_vector;
    dx[:,1] = (e.Constraints_DEP(x, fix_high, *tuple_1) - e.Constraints_DEP(x, fix_low, *tuple_1) )/(2*step_vector)
    
    fix_high = fix;
    fix_high[2] = fix[2] + step_vector; # --gamma
    fix_low = fix;
    fix_low[2] = fix[2] - step_vector;
    dx[:,2] = (e.Constraints_DEP(x, fix_high, *tuple_1) - e.Constraints_DEP(x, fix_low, *tuple_1) )/(2*step_vector)
    
    for i in range(length_x):
        
        x_high = x;
        x_high[i] = x[i] + step_vector; # --beta
        x_low = x;
        x_low[i] = fix[i] - step_vector;
        dx[:,i] = (e.Constraints_DEP(x_high, fix, *tuple_1) - e.Constraints_DEP(x_low, fix, *tuple_1) )/(2*step_vector)
        
        
    return dx

