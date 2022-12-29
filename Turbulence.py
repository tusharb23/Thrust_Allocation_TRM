#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 16:34:31 2022

@author: saumya
"""

import equation as e
import numpy as np
import AeroForcesDECOL

# Assuming Control Fixed --> i.e. B = 0
""" The velocity V used for the previous calculations is the axial Velocity of the aircraft in the aerodynamic
frame. This means that if the turbulence velocity is along the axis of the aircraft, it acts as a perturbation about the original Velocity.
Hence, now we need to obtain a relation between turbulence along the two perpendocular dircetions to the aircraft velocity to obtain the state space.
In that case, the equilibrium point about these two will be set at 0. So, from the previous stability analysis, we can establish, that the aircraft is stable for a perturbation about the Velocity. 
Now for vg and wg --> This is one approach, the other is to create trubulence components in three axis of the reference frame  of the aircraft and add it to the Velocity to see--> This will modify the angles as well."""

# Following the simplifications expressed in Cook, page 448 to 450
# ug--> deltaV = ug; i.e. V = V + ug
# vg--> deltabeta = vg
# wg--> deltaalpha = wg
# Similarly for pg, qg, and rg

def turbulence_ss(x, fix, CoefMatrix, atmo, g, PropWing):
    
    """The number of equations from equation.py are 10 but derivatives are obtained from the first eight
    equations. The number of columns is same as the number of fixed variables -(omega since only in last equation)
    plus the number of variables in x. This gives the total uncoupled Jacobian of the aircraft for analysis."""
    
    ug =0; vg=0; wg = 0;
    
    # Gusts centered about, ug =0, vg =0, wg =0---> If they are along the positive direction, their effect will be negative 
    step_vector = 0.001; # Let
    Velocity=(10,15,20,25,30,35)
    V = fix[0] # ---> To account for change in co-efficient matrix due to velocity change
    V1 = V - (ug + step_vector)
    V2 = V -(ug-step_vector)
    
    length_x = len(x);
    dx = np.zeros((8,3));
    Coef=AeroForcesDECOL.CoefInterpol(V, CoefMatrix, Velocity)
    tuple_1 = (CoefMatrix, atmo, g, PropWing)
    Coef=AeroForcesDECOL.CoefInterpol(V1, CoefMatrix, Velocity)
    tuple_V1 = (Coef, atmo, g, PropWing)
    Coef=AeroForcesDECOL.CoefInterpol(V2, CoefMatrix, Velocity)
    tuple_V2 = (Coef, atmo, g, PropWing)
    
    # Centered difference
    fix_high = np.zeros(np.size(fix));
    fix_high[0] = -(ug + step_vector); # --ug
    fix_low = np.zeros(np.size(fix));
    fix_low[0] = -(ug-step_vector);
    k = (e.Dynamics_DEP(x, fix+fix_high, *tuple_V1) - e.Dynamics_DEP(x, fix+fix_low, *tuple_V2) )/(2*step_vector)
    k = np.reshape(k,(10,1))
    dx[:,0] = k[:8,0]
    
    fix_high = np.zeros(np.size(fix));
    fix_high[1] = -(vg + step_vector)/V; # --ug
    fix_low = np.zeros(np.size(fix));
    fix_low[1] = -(vg-step_vector)/V;
    k = (e.Dynamics_DEP(x, fix+fix_high, *tuple_1) - e.Dynamics_DEP(x, fix+fix_low, *tuple_1) )/(2*step_vector)
    k = np.reshape(k,(10,1))
    dx[:,1] = k[:8,0]
    

    x_high = np.zeros(np.size(x));
    x_high[1] = -(wg + step_vector)/V; # --ug
    x_low = np.zeros(np.size(x));
    x_low[1] = -(wg-step_vector)/V;
    k = (e.Dynamics_DEP(x+x_high, fix, *tuple_1) - e.Dynamics_DEP(x+x_low, fix, *tuple_1) )/(2*step_vector)
    k = np.reshape(k,(10,1))
    dx[:,2] = k[:8,0]
    
    return dx
    
# V, beta, alpha, p, q, r, phi, theta
