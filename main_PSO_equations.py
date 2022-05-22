#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 22 10:15:01 2022

@author: saumya
"""

# Import required inbuilt packages
import numpy as np
import math
import scipy.linalg
import scipy.io #input/output with matlab
from datetime import datetime
from scipy.integrate import odeint
import random as rnd
from random import uniform

# Import DECOL packages
import AeroForcesDECOL
import DECOLgeometry
import ReadFileUtils
import PattersonAugmented as PA
import PSO_implemnted as pso

def main():
    
    # Atmospheric conditions for H = 0m
    H = 0 # Altitude(m)
    a = 340.3 # Velocity of Sound at 0m
    rho = 1.225 # Density of air at sea level as obtained from si2py.txt -> Eric Nguyen Van
    atmo = [a, rho]

    # Used to obtain Aerodynamic Coefficients at different velocities 
    Velocities=(10,15,20,25,30,35)
    rho_vec=(1.225,1.225,1.225,1.225,1.225,1.225)
    # Mach = [v/a for v in list(Velocities)]
    Mach=np.ones((len(Velocities),1))*0.0001 # To have incompressible flow as the speds are considered very low
    
    # Define DECOL Parameters
    inop_eng = 0 # Number of engines that are inoperative
    g = DECOLgeometry.data(inop_eng, r=0.113 / 2, rf=0.1865 / 2, zw=0.045)
    
    # Constant Flight Parameters
    V = 23.5   # Velocity (m/s)
    M = V/a
    beta = 0 / 180 * math.pi
    gamma = 0 / 180 * np.pi  # math.atan(0/87.4)#/180*math.pi # 3% slope gradient # 6.88m/s vertical
    R = 0  # in meters the turn radius
    g.P_var = 8 * 14.4 * 4  # I*V*N_eng/2    

    # Constraints 
    ThrottleMax = 1  
    ThrottleMin = 0.0001  
    phiMin = -30*math.pi/180
    phiMax = 30*math.pi/180
    thetaMin = -30*math.pi/180
    thetaMax = 30*math.pi/180
                                                                          
    g.FlapDefl = 0  # Standard flap deflection (degrees) is 14 but for our purpose we consider 0 flap deflection
    g.VelFlap = 12.5  # Maximum velocity at which flap are deployed (m/s)                                                      

    # Used in the Patterson modulus for modelling stall in accordance with Antony Jameson's proposal
    g.alpha_max = 10 / 180 * np.pi
    g.alpha_max_fl = 10 / 180 * np.pi

    # FLight measured Cd0:
    g.CD0T = 0.0636         # Global one  extracted from flight not stab the file
    
    # --- List all .stab file from vsp aero and read the coeff ----
    path = 'DECOL_STAB/'  # 'home/e.nguyen-van/Documents/DECOLStability&Analysis/DECOLDATA/DECOLGeom_DegenGeom_6_3_18h14
    filenameNoFin = [path + '_FinLess_Vinf10000.stab',
                     path + '_FinLess_Vinf15000.stab',
                     path + '_FinLess_Vinf20000.stab',
                     path + '_FinLess_Vinf25000.stab',
                     path + '_FinLess_Vinf30000.stab',
                     path + '_FinLess_Vinf35000.stab']
    MatrixNoFin = ReadFileUtils.ReadStabCoef(filenameNoFin)
    # copy the matrix to avoid error and keep track
    Matrix = np.copy(MatrixNoFin)
    CoefMatrix=g.NicolosiCoef(Matrix[:, 1:], Mach)
    Coef_base = AeroForcesDECOL.CoefInterpol(V, CoefMatrix, Velocities)
    g.Matrix_no_tail_terms = AeroForcesDECOL.CoefInterpol(V, Matrix[:, 1:], Velocities)
    g.PolarFlDeflDeg = 5
    g.PolarAilDeflDeg = 5
    PropPath = "DECOL_FEM/"
    PropFilenames = {'fem':[PropPath+"_FinLess_Vinf10000.0"],
                     'AirfoilPolar':PropPath+"S3010_XTr10_Re350.txt",
                     'FlapPolar':PropPath+"S3010_XTr10_Re350_fl5.txt",
                     'AileronPolar':PropPath+"S3010_XTr10_Re350_fl5.txt"} # format for prop file : [[Cldist=f(M)],polar clean airfoil, polar flap, polar aile]
    PW = PA.PropWing(g, PropFilenames)
    #PW.AoAZero[:,-1] = PW.AoAZero[:,-1] + 3.2/180*np.pi #correction for angle of incidence of wing
    PW.AoAZero[:, 0] = PW.AoAZero[:, 0]*10**(-3)
    PW.CLslope[:, 0] = PW.CLslope[:, 0]*10**(-3)
    PW.AoAZero[:, 1] = PW.AoAZero[:, 1]*10**(-6)
    PW.CLslope[:, 1] = PW.CLslope[:, 1]*10**(-6)         #TO CORRECT UNITS AS DECOL.vsp3 is in mm
    PW.AoAZero[:, 2] = PW.AoAZero[:, 2]*10**(-3)
    PW.CLslope[:, 2] = PW.CLslope[:, 2]*10**(-3)
    PW.DeltaCL_a_0 = 1  # CL_alpha correction factor
    
    g.nofin = False
    g.DisplayPatterInfo = False
    
    # --- imposed conditions ---
    # fix = [V, beta, gamma, omega, H]
    if R == 0:
        omega = 0
    else:
        omega = V * math.cos(gamma) / R
    # Implement the PSO Algorithm
    fix = np.array([V, beta, gamma, omega])
    k = pso.minimize(fix, CoefMatrix, atmo, g, PW)
    print(k)
      
    
if __name__=='__main__':
    main()
    
        
    
    
        

