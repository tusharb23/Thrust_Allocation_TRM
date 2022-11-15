""" Optimization Algorithm Called as a function.
Created by Saumya Sarawagi based on work done by Eric Nguyen Van and David Planas Andres"""


import numpy as np
import math
import scipy.linalg
import scipy.io #input/output with matlab
from scipy.optimize import  minimize
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import random
# Import DECOL packages
import AeroForcesDECOL
import DECOLgeometry
import ReadFileUtils
import equation as e
import PattersonAugmented as PA

def opt_call(g, V, beta, gamma, R, Mach, Velocity, x0, bnds):
    # Extracting the aerodynamic coefficients
    path = 'DECOL_STAB/'  
    filenameNoFin = [path + '_FinLess_Vinf10000.stab',
                 path + '_FinLess_Vinf15000.stab',
                 path + '_FinLess_Vinf20000.stab',
                 path + '_FinLess_Vinf25000.stab',
                 path + '_FinLess_Vinf30000.stab',
                 path + '_FinLess_Vinf35000.stab']
    MatrixNoFin = ReadFileUtils.ReadStabCoef(filenameNoFin)
    CoefMatrix=g.NicolosiCoef(MatrixNoFin[:,1:], Mach)
    Coef=AeroForcesDECOL.CoefInterpol(V, CoefMatrix, Velocity)

# Defining the Propulsion & Wing Interaction
    g.PolarFlDeflDeg = 5
    g.PolarAilDeflDeg = 5
    PropPath = "DECOL_FEM/"
    PropFilenames = {'fem':[PropPath+"_FinLess_Vinf10000.0"],
                 'AirfoilPolar':PropPath+"S3010_XTr10_Re350.txt",
                 'FlapPolar':PropPath+"S3010_XTr10_Re350_fl5.txt",
                 'AileronPolar':PropPath+"S3010_XTr10_Re350_fl5.txt"} # format for prop file : [[Cldist=f(M)],polar clean airfoil, polar flap, polar aile]
    PW = PA.PropWing(g,PropFilenames)
    #PW.AoAZero[:,-1] = PW.AoAZero[:,-1] + 3.2/180*np.pi # Correction for angle of incidence of wing
    PW.AoAZero[:,0] = PW.AoAZero[:,0]*10**(-3)
    PW.CLslope[:,0] = PW.CLslope[:,0]*10**(-3)
    PW.AoAZero[:,1] = PW.AoAZero[:,1]*10**(-6)
    PW.CLslope[:,1] = PW.CLslope[:,1]*10**(-6)         # In order to express in S.I. units as DECOL.vsp3 yields in mm
    PW.AoAZero[:,2] = PW.AoAZero[:,2]*10**(-3)
    PW.CLslope[:,2] = PW.CLslope[:,2]*10**(-3)
    PW.DeltaCL_a_0 = 1 #CL_alpha correction factor"""

    g.nofin = False
    g.DisplayPatterInfo = False
    
# --- imposed conditions ---
# fix = [V, beta, gamma, omega, H]
    if R == 0:
        omega = 0
    else:
        omega = V * math.cos(gamma) / R

    fixtest = np.array([V, beta, gamma, omega])
# put everything in tuples for passing to functions
    diccons = (np.copy(fixtest), np.copy(Coef), g.atmo, g, PW)  # fix, CoefMatrix,Velocities, rho, g
    dicfobj = (np.copy(fixtest), g.rho, g)
    maxit = 100
    tolerance = 1e-5
# TRM is working
    t0 = datetime.now()
    k = minimize(e.fobjectivePower, np.copy(x0), args=dicfobj, method = 'trust-constr', bounds=bnds,
                         constraints={'type': 'eq', 'fun': e.Constraints_DEP, 'args': diccons},
                         options={'maxiter': maxit, 
                                  'disp': True}, tol=tolerance)
    t1 = datetime.now()
    print("Evaluation_time :" , t1-t0)
    return k.x
    


