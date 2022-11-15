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

def fitmodify(X, Velocities):
    # Here we modify the initial estimates for all the values after 20m/s as the optimization algorithm 
    X_new = np.zeros((17,11))
    X_new[:,0:5] = X[:,0:5]
    for i in range(17):
        f = np.polyfit(X[i,0:5],Velocities[0:5],1)
        for j in range(5,11):
            X_new[i,j]=np.polyval(f,(Velocities[j])) 
    return X_new


def Interpol(V, A1, A2, v1, v2):
    # Function to interpol any kind of variables in function of the velocity V
    # input : 
    # V : current velocity
    # A1, A2 : lower matrix, higher matrix
    # v1, v2 : velocities corresponding to matrices
    a = (A2-A1)/(v2-v1)
    x0 = A1 + a*(V-v1)
    return x0

def interpolateinitial(V, X, Velocities):
    
    # Modify x based on where the optimization algorithm fails
    # Here, assuming the algorithm fails after 20m/s
    #X = fitmodify(X, Velocities)
    #print(np.transpose(X))
    
    if V < Velocities[0]:
        return X[:,0]
    
    elif V > Velocities[-1]:
        return X[:,-1]
    else:
        exitcondition = 1
        length_v = len(Velocities)-1
        i = 0
        while exitcondition:
           
            if V == Velocities[i]:
                x0 = X[:,i]
                exitcondition = 0
            
            elif V > Velocities[i] and V < Velocities[i+1]:
                x0 = Interpol(V, X[:,i], X[:,i+1], Velocities[i], Velocities[i+1])
                exitcondition = 0  # exit
            
            else:
                i = i+1
                
            if i == length_v:  # security to exit the while
                print("Initial_value : Error in interpolating x0, returning 0")
                x0 = np.zeros((1,17))
                exitcondition = 0
    
    return x0
