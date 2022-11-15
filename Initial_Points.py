
""" This calls the opimization algorithm to obtain the curve for initial state value estimation """
""" - Modified Eric Nguyen Van and David Planas Andres 's work"""
""" For the pupose of Optimal Thrust Allocation """
""" - Modified by Saumya Sarawagi"""

# Import required inbuilt packages
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
import Optimization_Call as opt
import PattersonAugmented as PA


"""Create a function to compare when no DEP to define whether optimization required
    Analyze for different flight envelopes and graphs
    PAtterson Augmented
    Engines being Inoperative
    Time Analysis"""


# Atmospheric conditions for H = 0m
H = 0 # Altitude(m)
a = 340.3 # Velocity of Sound at 0m

# Used to obtain Aerodynamic Coefficients at different velocities 
Velocity=(10,15,20,25,30,35)
rho_vec=(1.225,1.225,1.225,1.225,1.225,1.225)
# Mach = [v/a for v in list(Velocities)]
Mach=np.ones((len(Velocity),1))*0.0001 # To have incompressible flow as the speds are considered very low

# Define DECOL Parameters
inop_eng = 0 # Number of engines that are inoperative
g = DECOLgeometry.data(inop_eng, r=0.113 / 2, rf=0.1865 / 2, zw=0.045)

# Constant Flight Parameters
V = 23.5  # Velocity (m/s)
M = V/a
beta = 0 / 180 * math.pi
gamma = 0 / 180 * math.pi  # math.atan(0/87.4)#/180*math.pi # 3% slope gradient # 6.88m/s vertical
R = 0  # in meters the turn radius
g.P_var = 8 * 14.4 * 4  # I*V*N_eng/2 
g.rho = 1.225 # Density of air at sea level as obtained from si2py.txt -> Eric Nguyen Van
g.atmo = [a, g.rho]
   

# Constraints 
ThrottleMax = 1  
ThrottleMin = 0.0001  
    # Control Surface Constraints 
deltaAMin = -30*math.pi/180 # Ailerons
deltaAMax = 30*math.pi/180
deltaEMin = -20*math.pi/180 # Elevators
deltaEMax = 20*math.pi/180
deltaRMin = -30*math.pi/180 # Rudder
deltaRMax = 30*math.pi/180
    # Angle Constraints
alphaMin = -5*math.pi/180
alphaMax = 25*math.pi/180
phiMin = -10*math.pi/180
phiMax = 10*math.pi/180
thetaMin = -30*math.pi/180
thetaMax = 30*math.pi/180
                                                                          
g.FlapDefl = 0  # Standard flap deflection (degrees) is 14 but for our purpose we consider 0 flap deflection
g.VelFlap = 12.5  # Maximum velocity at which flap are deployed (m/s)                                                      

# Used in the Patterson modulus for modelling stall in accordance with Antony Jameson's proposal
g.alpha_max = 10 / 180 * np.pi
g.alpha_max_fl = 10 / 180 * np.pi

# FLight measured Cd0:
g.CD0T = 0.0636         # Global one  extracted from flight not stab the file
                                                            
Velocities=(10,12,15,17,20,23,25,27,30,32,35)

X_initial = np.zeros((len(Velocities),17))
#compaare objective function value and constraints



# x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i] where delta i has 8 elelments
x0=np.array([5*math.pi/180, 0,0,0, 0.00, 0.0, 0.0, 0.0, 0.0]) # Assuming initial alpha to be 5degrees
bnds=( (alphaMin,alphaMax), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (phiMin,phiMax), (thetaMin,thetaMax), (deltaAMin,deltaAMax), (deltaEMin,deltaEMax), (deltaRMin, deltaRMax))
# Complete the vectors with engines:
eng_vec = np.array([0.4] * g.N_eng)
x0 = np.append(x0, eng_vec)
## General formulation:
bnds_eng = ((ThrottleMin, ThrottleMax), (ThrottleMin, ThrottleMax))
for i in range(int(g.N_eng / 2)):
    bnds = bnds + bnds_eng


for i in range(len(Velocities)):
    V = Velocities[i]
    temp = opt.opt_call(g, V, beta, gamma, R, Mach, Velocity, x0, bnds)
    X_initial[i,:] = temp
    
np.savetxt("Initial.txt",X_initial)
    


