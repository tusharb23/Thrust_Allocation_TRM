""" The main_DEP.py calls the necessary functions from other packages and executes the following functions:
    1) Get all the aerodynamic coefficients for the aircraft at the required velocity and flight parameters
    2) Include aero-propulsive interactions
    3) Calculate the thrust and moment needed from the flight equations
    4) Optimally distribute thrust among the engines and control surfaces """
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
import PattersonAugmented as PA


"""Create a function to compare when no DEP to define whether optimization required
    Analyze for different flight envelopes and graphs
    PAtterson Augmented
    Engines being Inoperative
    Time Analysis"""


# Atmospheric conditions for H = 0m
H = 0 # Altitude(m)
a = 340.3 # Velocity of Sound at 0m
rho = 1.225 # Density of air at sea level as obtained from si2py.txt -> Eric Nguyen Van
atmo = [a, rho]

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
beta = 10 / 180 * math.pi
gamma = 0 / 180 * math.pi  # math.atan(0/87.4)#/180*math.pi # 3% slope gradient # 6.88m/s vertical
R = 0  # in meters the turn radius
g.P_var = 8 * 14.4 * 4  # I*V*N_eng/2    

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
                                                            
# Extracting the aerodynamic coefficients
#
#Velocities=(10,12,15,17,20,23,25,27,30,32,35)
#
#P = np.zeros(np.size(Velocities))
#color = ['cyan', 'lightblue', 'deepskyblue', 'steelblue', 'dodgerblue', 'red', 'royalblue', 'navy', 'blue', 'slateblue', 'indigo']
#plt.Figure()
#for count in range(np.size(Velocities)):
#compaare objective function value and constraints

Power_minimum = 148*8;
x_min = np.zeros(17)
con = np.ones(10)*10e-1


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

# x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i] where delta i has 8 elelments
    #x0=np.array([5*math.pi/180, 0,0,0, 0.00, 0.0, 0.0, 0.0, 0.0]) # Assuming initial alpha to be 5degrees
   
"""x0 = np.array([random.uniform(alphaMin,alphaMax),random.uniform(-0.2,0.2), random.uniform(-0.2,0.2), random.uniform(-0.2,0.2), random.uniform(phiMin,phiMax), random.uniform(thetaMin,thetaMax), random.uniform(deltaAMin,deltaAMax), random.uniform(deltaEMin,deltaEMax), random.uniform(deltaRMin, deltaRMax)])
eng_vec = np.array([random.uniform(0, 1)] * g.N_eng)
x0 = np.append(x0, eng_vec)""" # --- Define x0 using the polynomial fit instead for a cruise scenario instead of randomly generated initial points  
    
bnds=( (alphaMin,alphaMax), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (phiMin,phiMax), (thetaMin,thetaMax), (deltaAMin,deltaAMax), (deltaEMin,deltaEMax), (deltaRMin, deltaRMax))
#Constrains Used by Eric & David
#phimax = 10  # in degree the max bank angle authorized
#alphamax = 25  # in degree to adjust if it is below 71m/s
#deltaRmax = 30  # in degree
#bnds=( (-5*math.pi/180,alphamax*math.pi/180), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (-phimax/180*math.pi,phimax/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-20/180*math.pi,20/180*math.pi), (-deltaRmax/180*math.pi,deltaRmax/180*math.pi))

# Complete the vectors with engines:
    #eng_vec = np.array([0.4] * g.N_eng)
    #x0 = np.append(x0, eng_vec)
## General formulation:
bnds_eng = ((ThrottleMin, ThrottleMax), (ThrottleMin, ThrottleMax))
for i in range(int(g.N_eng / 2)):
        bnds = bnds + bnds_eng
    
# --- imposed conditions ---
# fix = [V, beta, gamma, omega, H]
if R == 0:
        omega = 0
else:
        omega = V * math.cos(gamma) / R

fixtest = np.array([V, beta, gamma, omega])
# put everything in tuples for passing to functions
diccons = (np.copy(fixtest), np.copy(Coef), atmo, g, PW)  # fix, CoefMatrix,Velocities, rho, g
dicfobj = (np.copy(fixtest), rho, g)
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
print(k)
    
if np.sum(np.absolute((e.Constraints_DEP(k.x,*diccons))))<np.sum(np.absolute(con)) and  Power_minimum>k.fun:
            Power_minimum =k.fun
            x_min = k.x
            con = e.Constraints_DEP(k.x,*diccons)

def printx(x, fix, atmo, g, PW):
        V = fix[0]
        alpha = x[0]/math.pi*180
        beta = fix[1]/math.pi*180
        pqr = x[1:4]/math.pi*180
        phi = x[4]/math.pi*180
        theta = x[5]/math.pi*180
        da = x[6]/math.pi*180
        de = x[7]/math.pi*180
    
        print("\nState vector value:")
        print("V= {0:0.2f}m/s, alpha = {1:0.2f}\xb0, beta={2:0.2f}\xb0, phi={3:0.2f}\xb0, theta={4:0.2f}\xb0".format(V, alpha, beta, phi, theta))
        print("p={0:0.4f}\xb0/s q={1:0.4f}\xb0/s r={2:0.4f}\xb0/s".format(*pqr))
        print("da={0:0.2f}\xb0, de= {1:0.2f}\xb0".format(da,de))

        V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * fix[1] + g.wingsweep) - x[3] * g.PosiEng

        if g.IsPropWing:
            if V <= g.VelFlap:
                PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], g.FlapDefl, g, False, fix[1], x[1], V, x[3])
            else:
                PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], 0, g, False, fix[1], x[1], V, x[3])

        if g.nofin==False:
            print("dr = {0:0.2f}\xb0".format(x[8]/math.pi*180))

printx(k.x, fixtest, atmo,g,PW)
print(k.fun)
    
# check if constraints are validated
constraints_calc=e.Constraints_DEP(k.x,*diccons)
print("\nConstraints")
print(constraints_calc)



"""Engine = [1,2,3,4,5,6,7,8]
del_x = k.x[9:]
plt.plot(Engine,del_x, color[count],label =str(Velocities[count])+"m/s")
plt.title("Fraction of maximum Power at each engine")
plt.xlabel("Engine")
plt.ylabel("Delta_x")    
plt.legend(loc = 'upper right', fontsize = "xx-small")
plt.grid()
plt.show()
Z = [397.53995125454855, 438.57296670349535, 419.4142735202689, 440.18643949094616, 313.6531953322101, 410.14216342352347, 482.0370351986316, 600.5445635914646, 713.4563025540505, 859.1137508254681, 1121.2532233356421]
plt.scatter(Velocities, P, c='g', label = "Updated model")
plt.plot(Velocities, P, c='b')
plt.scatter(Velocities, Z, c='r', label="Old model")
plt.title("Power(W) variation with Velocity(m/s)")
plt.xlabel("Velocity m/s")
plt.ylabel("Power W")
plt.legend(loc = 'upper right', fontsize = "xx-small")
plt.grid()
plt.show()
print(Power_minimum)
print(x_min)
print(con)
Engine = [1,2,3,4,5,6,7,8]
del_x = x_min[9:]
plt.plot(Engine,del_x)
plt.title("Fraction of maximum Power at each engine")
plt.xlabel("Engine")
plt.ylabel("Delta_x")    
plt.legend(loc = 'upper right', fontsize = "xx-small")
plt.grid()
plt.show()"""
