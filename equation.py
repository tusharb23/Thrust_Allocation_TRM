""" The purpose of equation.py is to implement the control laws to obtain the thrust and moment components required."""
""" - Modified Eric Nguyen Van and David Planas Andres 's work"""
""" For the pupose of Optimal Thrust Allocation """
""" - Modified by Saumya Sarawagi"""


import numpy as np
import math
from numpy.linalg import inv
import AeroForcesDECOL as AeroForces
from Propeller_Data import Prop

def Constraints_DEP(x, fix, CoefMatrix, atmo, g, PropWing):
    """function defining constraints for power minimization
    inputs:
        -x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
        x is the state to determine
        length of x except the propulsion levels is 8
        -fix = [V, beta, gamma, omega]
        fix is the vector of parameters whom are fixed by the user
    """

    rho = atmo[1]

    # --- Now prepare variables for equations ---
    V=fix[0]
    alpha=x[0]
    beta=fix[1]
    gamma=fix[2]
    omega=fix[-1]
    p=x[1]
    q=x[2]
    r=x[3]
    phi=x[4]
    theta=x[5]
    I=np.array([ [g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz] ])
    
    # --- Compute aerodynamic forces ---
    # Here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
    sub_vect=np.array([alpha,beta,p,q,r])
    sub_vect=np.append(sub_vect,[x[6],x[7],x[8]]) # Rudder is allowed
    # The velocity seen by all engines is not same due to sweep and roll, as a result, the velocity vector acting on these are different
    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * beta + g.wingsweep) - r * g.PosiEng

    Fx_body=g.Thrust(x[-g.N_eng:],V_vect)
    Mx = g.Torque(x[-g.N_eng:],V_vect)
    # Convert Thrust in body frame to thrust in aerodynamic frame
    Tab = [[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)],
           [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(alpha)*np.sin(beta)],
           [-np.sin(alpha), 0, np.cos(alpha)]]
    Fx_aero = np.matmul(Tab,Fx_body)
    Mx_aero = Tab @ Mx
    
    # Convert thrust in Tc for patterson
    Tc = g.Thr
    F=AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    
    # Now sum up the constraints:
    sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank=np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi) 

    A=np.zeros(10+g.inop)
    """
    A0 = x
    A1 = y
    A2 = z
    A3 = l
    A4 = m
    A5 = n
    A6 = phi
    A7 = theta
    A8 = gamma
    A9 = Omega
    """
    A[0]=-9.81*np.sin(gamma)+F[0]/g.m+Fx_aero[0]/g.m # Velocity_derivative
    A[1]=(p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F[1]/(g.m*V) + Fx_aero[1]/(g.m*V) # Beta_derivative
    A[2]=-(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta)+ 9.81*cosbank/(V*np.cos(beta)) + F[2]/(g.m*V*np.cos(beta))+Fx_aero[2]/(g.m*V*np.cos(beta)) # Alpha_derivative
    A[3:6] = np.dot(inv(I), np.array([Mx_aero[0], Mx_aero[1], Mx_aero[2]])+F[3:6]-np.cross(np.array([p, q, r]), np.dot(I, np.array([p, q, r])))) #p,q,r_derivative
    A[6]=p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta) # phi_derivtive
    A[7]=q*math.cos(phi) -r*math.sin(phi) # theta_derivative
    A[8]=-np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta) # No derivative
    A[9]=-omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta) # No derivative

    #if beta  0, make all the thrust equal
    
    for i in range(g.inop):
        A[-1-i]=x[-1-i]
    """if beta == 0:                                                                                 #For obligating all the engines to have the same thrust
        #no DEP with original twin or N engines; all engines have the same thrust
        D = np.copy(A)
        for i in range(g.N_eng-g.inop-1):
            AAd = x[-g.N_eng]-x[-g.N_eng+i+1]
            D = np.append(D, [AAd])
        return D
    else:
        return A"""
    return A


def Dynamics_DEP(x, fix, CoefMatrix, atmo, g, PropWing):
    """function defining constraints for power minimization
    inputs:
        -x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
        x is the state to determine
        length of x except the propulsion levels is 8
        -fix = [V, beta, gamma, omega]
        fix is the vector of parameters whom are fixed by the user
    """
    # Since during perturbation of velocity, we need to recalculate the co-efficient matrix
    # This is accounted where the Jacobian is formed in the section for V

    rho = atmo[1]

    # --- Now prepare variables for equations ---
    V=fix[0]
    alpha=x[0]
    beta=fix[1]
    gamma=fix[2]
    omega=fix[-1]
    p=x[1]
    q=x[2]
    r=x[3]
    phi=x[4]
    theta=x[5]
    I=np.array([ [g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz] ])
    
    # --- Compute aerodynamic forces ---
    # Here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
    sub_vect=np.array([alpha,beta,p,q,r])
    sub_vect=np.append(sub_vect,[x[6],x[7],x[8]]) # Rudder is allowed
    # The velocity seen by all engines is not same due to sweep and roll, as a result, the velocity vector acting on these are different
    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * beta + g.wingsweep) - r * g.PosiEng

    Fx_body=g.Thrust(x[-g.N_eng:],V_vect)
    Mx = g.Torque(x[-g.N_eng:],V_vect)
    # Convert Thrust in body frame to thrust in aerodynamic frame
    Tab = [[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)],
           [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(alpha)*np.sin(beta)],
           [-np.sin(alpha), 0, np.cos(alpha)]]
    Fx_aero = np.matmul(Tab,Fx_body)
    Mx_aero = Tab @ Mx
    
    # Convert thrust in Tc for patterson
    Tc = g.Thr
    F=AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    
    # Now sum up the constraints:
    sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank=np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi) 
    A=np.zeros(10+g.inop)
    
    # Since perturbation about trim point - since perturbation about theta -- should not transfer to gamma
    #sin_gamma = np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta) # No derivative
    #omega = (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta) # No derivative

    A[0]=-9.81*np.sin(gamma)+F[0]/g.m+Fx_aero[0]/g.m # Velocity_derivative
    A[1]=(p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F[1]/(g.m*V) + Fx_aero[1]/(g.m*V) # Beta_derivative
    A[2]=-(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta)+ 9.81*cosbank/(V*np.cos(beta)) + F[2]/(g.m*V*np.cos(beta))+Fx_aero[2]/(g.m*V*np.cos(beta)) # Alpha_derivative
    A[3:6] = np.dot(inv(I), np.array([Mx_aero[0], Mx_aero[1], Mx_aero[2]])+F[3:6]-np.cross(np.array([p, q, r]), np.dot(I, np.array([p, q, r])))) #p,q,r_derivative
    A[6]=p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta) # phi_derivtive
    A[7]=q*math.cos(phi) -r*math.sin(phi) # theta_derivative 
   
    return A


def fobjectivePower(x, fix, rho, g):
    
    V_vect = np.ones(g.N_eng) * fix[0] * np.cos((-np.sign(g.PosiEng)) * fix[1] + g.wingsweep) - x[3] * g.PosiEng
    Power = np.sum(g.Propeller.PropPower(x[-g.N_eng:], V_vect))
    return Power


    
