""" The AeroForcesDECOL.py calculates the aerodynamic forces and moments on the DECOL aircraft"""
""" - Modified Eric Nguyen Van and David Planas Andres 's work"""
""" For the pupose of Optimal Thrust Allocation """
""" - Modified by Saumya Sarawagi"""

import numpy as np


def Interpol(V, A1, A2, v1, v2):                                                                                        
    # Function to interpol any kind of variables in function of the velocity V
    # Input : 
    # V : current velocity
    # A1, A2 : lower matrix, higher matrix
    # v1, v2 : velocities corresponding to matrices
    a=(A2-A1)/(v2-v1)
    b=A1-a*v1
    Areturn=a*V+b
    return Areturn


def InterpolRho(V, rho, v):                                                                                            
    # Function to interpol Rho, since altitude is function of Velocity
    if V<v[0] :
        return rho[0]
    
    elif V>v[-1]:
        return rho[-1]
    else:
        exitcondition=1
        length_v=len(v)-1
        i=0
        while exitcondition :
           
            if V==v[i]:
                rhoreturn=rho[i]
                exitcondition=0
            
            elif V>v[i] and V<v[i+1]:
                rhoreturn=Interpol(V, rho[i], rho[i+1], v[i], v[i+1])
                exitcondition=0 #exit
            
            else:
                i=i+1
                
            if i==length_v: #security to exit the while
                print("AeroForces : Error in interpolating rho, returning 0")
                rhoreturn=0
                exitcondition=0
    
    return rhoreturn


def CoefInterpol( V, A, v):                                                                                             
    # A, function takes a numpy array composed of all matrices A [A1; A2; ...], all array types!!
    # v, an array of corresponding velocity
    # V, the velocity at which to compute the coef
    # Size of each matrix 6 row, m column
    row=6
    nmatrix=len(A[:,1])/6
    if nmatrix == float:
        # Error
        print('Aero forces: general coef matrix not a multiple of 6')
        return 0
    
    if V<v[0] :
        # Ill posed problem, the velocity is below smallest velocity for coef
        # Use first matrix but print warning
        print("WARNING : velocity, V = {0:0.2f} is below first ref velocity for coef v = {1:0.2f}".format(V,v[0]))
        print("Continue with first matrix")
        return A[0:row,:]
    
    elif V>v[-1]:
        # Same idea, the velocity is greater than max used to determine coef. Results in flight faster than cruise
        # Use last matrix but print warning
        print("WARNING : velocity, V = {0:0.2f} is higher than last ref velocity for coef v = {1:0.2f}".format(V,v[-1]))
        print("Continue with last matrix")
        return A[-6:]
    
    elif V<0.2:                                                                                                         
        # Hard coded use the first matrix with flaps
        return A[0:row,:]
        
    else : # Otherwise interpolate
        exitcondition=1
        length_v=len(v)-1
        i=1
        while exitcondition :
           
            if V==v[i-1]:
                Areturn=A[(i-1)*row:(i-1)*row+row]
                exitcondition=0
                
            if V==v[i]:
                Areturn=A[i*row:i*row+row]
                exitcondition=0
            
            elif V>v[i] and V<v[i+1]:
                Areturn=Interpol(V, A[i*row:i*row+row,:], A[(i+1)*row:(i+1)*row+row,:], v[i], v[i+1])
                exitcondition=0 #exit
                
            elif V<v[i] and V>v[i-1]:
                Areturn=Interpol(V, A[(i-1)*row:(i-1)*row+row,:], A[(i)*row:(i)*row+row,:], v[i-1], v[i])
                exitcondition=0 #exit
            
            else:
                i=i+1
                
            if i==length_v+1: #security to exit the while
                print("!!! FAILURE !!! AeroForces : Error in interpolation, returning 0")
                Areturn=0
                exitcondition=0

    return Areturn


def CalcForce_aeroframe_DEP(V, CoefMatrix, x, dx, atmo, g, PropWing):
    """ Function to compute aerodynamic forces in the velocity frame (aero frame)                                       
    for the DEP configuration. The propulsion force and moments are not computed here
    Since V is fixed, Coef Matrix must be calculated before
    Can handle DEP, with and without rudder, 2 or more engines
    """
    rho = atmo[1]
    a_sound = atmo[0]
    beta=x[1]                                                                                                          
    p=x[2]                                                                                                              
    r=x[4]

    #Compute aero forces
    # here x must be of the form (alpha, beta, p, q, r, da, de, dr)
    # set non dim for p,q,r
    nonDim=np.ones(len(x))
    nonDim[2]=g.b/(2*V)
    nonDim[3]=g.c/(2*V)
    nonDim[4]=g.b/(2*V)
    x=x*nonDim


    # If prop wing, aileron contributions are included in patterson
    if g.IsPropWing:
        CoefMatrix[3:6,5] = np.zeros(3)
    #    F=np.dot(CoefMatrix,x[0:7]) # commented form, modification to account for symmetric drag increase of side slip
    F=np.zeros((3))
    M=np.zeros((3))
    xsym=np.copy(x)
    xsym[1]=abs(xsym[1]) # make beta always positive since derivatives have already correct sign for drag and lift only
    xsym[5]=abs(xsym[5]) # make ailerons deflection always positive for drag increase and lift decrease
    xsym[-1]=abs(xsym[-1]) # make rudder deflection always positive for drag increase and lift decrease


    F[0]=np.dot(CoefMatrix[0,1:],xsym[1:])
    F[1]=np.dot(CoefMatrix[1],x)
    F[2]=np.dot(CoefMatrix[2],xsym)
    M=np.dot(CoefMatrix[3:6,:],x)
    
    # Drag force quadratic in alpha
    # DragQuad = F[0]+0.8029*x[0]**2 + x[0] * 0.12537 + 0.006434
    # Drag force quadratic in CL
    DragQuad = F[0] + 0.99747*x[0]**2 + 0.18081*x[0] + 0.02977


    #Additional Parameters.
    g.Hor_tail_coef_vol = (g.Sh*g.lv) / (g.S*g.c)
    g.eta_downwash = 0.8
    # g.X_CA_wb = g.x_cg/g.c - (   (CoefMatrix[4,0] + g.aht*g.Hor_tail_coef_vol *0.8)/ (CoefMatrix[2,0]-g.aht)      )    #--> Calculated so that downwash *ratio of dynamic pressures is 0.8
    g.SM_CA_wingfuselage =  (CoefMatrix[4,0] + g.aht*g.Hor_tail_coef_vol *g.eta_downwash)/ (CoefMatrix[2,0]-g.aht)      # It bothers me that if you calculate the CA_wingbody does not really match with geometry,
                                                                                                                        # maybe because center of gravity is not well estimated? IN DECOL for having a good value for
                                                                                                                        # CA_wingbody (25% of chord should give x_ca_wb=0.6725),
                                                                                                                        # Center of gravity should be around  0.621 m, is 0.713 now. Anyway calculation of
                                                                                                                        # Cm alpha only cares about their difference, or static margin.
    # Cl_alpha with interaction calculus
    alpha_1 = 0 * np.pi/180
    alpha_2 = 2 * np.pi/180
    CL_alpha_interaction = (   (PropWing.CalcCoef(dx,V/a_sound, atmo, alpha_2,x[5],g.FlapDefl,g,beta,p,V,r)[0] + g.aht*alpha_2 ) - (PropWing.CalcCoef(dx,V/a_sound, atmo, alpha_1,x[5],g.FlapDefl,g,beta,p,V,r)[0] + g.aht*alpha_1 ) )/ (alpha_2-alpha_1)


    # Calculates the new Cm_alpha_aero_interaction
    g.Cm_alpha_aero_interaction =( 1 + ((CL_alpha_interaction - CoefMatrix[2,0]) * g.SM_CA_wingfuselage)/CoefMatrix[4,0] ) * CoefMatrix[4,0]   #this formula is valis supossing that downwash, c.gravity and tail-wing pressure ratio
    # Does not change when implementing DEP


# Function computes the Aerodynamic Force and Moment for the DEP aircraft for both cases : with and without propuslision wing interaction
    if g.IsPropWing:
        
        CoefMatrix[4,0]=g.Cm_alpha_aero_interaction
        M=np.dot(CoefMatrix[3:6,:],x)

        if V<=g.VelFlap :
            # CLCl = PropWing.CalcCoef(dx,V, atmo, x[0],x[5],g.FlapDefl,g)
            CLCl = PropWing.CalcCoef(dx, V/a_sound, atmo, x[0], x[5], g.FlapDefl, g,beta,p,V,r)
            # Compute effects of other variables
            F[2]=np.dot(CoefMatrix[2][1:],xsym[1:])
            
            if len(CLCl)>2 and g.IsPropWingDrag:
                # Drag is computed by patterson, add contribution of other variables (than alpha and dx)
                Fbody=np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T,F[1],-F[2]-CLCl[0]]) # add alpha=0 coefficients
                Moment=M+np.array([CLCl[1],g.Cm0_fl,CLCl[4]])                                                           

            else:
                # Add CL only, due to prop blowing
                Fbody=np.array([-DragQuad-g.Cd0_fl-g.CD0T,F[1],-F[2]-CLCl[0]]) # add alpha=0 coefficients
                # Add roll effect
                Moment=M+np.array([CLCl[1],g.Cm0_fl,CLCl[4]])                                                           

        else :
            CLCl = PropWing.CalcCoef(dx,V/a_sound, atmo, x[0],x[5],0,g,beta,p,V,r)
            # Compute effects of other variables
            F[0]=np.dot(CoefMatrix[0][1:],xsym[1:])
            F[2]=np.dot(CoefMatrix[2][1:],xsym[1:])
            # By default lift and drag are computed here
            Fbody=np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T,F[1],-F[2]-CLCl[0]]) # add alpha=0 coefficients
            # Add roll effect
            Moment=M+np.array([CLCl[1],g.Cm0,CLCl[4]])  
                                                                
    else:
        if V<=g.VelFlap :
            Fbody=np.array([-DragQuad-g.Cd0_fl,F[1],-F[2]-g.CL0_fl-g.CL0]) # add alpha=0 coefficients                   
            Fbody[0]=Fbody[0]-g.CD0T                                                                                    
            Moment=M+np.array([0,g.Cm0_fl,0]) 
                                                                          
        else:
            Fbody=np.array([-DragQuad-g.CD0T,F[1],-F[2]-g.CL0]) # add alpha=0 coefficients                              
            Moment=M+np.array([0,g.Cm0,0])                                                                              
 
    # Add contribution of horizontal tail
    if g.IsPropWing:                        
        Fbody[2]=Fbody[2]-g.aht*x[0]

    g.TotalDrag = abs(Fbody[0])
    g.Lift = abs(Fbody[2])
    g.lift2drag = abs(Fbody[2]/Fbody[0])
    Fbody=0.5*V**2.0*rho*g.S*Fbody
    Moment=0.5*V**2.0*rho*g.S*Moment*np.array([g.b,g.c,g.b])
    
    return np.append(Fbody, Moment)
