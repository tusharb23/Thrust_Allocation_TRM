""" The AeroForcesDECOL.py calculates the aerodynamic forces and moments on the DECOL aircraft"""
""" - Modified Eric Nguyen Van and David Planas Andres 's work"""
""" For the pupose of Optimal Thrust Allocation """
""" - Modified by Saumya Sarawagi"""

import numpy as np


def Interpol(V, A1, A2, v1, v2):
    # Function to interpol any kind of variables in function of the velocity V
    # input : 
    # V : current velocity
    # A1, A2 : lower matrix, higher matrix
    # v1, v2 : velocities corresponding to matrices
    a = (A2-A1)/(v2-v1)
    b = A1-a*v1
    Areturn = a*V+b
    return Areturn

def InterpolRho(V, rho, v):
    # function to interpol Rho, since altitude is function of Velocity
    if V < v[0]:
        return rho[0]
    
    elif V > v[-1]:
        return rho[-1]
    else:
        exitcondition = 1
        length_v = len(v)-1
        i = 0
        while exitcondition:
           
            if V == v[i]:
                rhoreturn = rho[i]
                exitcondition = 0
            
            elif V > v[i] and V < v[i+1]:
                rhoreturn = Interpol(V, rho[i], rho[i+1], v[i], v[i+1])
                exitcondition = 0  # exit
            
            else:
                i = i+1
                
            if i == length_v:  # security to exit the while
                print("AeroForces : Error in interpolating rho, returning 0")
                rhoreturn = 0
                exitcondition = 0
    
    return rhoreturn

def CoefInterpol( V, A, v):
    # A, function takes a numpy array composed of all matrices A [A1; A2; ...], all array types!!
    # v, an array of corresponding velocity
    # V, the velocity at which to compute the coef
    # size of each matrix 6 row, m column
    row = 6
    nmatrix = len(A[:, 1])/6
    if nmatrix == float:
        # error
        print('Aero forces: general coef matrix not a multiple of 6')
        return 0
    
    elif V < v[0]:
        # ill posed problem, the velocity is below smallest velocity for coef
        # use first matrix but print warning
        print("WARNING : velocity, V = {0:0.2f} is below first ref velocity for coef v = {1:0.2f}".format(V, v[0]))
        print("Continue with first matrix")
        return A[0:row, :]
    
    elif V > v[-1]:
        # same idea, the velocity is greater than max used to determine coef. Results in flight faster than cruise
        # use last matrix but print warning
        print("WARNING : velocity, V = {0:0.2f} is higher than last ref velocity for coef v = {1:0.2f}".format(V, v[-1]))
        print("Continue with last matrix")
        return A[-6:]
        
    else : #otherwise interpolate
        exitcondition = 1
        length_v = len(v)-1
        i = 0
        while exitcondition:
           
            if V == v[i]:
                Areturn = A[i*row:i*row+row]
                exitcondition = 0
            
            elif V > v[i] and V < v[i+1]:
                Areturn = Interpol(V, A[i*row:i*row+row, :], A[(i+1)*row:(i+1)*row+row, :], v[i], v[i+1])
                exitcondition = 0  # exit
            
            else:
                i = i+1
                
            if i == length_v+1:  # security to exit the while
                print("!!! FAILURE !!! AeroForces : Error in interpolation, returning 0")
                Areturn = 0
                exitcondition = 0

    return Areturn


def CalcForce_aeroframe_DEP(V, CoefMatrix, x, Tc, atmo, g, PropWing):
    """ Function to compute aerodynamic forces in the velocity frame (aero frame)
    for the DEP configuration. The propulsion force and moments are not computed here
    Since V is fixed, Coef Matrix must be calculated before
    Can handle DEP, with and without rudder, 2 or more engines
    """
    rho = atmo[1]
    a_sound = atmo[0]
    beta = x[1]
    p = x[2]
    r = x[4]

    #Compute aero forces
    # here x must be of the form (alpha, beta, p, q, r, da, de, dr) 
    # set non dim for p,q,r
    nonDim = np.ones(len(x))
    nonDim[2] = g.b/(2*V)
    nonDim[3] = g.c/(2*V)
    nonDim[4] = g.b/(2*V)
    x = x*nonDim


    #If prop wing, aileron contributions are included in patterson; Cl_delta_a ; Cm_delta_a ; Cn_delta_a
    if g.IsPropWing and g.isail:
        CoefMatrix[3:6, 5] = np.zeros(3)
        dail = x[5]
    else:
        dail = 0

    #    F=np.dot(CoefMatrix,x[0:7]) # commented form, modification to account for symmetric drag increase of side slip
    F = np.zeros(3)
    M = np.zeros(3)
    xsym = np.copy(x)


    xsym[1] = abs(xsym[1])  # make beta always positive since derivatives have already correct sign for drag and lift only
    xsym[5] = abs(xsym[5])  # make ailerons deflection always positive for drag increase and lift decrease

    if g.nofin == False:
              xsym[-1] = abs(xsym[-1])  # make rudder deflection always positive for drag increase and lift decrease


    F[0] = np.dot(CoefMatrix[0, 1:], xsym[1:])               # (                CD_BETA*|BETA| + CD_P*P^ +  CD_Q*Q^ + CD_R*R^ + CD_DA*|DA| +  CD_DE*DE   + CD_DR*|DR|)  !not alpha
    F[1] = np.dot(CoefMatrix[1], x)                          # ( CY_ALFA*ALFA + CY_BETA* BETA  + CY_P*P^ +  CY_Q*Q^ + CY_R*R^ + CY_DA*DA   +  CY_DE*DE   + CY_DR*DR)
    F[2] = np.dot(CoefMatrix[2], xsym)                       # ( CL_ALFA*ALFA + CL_BETA*|BETA| + CL_P*P^ +  CL_Q*Q^ + CL_R*R^ + CL_DA*|DA| +  CL_DE*DE + CD_DR*|DR|)  !calculated again later if interaction
    M = np.dot(CoefMatrix[3:6, :], x)

    DragQuad = F[0] + g.Cda*x[0]**2 + g.Cdb * x[0] + g.Cdc

    g.Cm_alpha_aero_interaction = Cm_alpha(V, CoefMatrix, x, Tc, atmo, g, PropWing)

    if g.IsPropWing:

        CoefMatrix[4,0] = g.Cm_alpha_aero_interaction                                                                   #Use new CM_ALPHA
        M = np.dot(CoefMatrix[3:6, :], x)                                                                               #Redo calculus with new Cm_alpha
        F[2] = np.dot(CoefMatrix[2][1:], xsym[1:]) + g.aht*x[0]-g.CL0_HT                                                #F[2] Calculated without alpha: CL_BETA*BETA + CL_P*P + CL_Q*Q + CL_R*R + CL_DA*|DA| + CL_DE*DE + CL_DR*|DR|   )
                                                                                                                        #Terms for horizontal tail added (alpha, and 0-alpha term) to modify x[0] and g.CL0_HT to take into account slisptream
        if V <= g.VelFlap:

            CLCl = PropWing.CalcCoef(Tc, V/a_sound, atmo, x[0], dail, g.FlapDefl, g, beta, p, V, r)


            if len(CLCl)>2 and g.IsPropWingDrag:
                # Drag is computed by patterson, add contribution of other variables (than alpha and dx)
                Fbody = np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T, F[1], -F[2]-CLCl[0]]) # add alpha=0 coefficients
                Moment = M+np.array([CLCl[1], g.Cm0 + g.Cm0_fl, CLCl[4]])


            else:
                Fbody = np.array([-DragQuad,F[1],-F[2]-CLCl[0]])  # add alpha=0 coefficients
                # add roll effect
                Moment = M+np.array([CLCl[1], g.Cm0 + g.Cm0_fl, CLCl[4]])

        else:
            CLCl = PropWing.CalcCoef(Tc, V/a_sound, atmo, x[0], dail, 0, g, beta, p, V, r)
            # by default lift and drag are computed here
            Fbody = np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T, F[1], -F[2]-CLCl[0]])  # add alpha=0 coefficients
            # add roll effect
            Moment = M+np.array([CLCl[1], g.Cm0, CLCl[4]])
    else:

        if V <= g.VelFlap:
            Fbody = np.array([-DragQuad, F[1], -F[2]-g.CL0 - g.CL0_fl])
            Moment = M+np.array([0, g.Cm0 + g.Cm0_fl, 0])
        else:
            Fbody = np.array([-DragQuad, F[1], -F[2]-g.CL0])  # add alpha=0 coefficients
            Moment = M+np.array([0, g.Cm0, 0])
 
    g.TotalDrag = abs(Fbody[0])
    g.Lift = abs(Fbody[2])
    g.lift2drag = abs(Fbody[2]/Fbody[0])
    Fbody = 0.5*V**2.0*rho*g.S*Fbody
    Moment = 0.5*V**2.0*rho*g.S*Moment*np.array([g.b, g.c, g.b])
    
    # Both are in the serodynamic frame
    return np.append(Fbody, Moment)

def Cm_alpha(V, CoefMatrix, x, Tc, atmo, g, PropWing):
    
    # From OpenVSP

    rho = atmo[1]
    a_sound = atmo[0]
    beta = x[1]
    p = x[2]
    q = x[3]
    r = x[4]
    alpha = x[0]
    CL_alpha_no_int = CoefMatrix[2, 0]


    Fx_vec = Tc * (2*rho*g.Sp*V**2)
    Fx = np.sum(Fx_vec)


    # here x must be of the form (alpha, beta, p, q, r, da, de, dr)
    # set non dim for p,q,r
    nonDim = np.ones(len(x))
    nonDim[2] = g.b/(2*V)
    nonDim[3] = g.c/(2*V)
    nonDim[4] = g.b/(2*V)
    x = x*nonDim


    if g.IsPropWing and g.isail:
        CoefMatrix[3:6, 5] = np.zeros(3)
        dail = x[5]
    else:
        dail = 0

    #Additional Parameters.

    #Cl_alpha with interaction calculus
    alpha_1 = 0 * np.pi/180
    alpha_2 = 2 * np.pi/180
    CL_alpha_interaction = ((PropWing.CalcCoef(Tc, V/a_sound, atmo, alpha_2, dail, g.FlapDefl, g, beta, p, V, r)[0] + g.aht*alpha_2) - (PropWing.CalcCoef(Tc, V/a_sound, atmo, alpha_1, dail, g.FlapDefl, g,beta, p, V, r)[0] + g.aht*alpha_1)) / (alpha_2-alpha_1)

    g.Vef2Vinf_2 = PropWing.Augmented_velocity_wing(Tc, V/a_sound, atmo, x[0], dail, g.FlapDefl, g,beta,p,V,r)          #(V_ef/V_inf)^2  (is not vtail/Vinf , keep that in mind)

    g.eta_downwash = 0.858675
    g.X_CA_wb = g.x_cg/g.c - ((CoefMatrix[4, 0] + g.aht*g.Hor_tail_coef_vol * 0.8) / (CoefMatrix[2, 0]-g.aht))   #--> Calculated so that downwash *ratio of dynamic pressures is 0.8
    #deberia ser 11.242 + 0.25 * 2.3244 = 11.82 , para que eso pase g.x_cg debe ser 11.52 , pero claro, entonces cambiar cm-alfa , debe ser un proceso iterativo.

    g.SM_CA_wingfuselage = (CoefMatrix[4, 0] + g.aht*g.Hor_tail_coef_vol * g.eta_downwash) / (CoefMatrix[2, 0]-g.aht)

    #it bothers me that if you calculate the CA_wingbody does not really match with geometry,
    #maybe because center of gravity is not well stimated? IN ATR for having a good value for
    #CA_wingbody, center of gravity should be around 11.4 m, is 12.41 now. Anyway calculus of
    # Cm alpha only cares about static margin

    #CALCULUS OF NEW CM_ALPHA_AERODYNAMIC

    Cm_alpha_interaction = (1 + ((CL_alpha_interaction - CoefMatrix[2, 0]) * g.SM_CA_wingfuselage)/CoefMatrix[4, 0]) * CoefMatrix[4, 0]
    #this formula is valid supossing that downwash, c.gravity and tail-wing pressure ratio does not change when implementing DEP

    # Slipstream velocity to free stream velocity

    VarVtoV = (1+Fx/(0.5*g.N_eng*rho*g.Sp*V**2))**0.5 - 1  # Same value than before with Patterson

    # Contracted slisptream diameter

    D_s = g.Dp * ((V + 0.5 * V * VarVtoV)/(V + V * VarVtoV)) ** 0.5

    # Epsilon without inflow effects

    eps = g.eps0 + (g.deps_dalpha / CL_alpha_no_int) * CL_alpha_interaction * alpha

    # H calculus :   Vertical Distance of Slipstream Center Line to Horizontal Tail

    h = g.z_h_w - g.lh*np.sin(alpha) - (g.xp + 0.25 * g.c)*np.sin(alpha) + g.lh2*np.sin(g.K_e * eps) + g.FlChord * g.c * np.sin(g.FlapDefl) + 0.25 * (g.xp + 0.25 * g.c) * np.sin(PropWing.alpha0_fl * g.FlapDefl)

    # PropWing.alpha0_fl has to be the change in airfoil section zero-lift angle-of-attack due to flap deflection (obtained using AVL)
    # PropWing.alpha0_fl is the change in alpha_0 for every radian deflection of flap. It is necessary to multiply by flaps deflection, in radians.
    # g.FlapDefl is in radians

    # Epsilon with inflow effects

    if h < 1.25 * D_s:
        extra_eps = (g.var_eps * VarVtoV) * np.pi / 180
        eps = eps + extra_eps

    # Dynamic pressure ratio in the horizontal tail

    V2 = (1 + VarVtoV) * V

    if (1 - (2*h / D_s)**2) > 0 :
        bs = D_s * (1 - (2*h / D_s)**2) ** 0.5
        Sh_s = 2 * bs * g.c_ht

        dpratio = ((Sh_s / g.Sh) * (1 + VarVtoV)**2 + (1-Sh_s/g.Sh))
    else:
        dpratio = (1+VarVtoV)**2

    #BASICAMENTE PARECE QUE LA PARTE MOJADA DE LA COLA HORIZONTAL POR EL TUBO DE SLISPTREAM ES 0. ELLO IMPLICARIA QUE LA PRESION DINAMICA NO CAMBIA ?
    #MIRAR EN EL PAPER; VAYA MIERDA; EL SLIPSTREAM TAMPOCO CAMBIA O QUE...? EN DECOL SI AFECTARA QUAND MEME

    # TAIL-OFF PITCHING MOMENT

    """
    Cm_alpha = 
    Cm_beta = 
    Cm_p = 
    Cm_q = 
    Cm_r =
    Cm_de =
    Cm_da =
    Cm_dr =
    """

    #Cm_0 = Cm0 + Cm_alpha*alpha + Cm_beta*beta + Cm_p*p +Cm_q*q + Cm_r*r

    c_flaps = g.c * np.sqrt((1+g.FlChord*(np.cos(g.FlapDefl)-1))**2 + (g.FlChord*np.sin(g.FlapDefl))**2)

    Cm_s_0 = g.N_eng * ((D_s * g.c)/g.S) * g.cm_0_s * ((g.Vef2Vinf_2*V/V)**2-1)

    Cm_s_df = (c_flaps/g.c)*(-0.25+0.32*g.FlChord / c_flaps) * (1+0.2*(1-np.sqrt(2) * np.sin(g.FlapDefl))) * g.CL0_fl


    if g.FlapDefl <= 30*np.pi/180:

       F = (0.5 * (g.FlapDefl*180/np.pi) / 30 + 0.25)*(c_flaps/g.c-1) + 0.05 * (g.FlapDefl*180/np.pi) / 30 + (g.Var_xac_fus / g.c) * (1-(g.FlapDefl*180/np.pi) / 30)
    else:

       F = -0.75 * (c_flaps/g.c-1) - 0.05

    Cm_s_alpha = F * g.CL0_fl

    Cm_tail_off = Cm_s_0 + Cm_s_df + Cm_s_alpha #+Cm_0

    # TAIL MOMENT
    Cm_tail = -(alpha + g.it - eps) * (g.aht2 * g.S/g.Sh) * dpratio * (g.Sh * g.lv)/(g.S * g.c)

    return Cm_alpha_interaction


