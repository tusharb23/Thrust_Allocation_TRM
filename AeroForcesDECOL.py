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
        CoefMatrix[3:6,5] = np.zeros(3)
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
              xsym[-1] = abs(xsym[-1]) # make rudder deflection always positive for drag increase and lift decrease


    F[0] = np.dot(CoefMatrix[0, 1:], xsym[1:])               # (                CD_BETA*|BETA| + CD_P*P^ +  CD_Q*Q^ + CD_R*R^ + CD_DA*|DA| +  CD_DE*DE   + CD_DR*|DR|)  !not alpha
    F[1] = np.dot(CoefMatrix[1], x)                          # ( CY_ALFA*ALFA + CY_BETA* BETA  + CY_P*P^ +  CY_Q*Q^ + CY_R*R^ + CY_DA*DA   +  CY_DE*DE   + CY_DR*DR)
    F[2] = np.dot(CoefMatrix[2], xsym)                       # ( CL_ALFA*ALFA + CL_BETA*|BETA| + CL_P*P^ +  CL_Q*Q^ + CL_R*R^ + CL_DA*|DA| +  CL_DE*DE + CD_DR*|DR|)  !calculated again later if interaction
    M = np.dot(CoefMatrix[3:6, :], x)

    DragQuad = F[0] + g.Cda*x[0]**2 + g.Cdb * x[0] + g.Cdc

    g.Cm_alpha_aero_interaction = Cm_alpha(V, CoefMatrix, x, Tc, atmo, g, PropWing)

    if g.IsPropWing:

        CoefMatrix[4,0] = g.Cm_alpha_aero_interaction                                                                   #Use new CM_ALPHA
        M = np.dot(CoefMatrix[3:6,:],x)                                                                                 #Redo calculus with new Cm_alpha
        F[2] = np.dot(CoefMatrix[2][1:],xsym[1:]) + g.aht*x[0]-g.CL0_HT                                                 #F[2] Calculated without alpha: CL_BETA*BETA + CL_P*P + CL_Q*Q + CL_R*R + CL_DA*|DA| + CL_DE*DE + CL_DR*|DR|   )
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
            Fbody = np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T, F[1], -F[2]-CLCl[0]]) # add alpha=0 coefficients
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
    Moment = 0.5*V**2.0*rho*g.S*Moment*np.array([g.b,g.c,g.b])

    return np.append(Fbody, Moment)





def Cm_alpha(V, CoefMatrix, x, Tc, atmo, g, PropWing):
    """ Function to compute longitudinal stability issues regarding the influence on the horizontal tail
    when there is interaction. To compute for the slipstream velocity and down wash and see the influence on
    the horizontal tail.

    Normally
    Cm_aero = Cm0 + Cm_alpha * alpha + Cm_q * q + Cm_de * de + Cm_u * u + Cm_d(alfa) * d(alfa) + Cm_d(de) * d(de)

    Propulsion moments are already accounted in the equation modulus.

    Cm_u = 0 as M is regime is always incompressible

    Cm_q is calculated in OpenVSP, more about how to calculate it in Page 345

    Cm_delta_e Page 345

    Cm_de

    Cm_d(de)   Page 345

    Cm_d(alfa) ----> Teoría del retardo de estela




    Para proceder vamos a:

    Calcular el momento aerodinamico alrededor del centro aerodinamico conjunto ala-fuselaje. Este momento NO depende
    del angulo de ataque. Simulacion de OpenVSP con el ala y el fuselaje a 0 grados.

    Calcular el momento aerodinamico alrededor del centro aerodinamico de la cola horizontal. Este momento NO depende
    del angulo de ataque. Simulacion de OpenVSP con la cola horizontal a 0 grados.



    Luego habra que calcular momentos de la sustentacion y la resistencia en el conjunto ala fuselaje
        Para ello:
                    Sustentacion en el cjto ala fuselaje ------>  (Lo tienes)

                    Resistencia en el cjto ala fuselaje -------> (lo tienes aunque hay que adaptarlo para que no
                    cuente solo la del aumento de presion dinamica en la parasita)

                    Centro aerodinamico en el cjto ala fuselaje. El fuselaje casi no sustenta, asi que se puede
                    considerar que sera el mismo que el del ala. Los tienes calculados, sino de todas formas OpenVSP deberia
                    dartelo

                    Calcular la sustentacion y la resistencia en el estabilizador horizontal.
                            Para ello : Necesitas un CL0 -----> de OpenVSP
                                        Necesitas un Cl_alfa -----> de OpenVSP
                                        Necesitas un angulo de ataque
                                        Necesitas una velocidad en la cola

                    Centro aerodinamico del estabilizador horizontal

                    Teoria para el calculo de la deflexion de estela

                    Teoria para el calculo de la presion dinamica en la cola.


    Y después calculamos directamente todo el momento aerodinamico EN EJES CUERPO Y LO CAMBIAMOS A EJES VELOCIDAD!! como:

    M = Cmac_wb + Cmac_HT + cg-ca ^ CL + cg-ca ^ CD + cg-caHT ^ CL_HT + cg-caHT ^ CD_HT = 0

        Si descompones este momento, habra elementos que:

           Dependan de alfa: Formaran Cm_alpha
           Dependan de nada: formaran Cm0
           Dependan de delta_e: Lo incluimos
           Dependeran de q:  Lo incluimos?
           Dependeran de d(alfa): Lo incluimos?  Con toeria del retardo




        cg-ca = SM

        cg-caHT = g.lv

        za =

        L_HT =    ( -0.0068 + 0.7798 *   ALPHA_HT   )  * (PRESION DINAMICA COLA)  /  (PRESION DINAMICA AVION)

        D_HT =    ( (POLAR) * ALPHA _HT  )   * (PRESION DINAMICA COLA)  /  (PRESION DINAMICA AVION)
        
    """

    # From OpenVSP


    rho = atmo[1]
    a_sound = atmo[0]
    beta = x[1]
    p = x[2]
    r = x[4]

    # here x must be of the form (alpha, beta, p, q, r, da, de, dr)
    # set non dim for p,q,r
    nonDim = np.ones(len(x))
    nonDim[2] = g.b/(2*V)
    nonDim[3] = g.c/(2*V)
    nonDim[4] = g.b/(2*V)
    x = x*nonDim


    if g.IsPropWing and g.isail:
        CoefMatrix[3:6,5] = np.zeros(3)
        dail = x[5]
    else:
        dail = 0


    #Additional Parameters.



    g.Vef2Vinf_2 = PropWing.Augmented_velocity_wing(Tc, V/a_sound, atmo, x[0], dail, g.FlapDefl, g,beta,p,V,r)          #(V_ef/V_inf)^2  (is not tail/Vinf , keep that in mind)

    g.eta_downwash = 0.8
    # g.X_CA_wb = g.x_cg/g.c - (      (CoefMatrix[4,0] + g.aht*g.Hor_tail_coef_vol *0.8)/ (CoefMatrix[2,0]-g.aht)   )   #--> Calculated so that downwash *ratio of dynamic pressures is 0.8



    g.SM_CA_wingfuselage =  (CoefMatrix[4,0] + g.aht*g.Hor_tail_coef_vol *g.eta_downwash)/ (CoefMatrix[2,0]-g.aht)

    #it bothers me that if you calculate the CA_wingbody does not really match with geometry,
    #maybe because center of gravity is not well stimated? IN ATR for having a good value for
    #CA_wingbody, center of gravity should be around 11.4 m, is 12.41 now. Anyway calculus of
    # Cm alpha only cares about static margin


    #Cl_alpha with interaction calculus
    alpha_1 = 0 * np.pi/180
    alpha_2 = 2 * np.pi/180
    CL_alpha_interaction = ((PropWing.CalcCoef(Tc, V/a_sound, atmo, alpha_2, dail, g.FlapDefl, g, beta, p, V, r)[0] + g.aht*alpha_2) - (PropWing.CalcCoef(Tc, V/a_sound, atmo, alpha_1, dail, g.FlapDefl, g,beta, p, V, r)[0] + g.aht*alpha_1)) / (alpha_2-alpha_1)


    #CALCULUS OF NEW CM_ALPHA_AERODYNAMIC

    Cm_alpha_interaction = (1 + ((CL_alpha_interaction - CoefMatrix[2,0]) * g.SM_CA_wingfuselage)/CoefMatrix[4,0]) * CoefMatrix[4,0]
    #this formula is valid supossing that downwash, c.gravity and tail-wing pressure ratio does not change when implementing DEP





    return Cm_alpha_interaction


