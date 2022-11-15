""" The DECOLgeometry.py models the model DECOL aircraft defining all its geometrical parameters. """
""" - Modified Eric Nguyen Van and David Planas Andres 's work"""
""" For the pupose of Optimal Thrust Allocation """
""" - Modified by Saumya Sarawagi"""

import math
import numpy as np
from Propeller_Data import Prop


""" The Thrust model is based on the momentum theory and the actuator disc theory.
    It uses data from SeligProp86 to obtain the values of the Thrust and Power co-efficients.
    The model interpolates values from the SeligProp86
    TO DO:
        i Analyse non-linearities
        ii Take care of Max_Power
        iii Get the values of J from deltx and find corresponding power and thrust - DONE in SeligPropellerRead"""
        
""" Advance Ratio (J) = Velocity of flow at the propeller / Diameter(D) * Rotational Speed (n)
    n = n_max*dx, where n_max = 10000RPM"""
    

class data:
    
    # Mass & Moments of Inertia
    x_cg = 0.713 # Location of center of gravity (m)
    m = 8.25 # kg 
    Ix = 1.1 # kg/m^2
    Iy = 1.2 # kg/m^2
    Iz = 2.0 # kg/m^2
    Ixz = 0 # kg/m^2
    
    # Geometrical parameters of the DECOL aircraft
    S = 0.5 # Wing area (m^2)
    b = 2 # Wing Span (m)
    c = 0.25 # Wing Cord (m)
    AR = 8 # Wing Aspect Ratio
    wingsweep = 0                                                                                               
    dihedral = 0
    FusWidth = 0.170 # Width of the fuselage
    
    # Vertical Tail Design
    zf = 0.165+0.085 # 3.16m, z position of the MAC of the fin
    fswept = 0/180*math.pi # Swept angle of vertical tail
    ftaper = 1.0 # Taper Ratio
    fAR = 1.8 # Aspect ratio    
    lv = 1.809-x_cg  # m Distance from center of gravity to center of pressure of horizontal tail
    lemac = 0.610124  # Distance from the tip to the leading edge of the MAC (here MAC matches with root chord)
    VolH = 0.8
    bh = 0.611                # Horizontal tail wingspan
    Sh = 0.0877  # Horizontal tail surface
    Hor_tail_coef_vol = (Sh*lv) / (S*c)     # 0.7552  Volume coef of Horizontal tail
    taudr = 0.24  # Ruder efficiency factor see nicolosi paper and dela-vecchia thesis
    it = 2.1/180*np.pi             # Angle between zero lift line and the fuselage or horizontal tail tilt angle 
    Var_xac_fus = -0.019914 # Variation in the aerodynamic centre. Compute by difference of the neutral points on OpenVSP between wing and wing + fuselaje (without hor tail)                                                                                            

    # Flap & Aileron definition
    isflap = True                                                                                                       
    FlPosi = 0.0425 # With respect to wingspan, the start position of flap [0,0.5]
    FlRatio = 0.2325*2 # The total flap length to wingspan ratio
    FlSpan = 0.450 # Flap Span (m)
    FlChord = 0.25 # Flap Cord (% of c)
    isail = True                                                                                                        
    AilDiff = 0.5 # Ratio of aileron differential: deflection down is only half as important
    AilPosi = 0.275 # [0,0.5]
    AilRatio = 0.225 # One aileron should be consistent with b/2 >= AilRatio + AilPosi
    AilSpan = 0.450 # Aileron Span (m)
    AilChord = 0.25 # Same as flap
    
    # Power
    P_a = 576
    hp = 0 # Rotor term
    N_eng = 8  # Number of engines
    TipClearance = True
    dprop = 0.1 # Spacing between the propellers
    dfus = 0
    prop_eff = 0.7 # Propulsion efficiency - considering it to be constant -  NOT IN REALITY
    ip = (-2.41)/180*np.pi # Propeller incidence angle with respect to zero lift line
    h_m = 0.05 # Distance from leading edge to propeller. Propeller is forward      
    x_m = 0.104+h_m # Distance of the leading edge from the CG + h_m      
    z_m = 0.075  # IGUAL ES 0.075 approximate vertical distance between engine and CG. Propellers is above. Computed with OpenVSP
    z_h_w = 0.025  # vertical distance from the horizontal tail 1/4 point to the propeller axis. Computed with OpenVSP
    lh = 1.16575    # Horizontal distance between the aerodynamic centers of horizontal tail and wing (0.25 of their chord in root is enough) Computed with OpenVSP.
    lh2 = 0.942  # Horizontal distance from the wing trailing edge to the horizontal tail leading edge. Computed with OpenVSP
    K_e = 1.44   # Down wash factor, see Modeling the Propeller Slipstream Effect on Lift and Pitching Moment, Bouquet, Thijs; Vos, Roelof
    c_ht = 0.145  # Average chord of the horizontal tail
    var_eps = 1.5  # parameter for inflow in slisptream. See Modeling the Propeller Slipstream Effect on Lift and Pitching Moment, Bouquet, Thijs; Vos, Roelof
    cm_0_s = -0.0512  # zero lift pitching moment of the wing section at the propeller axis location. From the xlfr5 file, alpha = 0°  

    # Propeller Data for Accurate Model
    p_max = 10*14.8; # To be verified
    AvGearEff=0.80 # average effiency of engine+gearbox
    D = 8*0.0254 # Diameter in m
    Pitch = 6*0.0254 # Pitch in m
    n_max = 16280/60 # Maximum rotations per second considering a nominal voltage of 14.8V
    

    # Propeller Wing Interaction
    IsPropWing=True
    IsPropWingDrag=True
    nofin = False
    
    # ---Unique coeff ---
    aht = 0.5448
    aht2 = 0.7049       # Horizontal tail lift coefficient, for the tail analysed alone. Dimensioned with S.
    #  Cm_de = -2 # per rad, is constant for DECOL             You can use the one from STAB file, or this one
    Cm_alpha_fus = 0.015*180/np.pi

    # alpha=0 coeff

    # without flaps
    CD0T = 0.0636         # global one  extracted from flight not stab the file
    CDO_wo_VT = 0.0627
    CL0 = 0.44768
    CL0_HT = -0.002489
    Cm0 = 0.004684
    CDMax90 = 1.25

    # Drag polar without patterson. Interpolated from VSP v26, updated VSPAERO
    Cda_fl_0 = 1.6372
    Cdb_fl_0 = 0.2934
    Cdc_fl_0 = 0.0637

    # with flaps down 15°
    Cd0_fl_15 = 0.026699          # Extra drag due to flap for prop wing interaction if no patterson used
    CL0_fl_15 = 0.27142           # Extra lift due to flaps when no patterson used
    Cm0_fl_15 = 0.096931          # Extra pitch moment due to flaps

    Cda_fl_15 = 1.4937      # Coefficients for calculus of CD (CD=Cda * alpha ** 2 + Cdb * alpha + Cdc) for flaps = 15 °
    Cdb_fl_15 = 0.4327
    Cdc_fl_15 = 0.0905

    # Down-Wash parameters

    # No flaps
    eps0_flaps0 = 2.274 * np.pi/180   # Downwash at 0 angle of attack in no flaps configuration
    deps_dalpha_flaps0 = 0.281        # Derivative of downwash with respect to alpha, no flaps conf

    # 15° flaps
    eps0_flaps15 = 4.07 * np.pi/180   # Downwash at 0 angle of attack in 15° flaps configuration
    deps_dalpha_flaps15 = 0.281       # Derivative of downwash with respect to alpha, 15° flaps conf

    #Airfoil characteristics
    Cd0_laminar = 0.01252 
    Cd0_turbulent = 0.01252


    #wing tilt angle, angle between reference line of fuselaje and reference line of profile
    alpha_i = 3.2/180*np.pi
    
    #airfoil zero lift angle       angle between reference line of profile and zero lift line of airfoil. Negative means that airfoil lifts with 0 local angle of attack
    alpha_0 = -2.41/180*np.pi  # update with vsp file

    
    # Input file name
    Files = ['cldistribution', 'polar', 'flappolar', 'aileronpolar']  # best to replace the value
    alphaVSP = 0/180*np.pi
    PolarFlDeflDeg = 5


    #unique data go here:
    def CalcKf(self, bv, r):
        return 1.4685*(bv/(2*r))**(-0.143)
    
    def CalcKw(self, zw, rf):
        return -0.0131*(zw/rf)**2-0.0459*(zw/rf)+1.0026
        
    def CalcKh(self, zh, bvl, Sh, Sv):
        x = zh/bvl
        Khp = 0.906*x**2-0.8403*x+1.1475
        Khs = math.exp(0.0797*math.log(Sh/Sv)-0.0079)
        return 1+Khs*(Khp-1)
    
    def CalcKdr(self, Kf, Av):
        # Kdr=(1+((Kf-1)/2.2))*(1.33-0.09*Av) # for T-tail formula
        Kdr = 1.07*(1+(Kf-1)/2.2)  # for a body mounted H tail
        return Kdr
            
    def SetEngineNumber(self):
        # To apply Patterson Theory, wing tip clearance should be True
        self.Dp = (self.b/2-self.FusWidth/2)/(self.N_eng/2+self.dfus+(self.N_eng/2-1)*self.dprop)
        self.Sp = self.Dp**2/4*math.pi
        self.xp = self.Dp/2
        self.step_y=self.Dp+self.dprop*self.Dp
        self.PosiEng = np.arange(self.FusWidth/2+self.Dp*(self.dfus+0.5),self.b/2,self.step_y)
        self.PosiEng = np.append(-self.PosiEng,self.PosiEng)
        self.PosiEng = np.sort(self.PosiEng)
        return   
    

  
    def __init__(self, inop_eng, bv=0.329, r=0.057/2, zw=0.073, rf=0.085, zh=0, bvl=0.348, Sh=0.0877, Sv=0.06, TipClearance = True, dfus=0, dprop=0.1):
       
        self.Propeller = Prop(self.D, self.Pitch, self.n_max, self.p_max)
        self.VTsize=1 # We do not modify the vertical taile size
        self.SetEngineNumber()
        self.inop = inop_eng
        self.Sv = Sv
        self.SvBase = Sv
        self.bv = bv
        self.r = r
        self.Av = bv**2/Sv
        self.Sh = Sh
        self.zw = zw
        self.rf = rf
        self.zh = zh
        self.bvl = bv+r
        
        #Update drag
        self.Cdvt = 0.01374*Sv/self.S
        
        #Nicolosi csts
        self.Kf = self.CalcKf(bv, r)
        self.Kw = self.CalcKw(zw, rf)
        self.Kh = self.CalcKh(zh, bvl, Sh, Sv)
        self.Kdr = self.CalcKdr(self.Kf, self.Av)
        self.taudrBase = 0.018/Sv*0.8  #Rudder surface/fin surface * efficiency of rudder, rudder surface proportional to vtsize


        #Aero coefficients for DragQuad
        self.FlapDefl = 0   # For our problem statement we consider flap and aileron deflection to be zero
        self.Cdo_fl, self.CL0_fl, self.Cm0_fl, self.Cda, self.Cdb, self.Cdc = self.AeroCoefs(self.FlapDefl)
        

    def AeroCoefs(self, FlapDefl):
         if FlapDefl == 0:
             self.Cd0_fl = 0
             self.CL0_fl = 0
             self.Cm0_fl = 0
             self.Cda = self.Cda_fl_0
             self.Cdb = self.Cdb_fl_0
             self.Cdc = self.Cdc_fl_0
             self.eps0 = self.eps0_flaps0
             self.deps_dalpha = self.deps_dalpha_flaps0
         elif FlapDefl == 15:
             self.Cd0_fl = self.Cd0_fl_15
             self.CL0_fl = self.CL0_fl_15
             self.Cm0_fl = self.Cm0_fl_15
             self.Cda = self.Cda_fl_15
             self.Cdb = self.Cdb_fl_15
             self.Cdc = self.Cdc_fl_15
             self.eps0 = self.eps0_flaps15
             self.deps_dalpha = self.deps_dalpha_flaps15
         else:
             print("Chose an allowable value for flaps deflection, options are: No flaps (0°) or 15°")
             exit()
         return self.Cd0_fl, self.CL0_fl, self.Cm0_fl, self.Cda, self.Cdb, self.Cdc
        
        
    def NicolosiCoef(self, MCoef, Mach):
        # function to compute Cy, Cl and Cn derivatives using VeDSC methods
        # replaces the coefficients in the matrix Mcoef by those computed by VeDSC
        # to call only once, it computes for every velocities/Mach number. Then linear interpol
        # Cy_beta, Cy_p = 0, Cy_r = 0, Cl_beta = 0, Cl_r = 0, Cn_beta= 0, Cn_p = 0, Cn_n = 0
        
        MVeDSC = np.copy(MCoef)  # create a copy of the coefficients matrix
        if self.nofin == False:
            # add a column to account for rudder
            dimZero = len(MVeDSC[:, 0])
            MVeDSC = np.hstack((MVeDSC, np.zeros((dimZero, 1))))
        
        K = self.Kf*self.Kh*self.Kw
        for i in range(len(Mach)):
            Av = self.bv**2/self.Sv
#            print(Av)
            cla = 2*math.pi # local lift curve slope coefficient of fin (perpendicular to leading edge) per rad
            eta = cla/(2*math.pi)
            av = -(Av * cla * math.cos(self.fswept) * 1/math.sqrt(1-Mach[i]**2*(math.cos(self.fswept))**2))/(Av*math.sqrt(1+4*eta**2/(Av/math.cos(self.fswept))**2)+2*eta*math.cos(self.fswept))# mach number formula
#            print(av)
            VeDSC_Coef = np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
            VeDSC_Coef[0, 0] = K*av*self.Sv/self.S
            VeDSC_Coef[0, 1] = K*av*self.Sv/self.S*self.zf/self.b*2
            VeDSC_Coef[0, 2] = -K*av*self.Sv/self.S*2*self.lv/self.b #apparently x position of fin doesn't change
            VeDSC_Coef[0, 3] = -self.Kdr*av*self.taudr*self.Sv/self.S
            VeDSC_Coef[1, 0] = K*av*self.Sv/self.S*2*self.zf/self.b
            VeDSC_Coef[1, 1] = 0
            VeDSC_Coef[1, 2] = -K*av*self.Sv/self.S*self.zf/self.b*self.lv/self.b*2.0
            VeDSC_Coef[1, 3] = -self.Kdr*av*self.taudr*self.Sv/self.S*2*self.zf/self.b
            VeDSC_Coef[2, 0] = -K*av*self.lv/self.b*self.Sv/self.S
            VeDSC_Coef[2, 1] = -K*av*self.lv/self.b*self.Sv/self.S*self.zf/self.b*2.0
            VeDSC_Coef[2, 2] = K*av*self.lv/self.b*self.Sv/self.S*self.lv/self.b*2.0
            VeDSC_Coef[2, 3] = self.Kdr*av*self.taudr*self.lv/self.b*self.Sv/self.S
            # Coefficients are computed now access the right matrix and replace them
            VarPosi = (1, 2, 4)
            EffPosi = (1, 3, 5)
            NumEff = 6  # number of force equations
#            print(VeDSC_Coef[2,2])
            for kk in range(len(EffPosi)):
                # Replace rudder coefficient
                if self.nofin == False:
                    MVeDSC[EffPosi[kk]+i*NumEff, -1] = VeDSC_Coef[kk, -1]
                # Now coefficients are from the finless DECOL. Add new coefficients to the matrix
                for jj in range(len(VarPosi)):
                    if VeDSC_Coef[kk, jj] != 0:
                        MVeDSC[EffPosi[kk]+i*NumEff, VarPosi[jj]] = MVeDSC[EffPosi[kk]+i*NumEff, VarPosi[jj]]+VeDSC_Coef[kk, jj]
#            print(VeDSC_Coef)
            
        return MVeDSC

    def Thrust(self, dx, V):
        
        self.Thr = self.Propeller.Thrust(dx, V)
        # Get thrust components along the 3 axis
        Thr_xb = 0
        Thr_yb = 0
        Thr_zb = 0
        Thr_xb = np.sum(self.Thr*math.cos(self.ip))
        Thr_zb = -np.sum(self.Thr*math.sin(self.ip)) # The direction of thrust will be along negative z-axis
        Thr_body = [Thr_xb, Thr_yb, Thr_zb]
        # The thrust is in the body frame so, it needs to be transformed to the aerodynamic frame in Equations
        return Thr_body  

    def Torque(self,dx,V):
        # Compute torque consdering ip and hence, moments along all 3 axis
        Moment= np.zeros((self.N_eng,3))
        for i in range(self.N_eng):
            a= np.array([ self.x_m , self.PosiEng[i] , self.z_m])
            b=np.array([ self.Thr[i]*np.cos(self.alpha_i + self.alpha_0+self.ip)  ,  0  ,  -self.Thr[i]*np.sin(self.alpha_i + self.alpha_0+self.ip)  ])
            Moment[i,:] = np.cross(a,b)
        Thrust_Moment = np.array(( np.sum(Moment[:,0]), np.sum(Moment[:,1]) , np.sum(Moment[:,2]) ) )  
        return Thrust_Moment
   
    def DragModel(self, alpha, alpha_patter):
        # Added drag to account for induced drag due to flap deflection
        self.Cd_fl=(alpha_patter-alpha)*self.Cdalpha * self.FlRatio
        return None