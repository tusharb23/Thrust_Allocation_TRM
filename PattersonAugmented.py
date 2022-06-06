# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:32:24 2018

Integrated prop-wing interaction model

Based on Patterson added the flap treatment

Edit 18.06.19: Adding a model for stall. The lift formula is:
    CL = CLslope * alpha if alpha < alpha_stall
    CL = CLslope * sin(alpha_stall) * cos(alpha)/cos(alpha_max) # from Jamesson
alpha_max must be given in aircraft class

Edit 06.06.19: Bugs fixed, integration for drag corrected and cdip, drag from 
    propeller wash removed. All computation are made with respect to local data.
    >> Propellers should be away from wingtip to estimate drag accurately <<

Edit 25/05/19: Rendering it completely non-dimensional:
    -Taking as input Tc = Thrust / (2*rho*plane.Sp*V**2)

Edit 24/07/2018 : Adding drag computation based on :
    -Induced drag formulation of potential flow theory Di=L*tan(w/V)
    -Friction drag increase due to turbulent transition accelerated 
    by propellers (taken into account by adding a baseline lift distribution 
    file with forced turbulent transition)

@author: Eric Nguyen Van
         david.planas-andres
"""

import numpy as np
import ReadFileUtils as Read  # utils to read Xfoil file
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

## ----- Surrogate model utility ------

class PropWing:
    '''


    This class defines a prop wing interation.
    
    It includes flap effect and differential thrust / non-uniform thrust distribution
    
    Uses information from class 'aircraft'
    
    The correct way to use it :
        -Instantiate it with a plane model and the a vector of file name containing:
            -Name of Cl distribution('.fem' from VSP) and stab file at diff Mach
            -polar file of airfoil without flap (Xfoil/XFLR5 style))
            -polar file of airfoil with flap, (Xfoil/XFLR5 style)
            -Don't forget to give aoa of VSP Cl distribution file plane.alphaVSP
            
        -Call :
            CalcCoef: gives CL, Cl, Cdi, Cd_wash and Cd0 (CD = Cdi+Cd_wash+Cd0)
            PlotDist: plot the lift distribution
            PatterJames: returns lift distribution under np.array(dtype('Yposi','Area','Chord','Cl'))
            
        Input to these function (dx, V, atmospher, aoa, dfl, Files, plane):
            
            dx : vector of engine power setting :
                -0<dx<1
                -size should be equal to number of engine
                
            V : flight velocity in m/s
            
            Atmosphere : vector of atmospher parameters:
                [sound velocity, air density] can be obtained with 'getAtmo' function of 'aircraft'. SI units
                
            aoa : angle of attack in rad
            
            dfl : flap deflection in ° (only symmetric) positive down
            
            Files : dictionnary containing the name of the files containing polars and cl distribution
            VLM cl distribution can be given at different mach number and the code will interpolate the result with "V" and "Atmospher" given.
            The compressibility is taken into account only in the VLM
            ['fem':["path+filenameMach1.fem","path+filenameMach2.fem",...],'AirfoilPolar':"path+filename.txt",'FlapPolar':"path+filename.txt"]
            
            plane : intense of 'aircraft' type class.
            
            
        
    Make sure to have the following in aircraft class:
        -thrust function
        -airfoil cd0 laminar(cst value)
        -airfoil cd0 turbulent(cst value)
        -alphaVSP : finite wing alpha_0 from VLM computation
        
    Error or differences : 
        -Small ofset in lift slope between VLM and patter even at CT=0. Should come from vsp aero
        -The formulation of Patterson introduces a divergence of the lift multiplier at L=0. A new formulation should be used
        -When adding 0.5*Cdip the drag still appears to have a better evolution/fitting mostly at high CL. Bad at lower than 0.75CL for T=1000N/m^2
            Could be some unconcidered effects (swirl recovery?) or just Patterson divergence.
        -Not really an error but the integration of wi still integrates around the discontinuity at 1/(y1-y) when y=y1
    
    Good points :
        -Increase of lift slope with Ct accurately captured (compared to vlm). Validates the model for lift increase (prop wash + delta V)
        -Drag precision is ok for low (less than 500N/m^2) prop loading
        
    Still doesn't take into account the wingtip propeller effects.
    
    '''
    # Self data, a few stuff to check the results
    RecomputeDrag = True
    alpha_ep = []
    aoa = 0
    SetPropWash = True
    LmFl = np.array([])
    beta = np.array([])
    Cd0_vec = np.array([])
    PlotDrag = False
    
    # definition of surrogate coefficients
    C0 = np.array([0.378269, 0.748135, -0.179986, -0.056464, -0.146746, -0.015255])
    C1 = np.array([3.071020, -1.769885, 0.436595, 0.148643, -0.989332, 0.197940])                                       #Patterson, page 94 (116)
    C2 = np.array([-2.827730, 2.054064, -0.467410, -0.277325, 0.698981, -0.008226])
    C3 = np.array([0.997936, -0.916118, 0.199829, 0.157810, -0.143368, -0.057385])
    C4 = np.array([-0.127645, 0.135543, -0.028919, -0.026546, 0.010470, 0.012221])
    
    #functions
    def __init__(self,plane,Files):
        # Will check the necessary data are given in aircraft
        print("PropWing interaction will use for friction drag, Cd0 laminaire : {0}, CD0 turbulent : {1}".format(plane.Cd0_laminar,plane.Cd0_turbulent));
        print("PropWing interaction will use zero lift angle : {0}".format(plane.alphaVSP))
        print("PropWing interaction will use propeller ip : {0}".format(plane.ip))
        
        #import the local lift distribution
        #if many files are given, assume a variation in Mach
        self.NumFiles = 1
        if len(Files['fem'])>1:
            print('Reading multiple files')
            self.NumFiles = len(Files['fem'])
            
            #Read first to have the format
            CLslope, AoAZero, Mach = Read.ReadSectionCLslope(Files['fem'][0])
            
            self.CLslope = np.zeros( (len(CLslope),len(CLslope[1,:]),self.NumFiles) )
            self.AoAZero = np.zeros( (len(CLslope),len(CLslope[1,:]),self.NumFiles) )
            self.M_vec = np.zeros((self.NumFiles))
            self.CLslope[:,:,0] = np.copy(CLslope)
            self.AoAZero[:,:,0] = np.copy(AoAZero)
            self.M_vec[0] = Mach
            for i in range(1,self.NumFiles):
                self.CLslope[:,:,i], self.AoAZero[:,:,i], self.M_vec[i]= Read.ReadSectionCLslope(Files['fem'][i])
            
            #Correction for any wing incidence angle in VSP
            self.AoAZero[:,-1,:]=self.AoAZero[:,-1,:] + plane.alpha_i
            
        else:
            self.CLslope, self.AoAZero, self.M_vec= Read.ReadSectionCLslope(Files['fem'][0])
            #Correction for any wing incidence angle in VSP
            self.AoAZero[:,-1]=self.AoAZero[:,-1] + plane.alpha_i
            
        #That's to manage airfoil drag after stall
        alphaDrag, self.StallDrag = Read.ReadAirfoilDrag(Files['AirfoilPolar'])
        self.alphaDrag = alphaDrag/180*np.pi
        self.StallDrag = interp1d(self.alphaDrag,self.StallDrag)
        
        # Read flap and aileron polars if any
        if plane.isflap == True:
            #assume no change for ailerons efficiency with respect to Mach number
            self.alpha0_fl = ((Read.ReadAlpha0(Files['FlapPolar']) - Read.ReadAlpha0(Files['AirfoilPolar']))/plane.PolarFlDeflDeg)/180*np.pi
        
        if plane.isail == True:
            self.alpha0_ail = (Read.ReadAlpha0(Files['AileronPolar']) - Read.ReadAlpha0(Files['AirfoilPolar']))/plane.PolarAilDeflDeg
            
      
            
    def Interpol(self,Input,M):
        
        if self.NumFiles<2:
            #No data for interpolation
            return np.copy(Input)
        
        BaseInput = np.copy(Input[:,:,0])
        MachInput = np.copy(Input)
        
        if M<=self.M_vec[0] :
            # use first coeff file
            return BaseInput
    
        elif M>=self.M_vec[-1]:
            #Use last coeff file
            BaseInput = np.copy(MachInput[:,:,-1])
            return BaseInput
        
        else:
            exitcondition=1
            length_v=len(self.M_vec)-1
            i=0
            while exitcondition :
               
                if M==self.M_vec[i]:
                    #if it's exactly on one file
                    BaseInput = np.copy(MachInput[:,:,i])
                    exitcondition=0
#                    print("Exactly equal")
#                    print(self.M_vec[i])
                
                elif M>self.M_vec[i] and M<self.M_vec[i+1]:
                    #linear interpolation
                    a=(MachInput[:,-1,i+1]-MachInput[:,-1,i])/(self.M_vec[i+1]-self.M_vec[i])
                    b=MachInput[:,-1,i]-a*self.M_vec[i]
                    Areturn=a*M+b
                    BaseInput[:,-1]=Areturn
                    exitcondition=0 #exit
                
                else:
                    i=i+1
                    
                if i==length_v: #security to exit the while
                    print("AeroForces : Error in interpolating dist Cl, returning dist at M=0")
                    exitcondition=0
        
        return BaseInput

    def BetaSurro(self,a, Mu, rho, SectMu):
        """
        This function computes the beta, corrective term of Patterson propeller
        lift model.
        It implements the surrogate model present in the paper "High lift prop
        system for nasa's sceptor"
        Input variables:
            a : aircraft ATR class, with automatic limited propeller distribution
            Mu : Vjet/V
            rho : actual air density
            SectMu : Vector saying which mu(or engine) each section is associated with. If 0 mu = 1 (no blowing)
        Outputs :
            beta : vector of beta value in the order of the deltax given
        """
        
        # definition of surrogate coefficients
#        C0 = np.array([0.378269, 0.748135, -0.179986, -0.056464, -0.146746, -0.015255])
#        C1 = np.array([3.071020, -1.769885, 0.436595, 0.148643, -0.989332, 0.197940])
#        C2 = np.array([-2.827730, 2.054064, -0.467410, -0.277325, 0.698981, -0.008226])
#        C3 = np.array([0.997936, -0.916118, 0.199829, 0.157810, -0.143368, -0.057385])
#        C4 = np.array([-0.127645, 0.135543, -0.028919, -0.026546, 0.010470, 0.012221])
        
        #Definition of surrogate vector
        Lratio = 0
        Rratio = 0
        
        #Retrieve the local chord:
        if self.NumFiles>1:
            LocalChord = self.CLslope[:,2,0]
        else:
            LocalChord = self.CLslope[:,2]
        
        beta=np.zeros(len(LocalChord))
        for i in range(len(beta)):
            Lratio = a.xp/LocalChord[i]                                                                                 #xp is distance between propeller and leading edge
            Rratio = a.Dp/(2*LocalChord[i])                                                                             #Dp is the radio of the propeller
            if SectMu[i] != 0:
                X = np.array([1, Lratio, Lratio**2, Lratio*Mu[int(SectMu[i])-1], Mu[int(SectMu[i])-1], Mu[int(SectMu[i])-1]**2])
            else:
                X = np.array([1, Lratio, Lratio**2, Lratio*1, 1, 1**2])
            # Compute the whole thing
            beta[i] = np.dot(self.C0,X) + np.dot(self.C1,X)*Rratio + np.dot(self.C2,X)*Rratio**2 + np.dot(self.C3,X)*Rratio**3 + np.dot(self.C4,X)*Rratio**4
        
        return beta

## ------- Utility re-organize the lift in custom nnp data ------------

    def ReOrganiseLift(self, lift):
        # reorganise lift distribution for plotting or other uses
        dtype=[('Yposi', np.float), ('Area', np.float), ('LocalChord', np.float), ('Cl', np.float), ('Cdw', np.float), ('Cd0', np.float), ('Vep_total', np.float), ('V_r_effects', np.float) ]
        structArray = np.zeros((len(lift[:,1]),),dtype=dtype)
        structArray['Yposi'] = lift[:, 0]
        structArray['Area'] = lift[:, 1]
        structArray['LocalChord'] = lift[:, 2]
        structArray['Cl'] = lift[:, 3]
        structArray['Cdw'] = lift[:, 4]
        structArray['Cd0'] = lift[:, 5]
        structArray['Vep_total'] = lift[:, 6]
        structArray['V_r_effects'] = lift[:,7]

        return np.sort(structArray, order='Yposi')
    
    def SumDistributedCoef(self, DistCoef, plane, V):
        ''' Takes as input the distributed coef
        Returns CL and Cl (lift and rolling moment coefficient)
        
        Recompute the induced velocity and sum the friction drag and prop wash.
        
        The function works with organised coefficients in a dictionnary :
            dtype=[('Yposi',np.float),('Area',np.float),('LocalChord',np.float),('Cl',np.float),('Cdw',np.float),('Cd0',np.float)]
            The data typically comes from a VLM, it should be ordered from -b/2 to b/2
        '''
        
        SortedCoef = self.ReOrganiseLift(DistCoef)

        Vep_total = SortedCoef['Vep_total']
        Vi = SortedCoef['V_r_effects']



        tempRoll = np.sum((-SortedCoef['Yposi']*SortedCoef['Cl']*SortedCoef['Area']*Vi**2))/(plane.b*plane.S*V**2)

        tempCL = np.sum(SortedCoef['Cl'] * SortedCoef['Area'] * Vi**2) / (plane.S * V**2)

        tempCdWash = np.sum(SortedCoef['Area'] * SortedCoef['Cdw'] * Vi**2) / (plane.S * V ** 2)

        tempCd0 = np.sum(SortedCoef['Area'] * SortedCoef['Cd0'] * Vi**2) / (plane.S * V ** 2)




        ### New integration for induced drag.
        # Compute DeltaCL at Panel seperation

        wiadim = np.zeros(len(SortedCoef['Yposi']))







        
        """ The validation cases have full wing flaps which create a large lift differential at the wingtip
        It is thought to be un-realistic based on the results.
        For the validation cases the lift derivative is brought to zero at the extreme segments to avoid drag divergence
        For normal use with flap not extending toward wingtip, the lift derivative has to be maintained """



        """
        #Fast computation no smoothing (better without flaps and rather low Tc):
        Diffcl0 =(SortedCoef['LocalChord'][0]*SortedCoef['Cl'][0] - SortedCoef['LocalChord'][0]*0)/(SortedCoef['Yposi'][0]-(-plane.b/2))
        Diffclend = (0-SortedCoef['LocalChord'][-1]*SortedCoef['Cl'][-1])/((plane.b/2)-SortedCoef['Yposi'][-1])
        Diffcl = np.hstack( (Diffcl0, np.diff((SortedCoef['LocalChord']*SortedCoef['Cl']))/np.diff((SortedCoef['Yposi'])), Diffclend) )
        """


        #Fast computation smoothing (for study with flaps / high Tc)

        #Choice 1: brings Cl to its negative symmetry (vortex)
        Cl = np.hstack((-SortedCoef['LocalChord'][0]*SortedCoef['Cl'][0], SortedCoef['LocalChord']*SortedCoef['Cl'], -SortedCoef['LocalChord'][-1]*SortedCoef['Cl'][-1]))
        
        #Choice 2: brings Cl to 0
#        Cl =np.hstack( (0, SortedCoef['LocalChord']*SortedCoef['Cl'], 0) )




        dY = SortedCoef['Yposi'][-1]-SortedCoef['Yposi'][-2]
        Yextended = np.hstack((SortedCoef['Yposi'][0]-dY, SortedCoef['Yposi'], SortedCoef['Yposi'][-1]+dY))
        Diffcl1 = np.hstack((0, np.diff(Cl)/np.diff(Yextended)))        # diff gives out Cl[i+1]-Cl[i],Cl[i+2] - Cl[i+1] ...
        Diffcl2 = np.hstack((np.diff(Cl)/np.diff(Yextended), 0))
        Diffcl3 = (Diffcl1+Diffcl2)/2
        Diffcl = (Diffcl3[:-1] + Diffcl3[1:])/2


        """
        #Choice 3 fast computation but over-smoothing, use only at large Tc > 0.3
        testgrad1=np.gradient(Cl[1:-1],SortedCoef['Yposi'])
        Diffcl=testgrad1
        deltaij=np.ones(len(testgrad1))
        """



        #Adjust position. Yposi is in the center of the slices, DiffclPosi is in the extremes
        DiffclPosi = np.hstack(((-plane.b/2), SortedCoef['Yposi'][1:] - np.diff(SortedCoef['Yposi'])/2, (plane.b/2)))




        #Compute Downwash distribution by integration; 2 things needed: Diffcl, DiffclPosi
        for i in range(len(wiadim)):
            wiadim[i] = np.trapz( Diffcl/(SortedCoef['Yposi'][i]-DiffclPosi), DiffclPosi)

        """
        #That's for solution 3:
        for i in range(len(wiadim)):
         den = SortedCoef['Yposi'][i]-SortedCoef['Yposi']
         deltaij[i]=0
         den[i]=1
         wiadim[i] = np.trapz(testgrad1*deltaij/(den),SortedCoef['Yposi'])
         deltaij[i]=1
         """

        


        wiadim = wiadim * Vi * 1 / (8 * np.pi)
        if self.PlotDrag == True:
            self.wiadim = wiadim  # save for later plotting


        # Compute new induced drag by integrating downwash wiadim

        Cdi = np.trapz(SortedCoef['LocalChord'] * Vi * SortedCoef['Cl'] * wiadim, SortedCoef['Yposi'])/(plane.S*V**2)

        self.Cdi_vec = np.zeros(len(SortedCoef['Yposi']))
        for i in range(len(SortedCoef['Yposi'])):
           self.Cdi_vec[i] = (SortedCoef['LocalChord'][i] * Vi[i] * SortedCoef['Cl'][i] * wiadim[i]) * (DiffclPosi[i+1] - DiffclPosi[i])/(plane.S*V**2)


        # Compute yaw moment due to asymetric induced velocity: sum cdi_local*ylocal

        tempYaw = np.trapz(SortedCoef['LocalChord'] * SortedCoef['Cl'] * wiadim * SortedCoef['Yposi']*Vi**2, SortedCoef['Yposi']) / (plane.S * plane.b*V**3)

        tempYaw_w = sum(SortedCoef['Area'] * SortedCoef['Cdw'] * (SortedCoef['Yposi'] * Vi**2)) / (plane.b * plane.S * V**2)

        
        if plane.DisplayPatterInfo:
            print('TempYaw = {0:0.5f}, TempYaw_w = {1:0.5f}'.format(tempYaw, tempYaw_w))
            plt.figure()
            plt.plot(SortedCoef['Yposi'], SortedCoef['LocalChord']*SortedCoef['Cl']*wiadim/(plane.c))
            plt.xlabel('Span (m)')
            plt.title('Cdi local')
            plt.grid()
            
            plt.figure()
            plt.plot(DiffclPosi,Diffcl)
            plt.title("Diffcl at panel seperation")
            plt.grid(True)


        
        return np.array([tempCL, tempRoll, Cdi, tempCd0, tempYaw+tempYaw_w, tempCdWash])
    
    def CalcCoef(self, dx, Mach, atmo, aoa, dail, dfl, plane, beta, p, V, r):
        '''
        Returns the coefficient as [CL, Cl, Cdi, Cd0, Cn, Cdw]
        '''
#        self.PlotDist(dx,atmo,aoa,dfl,plane,False) Already in main no need to activate it here

        results = self.SumDistributedCoef(self.PatterJames(dx, Mach, atmo, aoa, dail, dfl, plane, beta, p, V, r), plane, V)
        
        return results

    def PlotDist(self, Tc, Mach, atmo, aoa, dail, dfl, plane, IfSave, beta, p, V, r):
        self.PlotDrag = True  # only here for accompanying drag distribution
        data = self.PatterJames(Tc, Mach, atmo, aoa, dail, dfl, plane, beta, p, V, r)
        Dist = self.ReOrganiseLift(data)
        self.Coef = self.SumDistributedCoef(data, plane, V)

        CL_corrected = (Dist['Cl'] * Dist['V_r_effects']**2) / (V**2)

        self.PlotDrag = False
        plt.figure()                                                                                                    #  Create a new figure, or activate an existing figure.
        plt.plot(Dist['Yposi'], CL_corrected, linestyle='--', color='0.25', label='$T_c$ = {0:0.3f}'.format(Tc[0]))     #  Plot y versus x as lines and/or markers.
        ax = plt.gca()                                                                                                  #  Get the current Axes.
        ax.set_xlabel('Y (m)')                                                                                          #  Writes label for an axe
        ax.set_ylabel('Local $C_L$')
        ax.legend()
        plt.grid()
        plt.tight_layout()
        
        fig1 = plt.figure()
        ax1 = fig1.gca()
        ax1.plot(Dist['Yposi'], self.wiadim/(8*np.pi)*180/np.pi, label="$α_i$, $T_c$ = {0:0.3f}".format(Tc[0]), linestyle='-.', color='0.25')
        ax1.set_xlabel('Y (m)')
        ax1.set_ylabel('Downwash angle (°)')
        ax1.legend()
        ax1.grid()
        fig1.tight_layout()


        fig2 = plt.figure()
        ax2 = fig2.gca()
        plt.plot(Dist['Yposi'], self.Cdi_vec  , linestyle='--', color='0.25', label='$T_c$ = {0:0.3f}'.format(Tc[0]))
        ax2.set_xlabel('Y (m)')
        ax2.set_ylabel('Cd induced')
        ax2.legend()
        ax2.grid()
        fig2.tight_layout()


        fig3 = plt.figure()
        ax3 = fig3.gca()
        plt.plot(Dist['Yposi'], Dist['Cdw'] , linestyle='--', color='0.25', label='$T_c$ = {0:0.3f}'.format(Tc[0]))
        ax3.set_xlabel('Y (m)')
        ax3.set_ylabel('Cd wash')
        ax3.legend()
        ax3.grid()
        fig3.tight_layout()


        fig4 = plt.figure()
        ax4 = fig4.gca()
        plt.plot(Dist['Yposi'],  Dist['Cd0'] , linestyle='--', color='0.25', label='$T_c$ = {0:0.3f}'.format(Tc[0]))
        ax4.set_xlabel('Y (m)')
        ax4.set_ylabel('Cd0_extra ')
        ax4.legend()
        ax4.grid()
        fig4.tight_layout()




        fig5 = plt.figure()
        ax5 = fig5.gca()
        plt.plot(Dist['Yposi'], self.Cdi_vec + Dist['Cdw'] + Dist['Cd0'] , linestyle='--', color='0.25', label='$T_c$ = {0:0.3f}'.format(Tc[0]))
        ax5.set_xlabel('Y (m)')
        ax5.set_ylabel('Cd_wash + Cd induced + Cd0 ')
        ax5.legend()
        ax5.grid()
        fig5.tight_layout()



        fig6 = plt.figure()
        ax6 = fig6.gca()
        ax6.plot(Dist['Yposi'], Dist['Cl'], label="$α_i$, $T_c$ = {0:0.3f}".format(Tc[0]), linestyle='-.', color='0.25')
        ax6.set_xlabel('Y (m)')
        ax6.set_ylabel('CL not corrected')
        ax6.legend()
        ax6.grid()
        fig6.tight_layout()


        plt.show(block=True)   # added to plot correctly

        if IfSave:
            plt.savefig('./CurrentLiftRepartition.pdf')
        
        return

    def PatterJames(self, Tc, Mach, atmospher, aoa, dail, dfl, plane,beta,p,V,r):
        '''
        This function computes the prop-wing interaction lift increase and friction drag due to blowing
        '''
        #get atmosphere info
        rho=atmospher[1]
        self.aoa = aoa #store locally for drag
        
        #In this version of non-dim Augmented Patterson, Tc must be nn-dim: tc = T/(2*rho*Sp*V**2)
        """ The function could be modified here to take directly Tc as input
        Or a different thrust model can be entered in 'aircraft' class"""
        
        #solve equation 8-2 from "McCormick, Aerodynamics of V/STOL" for Vp/V
        # Vp: propeller axial induced velocity
        myw=np.zeros(len(Tc))
        self.mu=np.zeros(len(Tc))


        #get wing alpha0
        alpha0w = self.Interpol(self.AoAZero, Mach) # test with section zero lift angle
        alpha0w=alpha0w[:,-1] # keep only alpha0

        #Get the local slope
        NormCl = self.Interpol(self.CLslope, Mach)






        """ 
        In equation 3.20, 3.21 Vinf should be the equation seen by each engine, instead of considering uniform,
        here a different speed is seen by each engine due to r and sideslips angle, so thrust coefficient cannot be made
        to appear as in 3.21
        """





        av_alpha_0 = np.mean(alpha0w)

        V_vect = V * (np.cos((-np.sign(plane.PosiEng)) * beta + plane.wingsweep)) - r * plane.PosiEng
        Velocity = V * (np.cos((-np.sign(NormCl[:, 0])) * beta + plane.wingsweep)) - r * NormCl[:, 0]
        T = Tc * (2 * rho * V ** 2 * plane.Sp)



        for i in range(len(Tc)):
            if Tc[i] == 0:
                #No Thrust, no need to solve the equation
                myw[i] = 0
            else:
                coef = [1, 2*np.cos(aoa-av_alpha_0+plane.alpha_i+plane.ip), 1, 0, - (T[i] / (2 * rho * plane.Sp * V_vect[i]**2)) ** 2]
                roots = np.roots(coef)
                #get the real positive root
                for j in range(len(roots)):
                    if np.real(roots[j]) > 0:
                        myw[i] = np.real(roots[j])
            # test the negative thrust effects by simply setting negative roots
            if Tc[i] < 0:
                self.mu[i] = -2*myw[i]
            else:
                self.mu[i] = 2*myw[i]



             # Be careful, what you get is Vp /V_vect, not Vep/V_vect
             # According to momenthum theory, inmediatly after the propeller the speed is Var_V. Far down-wash, it is 2Var_V.
             # Here we take Vp = 2Var_V.













        # Compute AoA modification due to flaps
        if plane.isflap:
            alpha_fl = aoa - self.alpha0_fl * dfl
        else:
            alpha_fl = 0
        
        
        #compute AoA modification due to ailerons
        if plane.isail:
            #Aileron differential
            if dail>0:
                dail_l = -dail
                dail_r = dail*plane.AilDiff
            else:
                dail_l = -dail*plane.AilDiff
                dail_r = dail
            alpha_ail_l = aoa - self.alpha0_ail * dail_l
            alpha_ail_r = aoa - self.alpha0_ail * dail_r
        else:
            alpha_ail_l = aoa
            alpha_ail_r = aoa


        alpha_t       =  aoa        - alpha0w + plane.alpha_i + beta*plane.dihedral*np.sign(NormCl[:, 0]) + p * NormCl[:, 0]/Velocity[:]
        alpha_fl_t    = alpha_fl    - alpha0w + plane.alpha_i + beta*plane.dihedral*np.sign(NormCl[:, 0]) + p * NormCl[:, 0]/Velocity[:]
        alpha_ail_t_l = alpha_ail_l - alpha0w + plane.alpha_i + beta*plane.dihedral*np.sign(NormCl[:, 0]) + p * NormCl[:, 0]/Velocity[:]
        alpha_ail_t_r = alpha_ail_r - alpha0w + plane.alpha_i + beta*plane.dihedral*np.sign(NormCl[:, 0]) + p * NormCl[:, 0]/Velocity[:]

        #corresponding alpha max, assume aileron stall angle is alpha_t_max:
#        alpha_t_max = plane.alpha_max + alpha0w
#        alpha_fl_t_max = plane.alpha_max_fl + alpha0w - self.alpha0_fl * dfl
        
        alpha_t_max = plane.alpha_max * np.ones_like(alpha0w)
        alpha_fl_t_max = plane.alpha_max_fl * np.ones_like(alpha0w) - self.alpha0_fl * dfl
        
        
        #Determine if section is behind propeller, has flap or ailerons
        SectInProp = np.zeros(len(NormCl[:,1]))
        SectHasFlap = [False]*len(NormCl[:,1])
        SectHasAilLeft = [False]*len(NormCl[:,1])
        SectHasAilRight = [False]*len(NormCl[:,1])




                
        if plane.isflap:
             #Flap 1, negative y
            Fl1Tip = -plane.FlPosi*plane.b-plane.FlRatio*plane.b/2
            Fl1Root = -plane.FlPosi*plane.b
            
            # Flap 2, positive y
            Fl2Tip = plane.FlPosi*plane.b+plane.FlRatio*plane.b/2
            Fl2Root = plane.FlPosi*plane.b
        else:
            Fl1Tip = 0
            Fl1Root = 0
            Fl2Tip = 0
            Fl2Root =0
        
        if plane.isail:
            # aileron 1 negative y
            Ail1Tip = -plane.AilPosi*plane.b - plane.AilRatio*plane.b
            Ail1Root = -plane.AilPosi*plane.b
            
            #Aileron 2, positive y
            Ail2Tip = plane.AilPosi*plane.b + plane.AilRatio*plane.b
            Ail2Root = plane.AilPosi*plane.b
        else:
            Ail1Tip = 0
            Ail1Root = 0
            Ail2Tip = 0
            Ail2Root = 0
        
#        print([Fl1Tip,Fl1Root,Fl2Tip,Fl2Root])
                
        for i in range(len(SectInProp)):
            # Check is sect has prop, flap, ailerons left or right
            for a in range(len(Tc)):
                if NormCl[i,0] <= plane.PosiEng[a]+plane.Dp/2 and NormCl[i,0] >= plane.PosiEng[a]-plane.Dp/2:
                    SectInProp[i] = int(a+1) # label engines from 1 to N_eng
            
            #Flap 1 negative y
            if NormCl[i,0] <= Fl1Root and NormCl[i,0] >= Fl1Tip:
                SectHasFlap[i] = True

            # Flap 2, positive y
            elif NormCl[i,0] <= Fl2Tip and NormCl[i,0] >= Fl2Root:
                SectHasFlap[i] = True
                
            #Left aileron strict inequalities
            elif NormCl[i,0] < Ail1Root and NormCl[i,0] >= Ail1Tip:
                SectHasAilLeft[i] = True
            
            #Left right aileron strict inequalities
            elif NormCl[i,0] <= Ail2Tip and NormCl[i,0] > Ail2Root:
                SectHasAilRight[i] = True
                
        #For drag, no matter what, sections at -b/2 and b/2 are outside propeller influence
        SectInProp[int(len(SectInProp)/2)-1] = 0
        SectInProp[-1] = 0
                
        # Can compute the surrogate for beta
        BetaVec = self.BetaSurro(plane, self.mu+1, rho, SectInProp)
        self.Beta = BetaVec
            
        #lift and drag multiplier from patterson calculation
        LmFl = np.zeros(len(NormCl[:, 1]))
        alpha_ep_drag = np.zeros(len(NormCl[:, 1]))
        alpha_ep = np.zeros(len(NormCl[:, 1]))
        LocalCl = np.copy(NormCl)
        self.PWashDrag = np.zeros(int(len(LocalCl))).reshape(int(len(LocalCl)), 1)
        self.Cd0_vec = np.zeros(int(len(LocalCl))).reshape(int(len(LocalCl)), 1)
        self.LocalVelocity = np.zeros(int(len(LocalCl))).reshape(int(len(LocalCl)), 1)

        
        
        for i in range(len(SectInProp)):
            if SectInProp[i] != 0:
                #Take engine's number from 0
                j=int(SectInProp[i]-1)
                #Compute lift multiplier, aoa and local drag with flap
                if SectHasFlap[i]:
                    LmFl[i] = ( 1 - BetaVec[i]*self.mu[j]*np.sin(plane.ip+self.alpha0_fl * dfl)/(np.sin(alpha_fl_t[i]))) * (1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_t[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)**0.5-1
                    alpha_ep_drag[i] = (alpha_fl_t[i] - self.mu[j]*(plane.ip+self.alpha0_fl * dfl)) /(1+self.mu[j])-alpha_fl_t[i] + p * NormCl[i, 0]/Velocity[i]    #+ self.AoAZero[i,-1]
                    self.Cd0_vec[i] = plane.Cd0_turbulent*((1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_fl_t[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)-1)
                    #self.Cd0_vec[i] = plane.Cd0_turbulent*(1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_fl_t[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)
                    self.LocalVelocity[i] = (1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_fl_t[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2) * V

                    alpha_ep[i] = np.arctan((alpha_fl_t[i] - self.mu[j]*(plane.ip+self.alpha0_fl * dfl)) /(1+self.mu[j]))
                    if alpha_ep[i] < alpha_fl_t_max[i]:
                        LocalCl[i, -1] = LocalCl[i, -1] * alpha_fl_t[i]
                    else:
                        LocalCl[i, -1] = LocalCl[i, -1] * np.sin(alpha_fl_t_max[i])*np.cos(alpha_fl_t[i])/np.cos(alpha_fl_t_max[i])
                        if alpha_fl_t[i] < self.alphaDrag[-1]:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(alpha_fl_t[i])
                        else:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(self.alphaDrag[-1])+np.sin(alpha_fl_t_max[i])*np.sin(alpha_fl_t[i])/np.cos(alpha_fl_t_max[i])
                
                elif SectHasAilLeft[i]:
                    #Compute lift multiplier, aoa and local drag with aileron
                    LmFl[i] = ( 1 - BetaVec[i]*self.mu[j]*np.sin(plane.ip+self.alpha0_ail * dail_l)/(np.sin(alpha_ail_t_l[i]))) * (1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_ail_t_l[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)**0.5-1
                    alpha_ep_drag[i] = (alpha_ail_t_l[i] - self.mu[j]*(plane.ip+self.alpha0_ail * dail_l)) /(1+self.mu[j])-alpha_ail_t_l[i] + p * NormCl[i, 0]/Velocity[i]#+ self.AoAZero[i,-1]
                    self.Cd0_vec[i] = plane.Cd0_turbulent*((1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_ail_t_l[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)-1)
                    #self.Cd0_vec[i] = plane.Cd0_turbulent*(1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_ail_t_l[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)
                    self.LocalVelocity[i] = (1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_ail_t_l[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2) *V
                    
                    alpha_ep[i] = np.arctan((alpha_ail_t_l[i] - self.mu[j]*(plane.ip+self.alpha0_ail * dail_l)) /(1+self.mu[j]))
                    if alpha_ep[i] < alpha_t_max[i]:
                        LocalCl[i, -1] = LocalCl[i, -1] * alpha_ail_t_l[i]
                    else:
                        LocalCl[i, -1] = LocalCl[i, -1] * np.sin(alpha_t_max[i])*np.cos(alpha_ail_t_l[i])/np.cos(alpha_t_max[i])
                        if alpha_ail_t_l[i] < self.alphaDrag[-1]:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(alpha_ail_t_l[i])
                        else:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(self.alphaDrag[-1])+np.sin(alpha_t_max[i])*np.sin(alpha_ail_t_l[i])/np.cos(alpha_t_max[i])
                
                elif SectHasAilRight[i]:
                    #Compute lift multiplier, aoa and local drag with aileron
                    LmFl[i] = ( 1 - BetaVec[i]*self.mu[j]*np.sin(plane.ip+self.alpha0_ail * dail_r)/(np.sin(alpha_ail_t_r[i]))) * (1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_ail_t_r[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)**0.5-1
                    alpha_ep_drag[i] = (alpha_ail_t_r[i] - self.mu[j]*(plane.ip+self.alpha0_ail * dail_r)) /(1+self.mu[j])-alpha_ail_t_r[i] + p * NormCl[i, 0]/Velocity[i]#+ self.AoAZero[i,-1]
                    self.Cd0_vec[i] = plane.Cd0_turbulent*((1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_ail_t_r[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)-1)
                    #self.Cd0_vec[i] = plane.Cd0_turbulent*(1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_ail_t_r[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)
                    self.LocalVelocity[i] = (1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_ail_t_r[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2) *V
                    
                    alpha_ep[i] = np.arctan((alpha_ail_t_r[i] - self.mu[j]*(plane.ip+self.alpha0_ail * dail_r)) /(1+self.mu[j]))
                    if alpha_ep[i] < alpha_t_max[i]:
                        LocalCl[i, -1] = LocalCl[i, -1] * alpha_ail_t_r[i]
                    else:
                        LocalCl[i, -1] = LocalCl[i, -1] * np.sin(alpha_t_max[i])*np.cos(alpha_ail_t_r[i])/np.cos(alpha_t_max[i])
                        if alpha_ail_t_r[i] < self.alphaDrag[-1]:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(alpha_ail_t_r[i])
                        else:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(self.alphaDrag[-1])+np.sin(alpha_t_max[i])*np.sin(alpha_ail_t_r[i])/np.cos(alpha_t_max[i])
                
                else:
                    #Compute lift multiplier, aoa and local drag on clean wing
                    LmFl[i] = (1 - BetaVec[i]*self.mu[j]*np.sin(plane.ip)/(np.sin(alpha_t[i]))) * (1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_t[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)**0.5-1
                    alpha_ep_drag[i] = (alpha_t[i] - self.mu[j] * plane.ip) /(1+self.mu[j])-alpha_t[i] + p * NormCl[i, 0]/Velocity[i]#+ self.AoAZero[i,-1]
                    self.Cd0_vec[i] = plane.Cd0_turbulent*((1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_t[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)-1)
                    #self.Cd0_vec[i] = plane.Cd0_turbulent*(1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_t[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2)
                    self.LocalVelocity[i] = (1 + 2*self.mu[j]*BetaVec[i]*np.cos(alpha_t[i] + plane.ip) + (BetaVec[i]*self.mu[j])**2) *V

                    
                    alpha_ep[i] = np.arctan((alpha_t[i] - self.mu[j] * plane.ip) /(1+self.mu[j]))
                    if alpha_ep[i] < alpha_t_max[i]:
                        LocalCl[i, -1] = LocalCl[i, -1] * alpha_t[i]
                    else:
                        LocalCl[i, -1] = LocalCl[i, -1] * np.sin(alpha_t_max[i])*np.cos(alpha_t[i])/np.cos(alpha_t_max[i])
                        if alpha_t[i] < self.alphaDrag[-1]:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(alpha_t[i])
                        else:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(self.alphaDrag[-1])+np.sin(alpha_t_max[i])*np.sin(alpha_t[i])/np.cos(alpha_t_max[i])
            else:
                LmFl[i] = 0
                alpha_ep_drag[i] = 0 + p * NormCl[i, 0]/Velocity[i]
                self.Cd0_vec[i] = 0
                self.LocalVelocity[i] = V * (np.cos((-np.sign(NormCl[i, 0])) * beta + plane.wingsweep)) - r * NormCl[i, 0]
                if SectHasFlap[i]:
                    alpha_ep[i] = alpha_fl_t[i]
                    # Check stall
                    if alpha_fl_t[i] < alpha_fl_t_max[i]:
                        LocalCl[i, -1] = LocalCl[i, -1] * alpha_fl_t[i]
                    else:
                        LocalCl[i,-1] = LocalCl[i, -1] * np.sin(alpha_fl_t_max[i])*np.cos(alpha_fl_t[i])/np.cos(alpha_fl_t_max[i])
                        if alpha_fl_t[i] < self.alphaDrag[-1]:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(alpha_fl_t[i])
                        else:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(self.alphaDrag[-1])+np.sin(alpha_fl_t_max[i])*np.sin(alpha_fl_t[i])/np.cos(alpha_fl_t_max[i])
                
                elif SectHasAilLeft[i]:
                    alpha_ep[i] = alpha_ail_t_l[i]
                    # Check stall
                    if alpha_ail_t_l[i] < alpha_t_max[i]:
                        LocalCl[i, -1] = LocalCl[i, -1] * alpha_ail_t_l[i]
                    else:
                        LocalCl[i, -1] = LocalCl[i, -1] * np.sin(alpha_t_max[i])*np.cos(alpha_ail_t_l[i])/np.cos(alpha_t_max[i])
                        if alpha_ail_t_l[i] < self.alphaDrag[-1]:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(alpha_ail_t_l[i])
                        else:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(self.alphaDrag[-1])+np.sin(alpha_t_max[i])*np.sin(alpha_ail_t_l[i])/np.cos(alpha_t_max[i])
                        
                elif SectHasAilRight[i]:
                    alpha_ep[i] = alpha_ail_t_r[i]
                    # Check stall
                    if alpha_ail_t_r[i] < alpha_t_max[i]:
                        LocalCl[i, -1] = LocalCl[i, -1] * alpha_ail_t_r[i]
                    else:
                        LocalCl[i, -1] = LocalCl[i, -1] * np.sin(alpha_t_max[i])*np.cos(alpha_ail_t_r[i])/np.cos(alpha_t_max[i])
                        if alpha_ail_t_r[i] < self.alphaDrag[-1]:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(alpha_ail_t_r[i])
                        else:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(self.alphaDrag[-1])+np.sin(alpha_t_max[i])*np.sin(alpha_ail_t_r[i])/np.cos(alpha_t_max[i])
                
                else:
                    alpha_ep[i] = alpha_t[i]
                    if alpha_t[i] < alpha_t_max[i]:
                        LocalCl[i, -1] = LocalCl[i, -1] * alpha_t[i]
                    else:
                        LocalCl[i, -1] = LocalCl[i, -1] * np.sin(alpha_t_max[i])*np.cos(alpha_t[i])/np.cos(alpha_t_max[i])
                        if alpha_t[i] < self.alphaDrag[-1]:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(alpha_t[i])
                        else:
                            self.Cd0_vec[i] = self.Cd0_vec[i]+self.StallDrag(self.alphaDrag[-1])+np.sin(alpha_t_max[i])*np.sin(alpha_t[i])/np.cos(alpha_t_max[i])
        
        BlownCl = np.copy(LocalCl)
        BlownCl[:, -1] = LocalCl[:, -1]*(LmFl+1) * np.cos(-alpha_ep_drag) *self.DeltaCL_a_0
        self.PWashDrag[:, 0] = BlownCl[:, -1] * np.sin(-alpha_ep_drag)
        self.LmFl = LmFl
        self.alpha_t_max = alpha_t_max
        self.alpha_fl_t_max = alpha_fl_t_max
        self.alpha_ail_t_l = alpha_ail_t_l
        self.alpha_ep = alpha_ep

        Vel = np.zeros(int(len(LocalCl))).reshape(int(len(LocalCl)), 1)
        Vel[:, 0] = Velocity

        return np.hstack((np.hstack((BlownCl, self.PWashDrag)), self.Cd0_vec, self.LocalVelocity, Vel))



    """
    alpha_ep_drag = alpha_ep - alpha_t. Its used for calculating the extra drag given by the lift when lift is
    deflected an angle alpha_ep . Negative sign in self.PWashDrag[:, 0] = BlownCl[:, -1] * np.sin(-alpha_ep_drag) 
    as  formule is  Cd_i,w = CL * sen (alpha-alpha_ep)  Page 75 Eric Nguyen's thesis
    
    self.DeltaCL_a_0 CL_alpha correction factor, equal to 1. Defined in Main. For tunning
    
    
    self.StallDrag = interp1d(self.alphaDrag,self.StallDrag)   interpolates between the columns alpha and 
    CD of the file naca3318Pol.txt that contains aerodynamic info about the airfoil.
    Is used if we are in stall. Finally there are two cases, alpha slower than the higher angle of the file
    and alpha above. The higher angle of the file is 26.5 ...
    How to model drag in stall I believe it comes in Jamesson as well... check
    """









    def Augmented_velocity_wing(self, Tc, Mach, atmospher, aoa, dail, dfl, plane, beta, p, V, r):

        rho = atmospher[1]


        V_vect = V * (np.cos((-np.sign(plane.PosiEng)) * beta + plane.wingsweep)) - r * plane.PosiEng
        T = Tc * (2 * rho * V ** 2 * plane.Sp)

        myw = np.zeros(len(Tc))
        muu = np.zeros(len(Tc))

        for i in range(len(Tc)):
            if Tc[i] == 0:
                 #No Thrust, no need to solve the equation
                 myw[i] = 0
            else:
                 coef = [1, 2*np.cos(aoa-plane.alpha_0+plane.ip), 1, 0, - (T[i] / (2 * rho * plane.Sp * V_vect[i]**2)) ** 2]
                 roots = np.roots(coef)
                 #get the real positive root
                 for j in range(len(roots)):
                     if np.real(roots[j])>0:
                         myw[i]=np.real(roots[j])
            # test the negative thrust effects by simply setting negative roots
            if Tc[i] < 0:
                 muu[i] = 2*myw[i]
            else:
                 muu[i] = 2*myw[i]



        #get wing alpha0
        alpha0w = self.Interpol(self.AoAZero, Mach)  # test with section zero lift angle
        alpha0w = alpha0w[:, -1]  # keep only alpha0

        #Get the local slope
        NormCl = self.Interpol(self.CLslope, Mach)


        Velocity = V * (np.cos((-np.sign(NormCl[:, 0])) * beta + plane.wingsweep)) - r * NormCl[:, 0]
        alpha_t = aoa - alpha0w + plane.alpha_i + beta*plane.dihedral*np.sign(NormCl[:, 0]) + p * NormCl[:, 0]/Velocity[:]
        alpha_t_mean = np.mean(alpha_t)



        V_ef_TO_V_inf = (1 + 2*muu[0]*np.cos(alpha_t_mean + plane.ip) + (muu[0])**2 )**0.5
        Pressure_ratio = V_ef_TO_V_inf**2


        return Pressure_ratio