# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 16:29:42 2018

@author: e.nguyen-van
"""

from math import pi
from scipy.interpolate import interp1d
import numpy as np
import sys

class Propu:
    # common data go here
    rho=1.225
    AvGearEff=0.80 # average effiency of engine+gearbox
    
    #definition
    def __init__(self, PropDia, PropPas, Ct = 0.05, Cp = 0.07):
        #default value 11"10" at J=0.6 wind tunnel on
        # PropDia and PropPas in inch and float
        self.D=PropDia*0.0254
        self.Pas=PropPas
        
        # Try to find selig prop data file from here : https://m-selig.ae.illinois.edu/props/propDB.html
        DataFileAvailable=True
        DiaTest=0
        try :
            file=open("SeligProp"+str(PropDia)+str(PropPas)+".txt",'r')
            
        except OSError:
            print("No data file for this propeller, try different Diameter")
            # if a propeller with different diameter but same pitch is found, 
            #its data are used and imposed diameter used for non-dimensional analysis
            # Not the best but ok, because the forces are non-dimensionalised by the diameter. The opposite cannot be done
            DiatoTest=5
            for i in range(DiatoTest):
                try:
                    DiaTest=PropDia+i
                    file=open("SeligProp"+str(DiaTest)+str(PropPas)+".txt",'r')
                    i=DiatoTest
                    DataFileAvailable=True
                    break
                except OSError:
                    DataFileAvailable=False
            if DataFileAvailable==False:
                for i in range(-DiatoTest,0):
                    try :
                        DiaTest=PropDia+i
                        file=open("SeligProp"+str(DiaTest)+str(PropPas)+".txt",'r')
                        print(DiaTest)
                        DataFileAvailable=True
                        break
                    except OSError:
                        DataFileAvailable=False

        
        if DataFileAvailable==False:
            #conduct analysis with fixed Cp, Ct, know what you do, know where you are on the J interval
            print("No data file for this propeller, exit program in SeligPropellerRead")
            sys.exit()
            
        else:
            if DiaTest!=0:
                print("Found data file for prop diameter : {0:0.0f}\"".format(DiaTest))
            #read stuff, carefull format
            self.Ct=np.array([])
            self.Cp=np.array([])
            self.J=np.array([])
            Heading=file.readlines(1)
            Heading_split = Heading[0].split()
            Jposi=Heading_split.index('J')
            CTposi=Heading_split.index("CT")
            CPposi=Heading_split.index("CP")

            for line in file:
                words=line.split()
                if len(words)!=0:
                    self.J=np.append(self.J,float(words[Jposi]))
                    self.Ct=np.append(self.Ct,float(words[CTposi]))
                    self.Cp=np.append(self.Cp,float(words[CPposi]))
            file.close()
            
        self.cpJ3 = self.Cp/self.J**3 # for determination of J from deltax

    #Functions
    def rpm2rps(self, a):
        return a/60
    
    def rps2rpm(self, a):
        return a*60
    
    def rpm2rads(self, a):
        return a/60*2*pi
    
    def rps2rads(self, a):
        return a*2*pi
    
    def rads2rps(self, a):
        return a/(2*pi)
    
    def CalcJ(self, n, V):
        return V/(self.D*self.rpm2rps(n))
    
    def DetermineJ(self, dx, V, pmax):
        # !!! pmax per engine !!!
        f=interp1d(self.cpJ3, self.J)
        return f(dx*pmax/(self.rho*V**3*self.D**2))
    
    def getCt(self, n=0, V=0,J=0):
        # interpol Ct, works with one vector and a float or two float
        #Works either by entering n and V or J not all three
        #Examples : self.getCt(8000,15) or self.getCt(J=0.8)
        # n in rpm
        #V in m/s
        if np.size(J) >1:
            if np.size(V)>1 or np.size(n)>1:
                print("ERROR: in SeligPropeller.getCt, J and V or n array given, only one or the other allowed")
                sys.exit()
            else:
                f=interp1d(self.J,self.Ct)
            return f(J)
        
        elif np.size(n) > 0 :
            if np.size(V)!=np.size(n):
                print("ERROR: in SeligPropeller.getCt, size of V != size of n.")
                sys.exit()
            else:
                J = V/(self.D*self.rpm2rps(n))
                f=interp1d(self.J,self.Ct)
                return f(J)
        else:
            print("Wrong input sequence for getCt, enter either J or n and V different than 0")
            return 0
        
    def getCp(self, n=0, V=0, J=0):
        # interpol Cp, works with one vector and a float or two float
        #Works either by entering n and V or J not all three
        #Examples : self.getCt(8000,15) or self.getCt(J=0.8)
        #n in rpm
        #V in m/s
        if np.size(J) >1:
            if np.size(V)>1 or np.size(n)>1:
                print("ERROR: in SeligPropeller.getCt, J and V or n array given, only one or the other allowed")
                sys.exit()
            else:
                f=interp1d(self.J,self.Cp)
            return f(J)
        elif np.size(n) > 0 :
            if np.size(V)!=np.size(n):
                print("ERROR: in SeligPropeller.getCt, size of V != size of n.")
                sys.exit()
            else:
                J = V/(self.D*self.rpm2rps(n))
                f=interp1d(self.J,self.Cp)
                return f(J)
        else:
            print("Wrong input sequence for getCp, enter either J or n and V different than 0")
            return 0
    
    def Thrust(self,n,V):
        #compute thurst at given regime and velocity
        #n in rpm
        return self.getCt(n,V)*self.rho*self.D**4*self.rpm2rps(n)**2
    
    def PropPower(self,n,V):
        #compute thurst at given regime and velocity
        #n in rpm
        return self.getCp(n,V)*self.rho*self.D**5*self.rpm2rps(n)**3
    
    def getn(self, J, V):
        return V/(self.D*J)
            