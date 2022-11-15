#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:51:04 2022

@author: saumya
"""

from math import pi
from scipy.interpolate import interp1d
import numpy as np


class Prop:
    rho = 1.225    
    def __init__(self, D,Pitch, n_max, p_max):
       
        self.D=D
        self.Pitch=Pitch
        self.n_max = n_max
        self.p_max = p_max
        file=open("SeligProp86.txt",'r')
        
        # Read Files
        self.Ct=np.array([])
        self.Cp=np.array([])
        self.J=np.array([])
        self.eta=np.array([])
        Heading=file.readlines(1)
        Heading_split = Heading[0].split()
        Jposi=Heading_split.index('J')
        CTposi=Heading_split.index("CT")
        CPposi=Heading_split.index("CP")
        etaposi = Heading_split.index("eta")

        for line in file:
            words=line.split()
            if len(words)!=0:
                self.J=np.append(self.J,float(words[Jposi]))
                self.Ct=np.append(self.Ct,float(words[CTposi]))
                self.Cp=np.append(self.Cp,float(words[CPposi]))
                self.eta= np.append(self.eta,float(words[etaposi]))
          
        file.close()
        
        # The limit is tested with respect to a fraction of power, now dx represents the fraction of Power that can be used with resoect to the maximum power
        self.J[0]=0.0009 # Otherwise will show NaN data
        self.CpJ3 = self.Cp / self.J ** 3
        # Interpolation of coefiicients of Thrust and Power with the Advance ratio
        self.Ct_f_J=interp1d(self.J,self.Ct, kind = 'quadratic')
        self.Cp_f_J=interp1d(self.J,self.Cp, kind = 'quadratic')
        self.J_f_CpJ3 = interp1d(self.CpJ3, self.J, kind = 'quadratic')
   
    
    def DetermineJ(self, dx, V):

        self.result_J = np.zeros((len(V)))
        a=np.zeros((len(V)))
        for i in range(len(V)):
            a[i]= (dx[i]*self.p_max / (self.rho * V[i] ** 3 * self.D ** 2))
            if a[i]>self.CpJ3[1]:
                self.result_J[i] = self.J[1]
            elif a[i]<self.CpJ3[-2]:
                self.result_J[i] = self.J[-2]
            else:
                self.result_J[i] = self.J_f_CpJ3(a[i])
        self.n = V/(self.D*self.result_J);
        # Maybe add a check for n_max

    def getCt(self):
        
        results=np.array([])
        for i in range(len(self.result_J)):
            results = np.append(results, self.Ct_f_J(self.result_J[i]))
        return results
        
        
    def getCp(self):
        
        results=np.array([])
        for i in range(len(self.result_J)):
            results = np.append(results, self.Cp_f_J(self.result_J[i]))
        return results
    
    def Thrust(self,  dx, V):
        

        self.DetermineJ(dx, V)
        return self.getCt()*self.rho*self.D**4*self.n**2
    
    def PropPower(self,  dx, V):
        
        self.DetermineJ(dx, V)
        return self.getCp()*self.rho*self.D**5*self.n**3

    
            
    