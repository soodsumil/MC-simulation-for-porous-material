# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 18:46:51 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from dataclasses import dataclass, replace
from math import sqrt, asin,acos, pi, sin, cos, exp
import random



class Light_properties:
    def __init__(self, nta1 = 1.003, nta2 = 1.003,):
        self.nta1 = nta1
        self.nta2 = nta2
        
    def snells_law(self, position, direction, center):
        normal1 = np.asarray(self.position) - np.asarray(center)
        norm = LA.norm(normal1)
        normal = normal1/norm
        phi1 = acos(np.dot(-normal,np.asarray(direction)))
        return (asin(self.nta1*sin(phi1)/self.nta2))
    
    def fresnels_law(self, position, direction, center):
        normal1 = np.asarray(self.position) - np.asarray(center)
        norm = LA.norm(normal1)
        normal = normal1/norm
        phi1 = acos(np.dot(-normal,np.asarray(direction)))
        phi2 = asin(self.nta1*sin(phi1)/self.nta2)
        rho_parallel = (self.nta2*cos(phi1)-self.nta1*cos(phi2))/(self.nta2*cos(phi1)+self.nta2*cos(phi2))
        rho_perpendicular = (self.nta1*cos(phi1)-self.nta2*cos(phi2))/(self.nta1*cos(phi1)+self.nta2*cos(phi2))

        trans_parallel = (2*sin(phi2)*cos(phi1))/(sin(phi1+phi2)*cos(phi1-phi2))
        trans_perpendicular = (2*sin(phi2)*cos(phi1))/(sin(phi1+phi2))
        
        rho_avg = (rho_parallel**2 + rho_perpendicular**2)/2
        trans_avg = (trans_parallel**2 + trans_perpendicular**2)/2
        
        return [rho_avg,trans_avg]
    
    
    
    
    
    
    
    
    
    
    