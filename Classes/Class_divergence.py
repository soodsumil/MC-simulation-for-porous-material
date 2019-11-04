# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 16:04:37 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from dataclasses import dataclass, replace
from math import sqrt, asin,acos, pi, sin, cos, exp
import random

class Divergence:
    def __init__(self, position = [0, 0, 0]):
        self.position = position
        
    def cone_divergence(self):
        R_gamma1 = np.random.random(1)[0] #randomly getting the center of the sphere being hit
        R_theta1 = np.random.random(1)[0]
        gamma1 = acos(R_gamma1)
        theta1 = 2*pi*R_theta1
        new_direction= np.array([sin(gamma1)*cos(theta1), sin(gamma1)*sin(theta1), cos(gamma1)])
        return new_direction
    
    def circular_divergence(self):
        R_gamma1 = np.random.random(1)[0] #randomly getting the center of the sphere being hit
        R_theta1 = np.random.random(1)[0]
        gamma1 = 2*pi*R_gamma1
        theta1 = 2*pi*R_theta1
        new_direction= np.array([sin(gamma1)*cos(theta1), sin(gamma1)*sin(theta1),cos(gamma1)])
        return new_direction
    
    def rotation(self, direction):
        k1 = np.array([0,0,1])
        k2_ = np.asarray(direction)
        norm1 =LA.norm(k2_)
        k2 = k2_/norm1
        norm2 = np.cross(k2,k1)
        i2 = (np.cross(k2,k1))/LA.norm(norm2)
        j2 = np.cross(k2,i2)
        rot = np.array([[i2],[j2],[k2]])
        rot = rot.transpose() 
        return rot
    
# =============================================================================
# x = Divergence()
# dir = x.circular_divergence()*0.476
# rot = x.rotation(dir)
# print(rot)
# dirs = x.cone_divergence()
# q = np.matmul(rot,dirs)
# q = np.reshape(q,(1,3))
# print(q)
# =============================================================================
    
    
    
    
    
    