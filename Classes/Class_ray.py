# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 12:28:14 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from dataclasses import replace
from math import sqrt, asin,acos, pi, sin, cos, exp
from Class_divergence import Divergence

class Ray:
    def __init__(self, position, direction,  rotation = None):
        #position and direction as lists and rotation as numpy array
        self.position = position
        self.direction = direction
        self.rotation = rotation
        
    def propagate(self, distance: float):
        new_position = np.asarray(self.position) + np.asarray(self.direction) * distance
        new_position = new_position.tolist()
        new_ray = replace(self, position=new_position)
        return new_ray
    
    def representation(self, distance):
        new_direction = np.mutual(self.rotation, np.asarray(self.direction))
        new_position = np.asarray(self.position) + new_direction * distance
        new_ray = replace(self, position= new_position, direction=new_direction)
        return new_ray
    
    def specular_reflection(self,center: float): # center is inputed as list
        normal1 = np.asarray(self.position) - np.asarray(center)
        norm = LA.norm(normal1)
        normal = normal1/norm
        dot_prod = np.dot(normal,self.direction)
        reflected_dir_cos1 = np.array([(self.direction[0]-2*normal[0]*dot_prod), (self.direction[1] - 2*normal[1]*dot_prod), (self.direction[2] - 2*normal[2]*dot_prod)])
        norm1 = LA.norm(reflected_dir_cos1)
        reflected_dir_cos = reflected_dir_cos1/norm1
        return reflected_dir_cos
    
    def diffuse_reflection(self, center):
        normal1 = np.asarray(self.position) - np.asarray(center)
        norm = LA.norm(normal1)
        normal = normal1/norm
        K = Divergence()
        rot = K.rotation(normal)
        new_dir1 = K.cone_divergence()
        new_dir = np.matmul(rot,new_dir1)
        return new_dir
        
    
    def refraction(self, diameter, center, nta1, nta2):
        n = nta1/nta2
        dir = np.asarray(self.direction)
        normal1 = np.asarray(self.position) - np.asarray(center)
        norm = LA.norm(normal1)
        normal = normal1/norm
        phi1 = acos(np.dot(-normal,dir))
        if phi1 < pi/2:
            c1 = np.dot(-normal,dir)
            c2 = 1 - ((n**2)*(1-c1**2))
            refracted_dir_cos1 = np.array([(n)*(dir[0]) +(n*c1-sqrt(c2))*normal[0], (n)*(dir[1]) +(n*c1-sqrt(c2))*normal[1], (n)*(dir[2]) +(n*c1-sqrt(c2))*normal[2]])
            norm1 = LA.norm(refracted_dir_cos1)
            refracted_dir_cos = refracted_dir_cos1/norm1        
        else:
            c3 = np.dot(normal,dir)
            c4 = 1 - ((n**2)*(1-c3**2))
            refracted_dir_cos2 = np.array([(n)*(dir[0]) - (n*c3-sqrt(c4))*normal[0], (n)*(dir[1]) - (n*c3-sqrt(c4))*normal[1], (n)*(dir[2]) -(n*c3-sqrt(c4))*normal[2]])
            norm3 = LA.norm(refracted_dir_cos2)
            refracted_dir_cos = refracted_dir_cos2/norm3  
            
        return refracted_dir_cos

    def power(self, I, eta, ex, refraction = False):
        if refraction==False:
            I_new = I*(1-eta)
        else:
            I_new = I*(1-ex)
        return I_new    
        



        
        
        
        
        
        
        
        