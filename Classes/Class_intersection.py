# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 16:12:01 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from math import sqrt

class Intersection:
    def __init__(self, position, direction):
        self.position = position
        self.direction = direction
    
    def dis_to_point(self, length, center):
        p22 = np.asarray(self.position) + np.asarray(self.direction)*length
        p11 = np.asarray(self.position)
        a = np.asarray(center)
        n1 = p11-a
        n2 = np.divide(p22 - p11, np.linalg.norm(p22 - p11))                
        n3 = a-p22
        s = np.dot(n1, n2)
        t = np.dot(n3, n2)
        h = np.maximum.reduce([s, t, 0])
        c = np.cross(a - p11, n2)
        w = LA.norm(c)
        return sqrt(w**2 + h**2)
    
    def line_sphere_intersection(self, diameter, center):
        p2 = np.asarray(self.position) + np.asarray(self.direction)*10
        p1 = np.asarray(self.position)
        center = np.asarray(center)
        D = diameter
        u = (p2[0] - p1[0])**2+(p2[1] - p1[1])**2+(p2[2] - p1[2])**2
        v = 2*(((p2[0] - p1[0])*(-center[0]+p1[0])) + ((p2[1] - p1[1])*(-center[1]+p1[1]))+ ((p2[2] - p1[2])*(-center[2]+p1[2])))
        w = (center[0]-p1[0])**2 + (center[1]-p1[1])**2+ (center[2]-p1[2])**2- (D/2)**2
        if (v**2-4*u*w)<0:
            return False
        else:
            t1 = (-v + sqrt((v**2)-4*u*w))/(2*u)
            t2 = (-v - sqrt((v**2)-4*u*w))/(2*u)
            ep1 = np.array([p1[0]+(p2[0]-p1[0])*t1,p1[1]+(p2[1]-p1[1])*t1,p1[2]+(p2[2]-p1[2])*t1])
            ep2 = np.array([p1[0]+(p2[0]-p1[0])*t2,p1[1]+(p2[1]-p1[1])*t2,p1[2]+(p2[2]-p1[2])*t2])
            return [ep1,ep2]
        
        
    def line_plane_intersection(self, distance, point, normal):
        p1 = self.position
        p2 = np.asarray(self.position) + np.asarray(self.direction)*distance
        dot_prod = np.dot((p2-p1), np.asarray(normal))
        if dot_prod == 0:
            return 'parallel'
        else:
            s = (np.dot( np.asarray(normal), (np.asarray(point)-p1)))/dot_prod
            if s<1 and s>0:
                return (np.asarray(self.position) + (p2-p1)*s)
            else:
                return False
            
            




    
    






   
    
