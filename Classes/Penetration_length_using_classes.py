# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 01:41:22 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from math import sqrt, asin,acos, pi, sin, cos, exp
import random
from Class_Import_export_data import Import_export_data
from Class_divergence import Divergence
from Class_overlap import Overlap
from Class_intersection import Intersection

class Penetration_length:
    def __init__(self,data):
        self.data= data    
    
    def penetration_length(self, diameter, iterations, ini_dim = -1, fin_dim = -1):
        x = Import_export_data()
        data1 = x.import_data(self.data)
        if ini_dim == -1:
            x1 = 500
            y1 = 500
            z1=500
            for a in data1:
                if a[0]<x1:
                    x1 = a[0]
                if a[1]<y1:
                    y1 = a[1]
                if a[2]<z1:
                   z1 = a[2]
            ini_dim = [x1,y1,z1]
        if fin_dim == -1:
            x1 = -500
            y1 = -500
            z1=-500
            for a in data1:
                if a[0]>x1:
                    x1 = a[0]
                if a[1]>y1:
                    y1 = a[1]
                if a[2]>z1:
                   z1 = a[2]
            fin_dim = [x1,y1,z1]
        
        data2 = x.slicing_data(data1, ini_dim, fin_dim)
        penetration_len = []
        
        for i in range(iterations):
            sphere_center = data2[:]
            d = False
            while d == False:
                g = random.choice(range(len(data2)-1))
                center = data2[g]
                if center[0]<(ini_dim[0] +diameter*5) or center[0]>(fin_dim[0] -diameter*5) or center[1]<(ini_dim[1] +diameter*5) or center[1]>(fin_dim[1] -diameter*5) or center[2]<(ini_dim[2] +diameter*1) or center[2]>(fin_dim[2] -diameter*1):
                    continue
                d = True
            
            del sphere_center[g]
            y = Divergence()
            dirc = y.circular_divergence()
            ep = np.asarray(center)+ np.asarray(dirc)*diameter
            rot = y.rotation(dirc)
            dir_ray_sec = y.cone_divergence()
            dir_ray = np.matmul(rot,dir_ray_sec)
            dir_ray = np.reshape(dir_ray,(1,3))
            ep = np.reshape(ep,(1,3))
            ep = np.squeeze(ep)
            dir_ray = np.squeeze(dir_ray)
            sphere_hitting = []
            for a in sphere_center:
                z = Overlap()
                hit = z.line_sphere(ep, dir_ray, a, diameter)
                if hit == True:
                    sphere_hitting.append(a)
            xx = 100
            if len(sphere_hitting) == 0:
                continue
            for b in sphere_hitting:
                dis = sqrt((center[0] - b[0])**2 + (center[1] - b[1])**2 + (center[2] - b[2])**2)
                if dis< xx:
                    xx = dis
                    closest_sphere = b
            r = Intersection(ep,dir_ray)
            eps = r.line_sphere_intersection(diameter, closest_sphere)
            penetration_length1 = sqrt((ep[0]-eps[0][0])**2 + (ep[1]-eps[0][1])**2 + (ep[2]-eps[0][2])**2)
            penetration_length2 = sqrt((ep[0]-eps[1][0])**2 + (ep[1]-eps[1][1])**2 + (ep[2]-eps[1][2])**2)
            if (penetration_length1 <= penetration_length2):
                penetration_length = penetration_length1
            else:
                penetration_length = penetration_length2
    
            penetration_len.append(penetration_length)
        
        return penetration_len
            
        
        
        
    
        
        
        







        
        

