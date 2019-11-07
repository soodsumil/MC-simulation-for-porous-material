# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:05:53 2019

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
from Class_ray import Ray
from Penetration_length_using_classes import Penetration_length
from Class_light_properties import Light_properties


class Monte_carlo:
    
    def monte_carlo_specular(self, iterations, diameter, absorption_coff, extinction_coff,refractive_index2, data, bed_dim, refractive_index1 = 1.003):
        D = diameter
        eta = absorption_coff
        ex = extinction_coff
        nta2 = refractive_index2
        nta1 = refractive_index1
        A = bed_dim
        E1 = []
        E2 = []
        pen_length = Penetration_length(data)
        pen_len = pen_length.penetration_length(D, iterations)
        
        for i in range(iterations):
            x = []
            dir_cos = []
            l = random.choice(pen_len)
            y = Divergence()
            dirc = y.cone_divergence()
            x.append([0,0,0])
            dir_cos.append(dirc)
            ep = dirc*l
            sphere_center = np.array([0,0,0])
            I = 1
            ii = 0
            iii = 0
            while I > 1E-10:
                overlap = False
                if len(x)>1:
                    l = random.choice(pen_len)
                    ep = np.array([x[ii][0]+dir_cos[ii][0]*l,x[ii][1]+dir_cos[ii][1]*l,x[ii][2]+dir_cos[ii][2]*l]) 
                if ep[2]>0 and ep[2] < A[2]:
                    y = Divergence()
                    center_dir_sec = y.cone_divergence()
                    rot = y.rotation(center_dir_sec)
                    center_dir = np.matmul(rot,center_dir_sec)
                    center_dir = np.reshape(center_dir,(1,3))
                    sphere_cen = ep + center_dir*D/2
                    sphere_cen = np.reshape(sphere_cen,(1,3))
                    sphere_center = np.vstack((sphere_center,sphere_cen))
                    for a in sphere_cen[1:]:
                        z = Overlap()
                        y = z.sphere_sphere(a,sphere_cen.tolist(), D)
                        if y == True:
                            overlap = True
                    for b in sphere_center[1:-2]:
                        r = Intersection(ep, dir_cos[ii])
                        eps = r.line_sphere_intersection(D,b)
                        if eps != False:
                            overlap = True
                            
                    if overlap==True:
                        iii = iii+1
                        sphere_center = np.delete(sphere_center,ii,0)
                        if iii > 100:
                            del x[ii]
                            del dir_cos[ii]
                            sphere_center = np.delete(sphere_center,ii-1,0)
                            ii = ii-1
                            iii = 0
                            continue
                        continue
                    
                    z = Light_properties(nta1,nta2)
                    rho = z.fresnels_law(ep, dir_cos[ii],sphere_center[ii+1].tolist())
                    rand = np.random.random(1)[0]
                    rho_avg = rho[0]
                    iii = 0
                    
                    if rand < rho_avg:
                        z = Ray(ep,dir_cos[ii])
                        reflected_dir = z.specular_reflection(sphere_center[ii+1].tolist())
                        dir_cos.append(reflected_dir)
                        x.append(ep)
                        I = z.power(I,eta,ex)
                        ii = ii+1
                        
                    else:
                        z = Ray(ep,dir_cos[ii])
                        refracted_dir_sec = z.refraction(D,sphere_center[ii+1].tolist(),nta1,nta2)
                        y = Intersection(ep,refracted_dir_sec)
                        eps = y.line_sphere_intersection(D,sphere_center[ii+1].tolist())
                        #s = sqrt((eps[0][1]-eps[1][0])**2 + (eps[0][1]-eps[1][1])**2 + (eps[0][2]-eps[1][2])**2)
                        q = Ray(eps[0],refracted_dir_sec)
                        refracted_dir = q.refraction(D,sphere_center[ii+1].tolist(),nta2,nta1)
                        dir_cos.append(refracted_dir)
                        x.append(eps[0])
                        I = z.power(I,eta,ex,True)
                        ii = ii+1
                        
                elif ep[2]< 0:
                    I = -I
                    E2.append(I)
                    break
                else:
                    E1.append(I)
                    break   
        E1_avg = sum(E1)/len(E1)
        E2_avg = sum(E2)/len(E2)                
        return [E1, E2,E1_avg, E2_avg]                        
                
                

       
    def monte_carlo_diffuse(self, iterations, diameter, absorption_coff, extinction_coff,refractive_index2, data, bed_dim, refractive_index1 = 1.003):
        D = diameter
        eta = absorption_coff
        ex = extinction_coff
        nta2 = refractive_index2
        nta1 = refractive_index1
        A = bed_dim
        E1 = []
        E2 = []
        pen_length = Penetration_length(data)
        pen_len = pen_length.penetration_length(D, iterations)
        
        for i in range(iterations):
            x = []
            dir_cos = []
            l = random.choice(pen_len)
            y = Divergence()
            dirc = y.cone_divergence()
            x.append([0,0,0])
            dir_cos.append(dirc)
            ep = dirc*l
            sphere_center = np.array([0,0,0])
            I = 1
            ii = 0
            iii = 0
            while I > 1E-10:
                overlap = False
                if len(x)>1:
                    l = random.choice(pen_len)
                    ep = np.array([x[ii][0]+dir_cos[ii][0]*l,x[ii][1]+dir_cos[ii][1]*l,x[ii][2]+dir_cos[ii][2]*l]) 
                if ep[2]>0 and ep[2] < A[2]:
                    y = Divergence()
                    center_dir_sec = y.cone_divergence()
                    rot = y.rotation(center_dir_sec)
                    center_dir = np.matmul(rot,center_dir_sec)
                    center_dir = np.reshape(center_dir,(1,3))
                    sphere_cen = ep + center_dir*D/2
                    sphere_cen = np.reshape(sphere_cen,(1,3))
                    sphere_center = np.vstack((sphere_center,sphere_cen))
                    for a in sphere_cen[1:]:
                        z = Overlap()
                        y = z.sphere_sphere(a,sphere_cen.tolist(), D)
                        if y == True:
                            overlap = True
                    for b in sphere_center[1:-2]:
                        r = Intersection(ep, dir_cos[ii])
                        eps = r.line_sphere_intersection(D,b)
                        if eps != False:
                            overlap = True
                            
                    if overlap==True:
                        iii = iii+1
                        sphere_center = np.delete(sphere_center,ii,0)
                        if iii > 100:
                            del x[ii]
                            del dir_cos[ii]
                            sphere_center = np.delete(sphere_center,ii-1,0)
                            ii = ii-1
                            iii = 0
                            continue
                        continue
                    
                    z = Light_properties(nta1,nta2)
                    rho = z.fresnels_law(ep, dir_cos[ii],sphere_center[ii+1].tolist())
                    rand = np.random.random(1)[0]
                    rho_avg = rho[0]
                    iii = 0
                    
                    if rand < rho_avg:
                        z = Ray(ep,dir_cos[ii])
                        reflected_dir = z.diffuse_reflection(sphere_center[ii+1].tolist())
                        reflected_dir = reflected_dir.tolist()
                        reflected_dir1 = []
                        for a in reflected_dir:
                            for b in a:
                                reflected_dir1.append(b)
                        dir_cos.append(reflected_dir1)
                        x.append(ep)
                        I = z.power(I,eta,ex)
                        ii = ii+1
                        
                    else:
                        z = Ray(ep,dir_cos[ii])
                        refracted_dir_sec = z.refraction(D,sphere_center[ii+1].tolist(),nta1,nta2)
                        y = Intersection(ep,refracted_dir_sec)
                        eps = y.line_sphere_intersection(D,sphere_center[ii+1].tolist())
                        #s = sqrt((eps[0][1]-eps[1][0])**2 + (eps[0][1]-eps[1][1])**2 + (eps[0][2]-eps[1][2])**2)
                        q = Ray(eps[0],refracted_dir_sec)
                        refracted_dir = q.refraction(D,sphere_center[ii+1].tolist(),nta2,nta1)
                        dir_cos.append(refracted_dir)
                        x.append(eps[0])
                        I = z.power(I,eta,ex,True)
                        ii = ii+1
                        
                elif ep[2]< 0:
                    I = -I
                    E2.append(I)
                    break
                else:
                    E1.append(I)
                    break   
        E1_avg = sum(E1)/len(E1)
        E2_avg = sum(E2)/len(E2)                
        return [E1, E2,E1_avg, E2_avg]  
            
                


