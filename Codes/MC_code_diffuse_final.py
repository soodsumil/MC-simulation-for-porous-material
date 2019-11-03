# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 22:33:19 2019

@author: Dell
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 18:38:55 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from math import sqrt, asin,acos, pi, sin, cos, tan, exp
import random


D = 0.476 #sphere diameter(cm)
eta = 0.4 #absorbtivity
nta = 20000 #refractive index
ex = 0.3 #extinction coefficient
nta_air = 1.0003
A = np.array([50,50,2]) #bed dimentions(cm)
E1 = [] #save intencities that gets transmitted
E2 = [] # saves intencities that are reflected

for i in range(5000): #number of iterations
#ground frame
    i1 = np.array([1,0,0])
    j1 = np.array([0,1,0])
    k1 = np.array([0,0,1])
    o = np.array([0,0,0])# origin
#List of the points being hit by this ray 
    x = np.array(o)
#direction cosines of the ray
    l = random.choice(penetration_len) #select a random penetration length from the data
    #l = (1+np.random.random(1)[0])*D
    R_gamma = np.random.random(1)[0]
    R_theta = np.random.random(1)[0]
    gamma = asin(R_gamma**0.5)
    theta = 2*pi*R_theta
    dir_cos = np.array([sin(gamma)*cos(theta),sin(gamma)*sin(theta),cos(gamma)])
    Dir_cosines = np.array([0,0,0])
    Dir_cosines = np.vstack((Dir_cosines,dir_cos)) #array to save all the direction cosines of the ray
    ip = np.array([dir_cos[0]*l,dir_cos[1]*l,dir_cos[2]*l]) # first hit point
    x = np.vstack((x,o)) #first point in x is redundant and is added just to make cade easier
    ep = ip # ep is the end point of each step
    sphere_center=np.array([0,0,0]) 
    I = 1 #Intensity or energy of the ray
    ii = 1  
    iii= 0
    
    while I > 1E-10:
        overlap = False
# skipping the first step when end point is already known       
        if len(x)> 2:
            l = random.choice(penetration_len)
            #l = (0.5+np.random.random(1)[0])*D #in case we don't have data for penetration length
            ep = np.array([x[ii][0]+Dir_cosines[ii][0]*l,x[ii][1]+Dir_cosines[ii][1]*l,x[ii][2]+Dir_cosines[ii][2]*l]) #end point of the ray for this iteration

        if ep[2]< A[2] and ep[2]>0:
#random cone and circumferential angles to get a random point of the selected particle
            R_gamma1 = np.random.random(1)[0] #randomly getting the center of the sphere being hit
            R_theta1 = np.random.random(1)[0]
            gamma1 = acos(R_gamma1)
            theta1 = 2*pi*R_theta1
# Rotation matrix to get the direction in the ground frame
            k2_ = np.array([(ep[0]-x[ii][0])/(D/2),(ep[1]-x[ii][1])/(D/2),(ep[2]-x[ii][2])/(D/2)])
            norm1 =LA.norm(k2_)
            k2 = k2_/norm1
            norm2 = np.cross(k2,k1)
            i2 = (np.cross(k2,k1))/LA.norm(norm2)
            j2 = np.cross(k2,i2)
            rot = np.array([[i2],[j2],[k2]])#rotation matrix to get direction cosines in the ground frame
            rot = rot.transpose()
#Direction cosines of the center of the sphere from the end point in grd frame
            dir_cos_sec = np.array([sin(gamma1)*cos(theta1), sin(gamma1)*sin(theta1), cos(gamma1)])
            Dir_cos_grd = np.matmul(rot,dir_cos_sec)
#coordinates of the center of the spheres
            sphere_center_grd = np.array([Dir_cos_grd[0]*D/2 + ep[0],Dir_cos_grd[1]*D/2 + ep[1],Dir_cos_grd[2]*D/2 + ep[2]])     
            sphere_center_grd = np.reshape(sphere_center_grd, (1, 3))
            sphere_center = np.vstack((sphere_center,sphere_center_grd))# adding sphere being hit to the list to keep track of the spheres being added, later add this to the main list
    
#code for checking the overlap
            #Code for checking whether the new sphere is overlaping with an existing sphere 
            dis_list = []
            B = []
            for a in sphere_center[1:-1]:
                dis = (sphere_center_grd[0][0] - a[0])**2 + (sphere_center_grd[0][1] - a[1])**2 + (sphere_center_grd[0][2] - a[2])**2       
                dis_list.append(dis)
            for a in dis_list:
                if a<((D**2)-0.001):
                    B.append(1)
                else:
                    B.append(0)
            if sum(B)>0.5:
                overlap = True
            # Code to check whether the emitted ray is hitting any of the existing spheres
            p11= np.array([ep[0],ep[1],ep[2]])
            p22= np.array([ep[0]+Dir_cosines[ii][0]*l,ep[1]+Dir_cosines[ii][1]*l,ep[2]+Dir_cosines[ii][2]*l])
            p11 = np.reshape(p11,(1,3))
            p22 = np.reshape(p22,(1,3))
            
            for a in sphere_center[1:-2]:
                a = np.array(a)
                y = np.divide(p22 - p11, np.linalg.norm(p22 - p11))                
                n1_ = np.squeeze(np.asarray(p11-a))
                n1 = np.array([n1_[0] - a[0],n1_[1] - a[1],n1_[2] - a[2]])
                n2 = np.squeeze(np.asarray(y))
                n3 = np.squeeze(np.asarray(a-p22))
                s = np.dot(n1, n2)
                t = np.dot(n3, n2)
                h = np.maximum.reduce([s, t, 0])
                c = np.cross(a - p11, y)
                w = LA.norm(c)
                Dis_center_to_line = sqrt(w**2 + h**2)
                if Dis_center_to_line < D/2:
                    overlap = True
#If any of the above conditions is true then do the tracing again and so deleting the previous step      
            if overlap==True:
                iii = iii+1
                sphere_center = np.delete(sphere_center,ii,0)
                if iii > 100:
                    x = np.delete(x,ii,0)
                    Dir_cosines = np.delete(Dir_cosines,ii,0)
                    sphere_center = np.delete(sphere_center,ii-1,0)
                    ii = ii-1
                    iii = 0
                    continue
                continue

#code for checking whether to do reflection or refraction
            normal = np.array([ep[0]-sphere_center[ii][0],ep[1]-sphere_center[ii][1],ep[2]-sphere_center[ii][2]])
            norm = LA.norm(normal)
            normal_dir_cos = np.array([normal[0]/norm, normal[1]/norm, normal[2]/norm])

            phi = acos(np.dot(-normal_dir_cos,Dir_cosines[ii]))
            phi2 = asin(nta_air*sin(phi)/nta)
            
            rho_parallel = (nta*cos(phi)-nta_air*cos(phi2))/(nta*cos(phi)+nta_air*cos(phi2))
            rho_perpendicular = (nta_air*cos(phi)-nta*cos(phi2))/(nta_air*cos(phi)+nta*cos(phi2))

            trans_parallel = (2*sin(phi2)*cos(phi))/(sin(phi+phi2)*cos(phi-phi2))
            trans_perpendicular = (2*sin(phi2)*cos(phi))/(sin(phi+phi2))
        
            rho_avg = (rho_parallel**2 + rho_perpendicular**2)/2
            trans_avg = (trans_parallel**2 + trans_perpendicular**2)/2
            
            rand1 = np.random.random(1)[0]

            iii = 0
            if rand1< rho_avg: # Diffuse reflection
                k3 = normal_dir_cos
                norm3 = np.cross(k3,k1)
                i3 = (np.cross(k3,k1))/LA.norm(norm2)
                j3 = np.cross(k3,i3)
                rot3 = np.array([[i2],[j2],[k2]])#rotation matrix to get direction cosines in the ground frame
                rot3 = rot.transpose()
                R_gamma2 = np.random.random(1)[0] #getting the center of the sphere being hit
                R_theta2 = np.random.random(1)[0]
                gamma2 = asin(1 - 2*R_gamma2)
                theta2 = 2*pi*R_theta2
                reflected_dir_cos1 = np.array([sin(gamma2)*cos(theta2), sin(gamma2)*sin(theta2), cos(gamma2)])
                reflected_dir_cos = np.matmul(rot,reflected_dir_cos1)
                reflected_dir_cos = np.reshape(reflected_dir_cos, (1, 3))
                Dir_cosines = np.vstack((Dir_cosines,reflected_dir_cos))
                x = np.vstack((x,ep))
                ii = ii+1
                I = I - I*eta # loss of intensity due to collision
            
            else: #refraction
                n = nta_air/nta
                cphi = np.dot(-normal_dir_cos,Dir_cosines[ii])
                phi = acos(np.dot(-normal_dir_cos,Dir_cosines[ii]))
                phi2 = asin(n*sin(phi))

                c1 = np.dot(-normal_dir_cos,Dir_cosines[ii])
                c2 = 1 - ((n**2)*(1-cphi**2))
                # refracted direction cosines
                refracted_dir_cos1 = np.array([(n)*(Dir_cosines[ii][0]) +(n*cphi-sqrt(c2))*normal_dir_cos[0], (n)*(Dir_cosines[ii][1]) +(n*cphi-sqrt(c2))*normal_dir_cos[1], (n)*(Dir_cosines[ii][2]) +(n*cphi-sqrt(c2))*normal_dir_cos[2]])
                norm1 = LA.norm(refracted_dir_cos1)
                refracted_dir_cos = np.array([refracted_dir_cos1[0]/norm1,refracted_dir_cos1[1]/norm1,refracted_dir_cos1[2]/norm1])
                
                # find the final point after the ray has refracted and travelled through the sphere
                p1 = ep
                p2 = np.array([p1[0]+refracted_dir_cos[0]*1,p1[1]+refracted_dir_cos[1]*1,p1[2]+refracted_dir_cos[2]*1])
                u = (p2[0] - p1[0])**2+(p2[1] - p1[1])**2+(p2[2] - p1[2])**2
                v = -2*((p2[0] - p1[0])*(sphere_center[ii][0]-ep[0])) + -2*((p2[1] - p1[1])*(sphere_center[ii][1]-ep[1]))+ -2*((p2[2] - p1[2])*(sphere_center[ii][2]-ep[2]))
                w = (sphere_center[ii][0]-ep[0])**2+(sphere_center[ii][1]-ep[1])**2+(sphere_center[ii][2]-ep[2])**2 - (D/2)**2
                t1 = (-v + sqrt(v**2-(4*u*w)))/(2*u)
                ep2 = np.array([ep[0]+ t1*refracted_dir_cos[0], ep[1]+ t1*refracted_dir_cos[1], ep[2]+ t1*refracted_dir_cos[2]])
                s = sqrt((ep2[0]-ep[0])**2 + (ep2[1]-ep[1])**2+ (ep2[2]-ep[2])**2)
                normal2 = np.array([ep2[0]-sphere_center[ii][0],ep2[1]-sphere_center[ii][1],ep2[2]-sphere_center[ii][2]])
                norm2 = LA.norm(normal2)
                normal_dir_cos2 = np.array([normal2[0]/norm2, normal2[1]/norm2, normal2[2]/norm2])
                b1 = nta/nta_air
                phi3 = acos(np.dot(normal_dir_cos2,refracted_dir_cos))
                phi4 = asin(b1*sin(phi3))
                print(phi4)                         
                c3 = np.dot(normal_dir_cos2,refracted_dir_cos)
                c4 = 1 - ((b1**2)*(1-c3**2))
                refracted_dir_cos2 = np.array([(b1)*(refracted_dir_cos[0]) - (b1*c3-sqrt(c4))*normal_dir_cos2[0], (b1)*(refracted_dir_cos[1]) - (b1*c3-sqrt(c4))*normal_dir_cos2[1], (b1)*(refracted_dir_cos[2]) -(b1*c3-sqrt(c4))*normal_dir_cos2[2]])
                norm3 = LA.norm(refracted_dir_cos2)
                refracted_dir_cos3 = np.array([refracted_dir_cos2[0]/norm3,refracted_dir_cos2[1]/norm3,refracted_dir_cos2[2]/norm3])
                print(acos(np.dot(normal_dir_cos2,refracted_dir_cos3)))
                #adding the final values to the list
                x = np.vstack((x,ep2))
                Dir_cosines = np.vstack((Dir_cosines,refracted_dir_cos3))
                ii = ii+1
                I = I*rho_avg
                
        elif ep[2]< 0:
            I = -I
            E2.append(I)
            #print("end2")
            break
        else:
            E1.append(I)
            #print("end1")
            break
    
    
    
    
