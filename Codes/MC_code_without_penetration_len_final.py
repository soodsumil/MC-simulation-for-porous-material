# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 01:47:53 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from math import sqrt, asin,acos, pi, sin, cos, exp
import random
import csv 

X = xx1 # xx1 is the data set with the origins of the all the particles in ground frame
particle_data = X[:] 
D = 0.476 #sphere diameter(cm)
eta = 0.4 #absorbtivity
nta = 2 #refractive index
ex = 0.3 #extinction coefficient
nta_air = 1.0003
dim = np.array([[10,10,0],[20,20,2]]) #bed dimentions start to end in the groung frame(cm)
#X = []#points that a all the rays hit
E1 = []  # energy that is transmitted through the bed
E2 = [] # energy being reflected back
E3 = [] # energy going out through the walls
for i in range(5000):
#initial ray
    # random initial point for sending the ray in the particle bed 
    sp_x = 13 +((random.random())*(dim[1][0]-15)) 
    sp_y = 13 +((random.random())*(dim[1][1]-15))
    sp = np.array([sp_x,sp_y,-5])
    dir_cos_ini = np.array([0,0,1])
    # Direction cosines of the ground frame
    i1 = np.array([1,0,0])
    j1 = np.array([0,1,0])
    k1 = np.array([0,0,1])
    #pts is a list to add all the points being hit by the ray and dir_cos is the list to save direction cosines
    pts = np.array([0,0,0])
    pts = np.vstack((pts,sp))
    dir_cos = np.array([0,0,0])
    dir_cos = np.vstack((dir_cos,dir_cos_ini))
    sphere_centers_hit = np.array([[0,0,0],[0,0,0]]) # all the spheres being hit by the ray
    I = 1
    ii = 1 # inicializing for use in the code
    
    while I>  1E-10:
        X = xx1 # list of all particles
#finding the sphere hit by this ray given the initial point and the direction cosines
        p1= np.array([pts[ii][0],pts[ii][1],pts[ii][2]]) #initial point outside the bed
        p2= np.array([pts[ii][0]+dir_cos[ii][0]*10,pts[ii][1]+dir_cos[ii][1]*10,pts[ii][2]+dir_cos[ii][2]*10]) # second point taken randomly on the ray at a distance of 10cm from initial point
        p1 = np.reshape(p1,(1,3))# reshape to make it compatible for further calculations
        p2 = np.reshape(p2,(1,3))
        sphere_hitting = np.array([[0,0,0],[0,0,0]]) # list being inicialized at each step to save all the spheres being hit by this ray so that later we can calculate the the first(nearest) hit point and sphere 
        particle_data_new = X[:] # data from which the sphere hit from lst iteration will be deleted to have proper further calculations
        if ii >1: # supressing first iteration
            index1,index2= np.where(X == c)  # finding the index of last hit sphere
            del particle_data_new[index1[1]] # deleting the last hit sphere
        # code for finding all the spheres being hit by finding the distance of the centers from this ray and the saving the particles with distance less than D/2(radius)
        for a in particle_data_new[1:]:
            b = np.array(a)
            y = np.divide(p2 - p1, np.linalg.norm(p2 - p1))
            n1 = np.squeeze(np.asarray(p1-b))
            n2 = np.squeeze(np.asarray(y))
            n3 = np.squeeze(np.asarray(b-p2))
            s = np.dot(n1, n2)
            t = np.dot(n3, n2)
            h = np.maximum.reduce([s, t, 0])
            c = np.cross(b - p1, y)
            w = LA.norm(c)
            Dis_center_to_line = sqrt(w**2 + h**2)
            if Dis_center_to_line < D/2:
                sphere_hitting= np.vstack((sphere_hitting,a))
        #checking if the ray is not hitting any particles then it has left the bed and thus checking which boundary did the ray leave
        if len(pts)>2:
            p3 = np.array([pts[ii][0]+dir_cos[ii][0]*D,pts[ii][1]+dir_cos[ii][1]*D,pts[ii][2]+dir_cos[ii][2]*D])
            if len(sphere_hitting)<3:
                if p3[2]<dim[0][2]: # reflection
                    E2.append(I)
                if p3[2]>dim[1][2]: # transmitted
                    E1.append(I)
                if p3[0]<dim[0][0] or p3[0]>dim[1][0] or p3[1]<dim[0][1] or p3[1]> dim[1][1]: # passing through the side boundries
                    E3.append(I)
                break
# code to find the closest sphere hit by the ray       
        y = 100# random initialization         
        for b in sphere_hitting[2:]:
            dis = sqrt((pts[ii][0] - b[0])**2+(pts[ii][1] - b[1])**2+(pts[ii][2] - b[2])**2)
            if dis< y:
                y = dis
                sphere_being_hit_center = b
       
        c = sphere_being_hit_center # saving the sphere being hit in a variable c 
        sphere_centers_hit = np.vstack((sphere_centers_hit,c))
#Code for finding the points on this spphere that are hit by the ray and then selecting the first hit point
        p1 = np.squeeze(np.asarray(p1))
        p2 = np.squeeze(np.asarray(p2))
        u = (p2[0] - p1[0])**2+(p2[1] - p1[1])**2+(p2[2] - p1[2])**2
        v = 2*(((p2[0] - p1[0])*(-sphere_being_hit_center[0]+p1[0])) + ((p2[1] - p1[1])*(-sphere_being_hit_center[1]+p1[1]))+ ((p2[2] - p1[2])*(-sphere_being_hit_center[2]+p1[2])))
        w = (sphere_being_hit_center[0]-p1[0])**2 + (sphere_being_hit_center[1]-p1[1])**2+ (sphere_being_hit_center[2]-p1[2])**2- (D/2)**2
        t1 = (-v + sqrt((v**2)-4*u*w))/(2*u)
        t2 = (-v - sqrt((v**2)-4*u*w))/(2*u)
        ep1 = np.array([p1[0]+(p2[0]-p1[0])*t1,p1[1]+(p2[1]-p1[1])*t1,p1[2]+(p2[2]-p1[2])*t1])
        ep2 = np.array([p1[0]+(p2[0]-p1[0])*t2,p1[1]+(p2[1]-p1[1])*t2,p1[2]+(p2[2]-p1[2])*t2])
        dis1 = sqrt((p1[0]-ep1[0])**2 + (p1[1]-ep1[1])**2 + (p1[2]-ep1[2])**2)
        dis2 = sqrt((p1[0]-ep2[0])**2 + (p1[1]-ep2[1])**2 + (p1[2]-ep2[2])**2)
        if (dis1 <= dis2):
            ep = ep1
        else:
            ep = ep2
#Code for finding whether tha ray is getting refracted or reflected
        if ep[2]< dim[1][2] and ep[2]>-0.25:
            normal = np.array([ep[0]-c[0],ep[1]-c[1],ep[2]-c[2]])
            norm = LA.norm(normal)
            normal_dir_cos = np.array([normal[0]/norm, normal[1]/norm, normal[2]/norm])
            phi = acos(np.dot(-normal_dir_cos,dir_cos[ii]))
            phi2 = asin(nta_air*sin(phi)/nta)
            rho_parallel = (nta*cos(phi)-nta_air*cos(phi2))/(nta*cos(phi)+nta_air*cos(phi2))
            rho_perpendicular = (nta_air*cos(phi)-nta*cos(phi2))/(nta_air*cos(phi)+nta*cos(phi2))
            trans_parallel = (2*sin(phi2)*cos(phi))/(sin(phi+phi2)*cos(phi-phi2))
            trans_perpendicular = (2*sin(phi2)*cos(phi))/(sin(phi+phi2))
            rho_avg = (rho_parallel**2 + rho_perpendicular**2)/2
            trans_avg = (trans_parallel**2 + trans_perpendicular**2)/2    
            rand1 = np.random.random(1)[0]
            
            if rand1< rho_avg: # spectral reflection
                dot_prod = np.dot(normal_dir_cos,dir_cos[ii])
                reflected_dir_cos1 = np.array([(dir_cos[ii][0] - 2*normal_dir_cos[0]*dot_prod), (dir_cos[ii][1] - 2*normal_dir_cos[1]*dot_prod), (dir_cos[ii][2] - 2*normal_dir_cos[2]*dot_prod)])
                norm1 = LA.norm(reflected_dir_cos1)
                reflected_dir_cos = np.array([reflected_dir_cos1[0]/norm1,reflected_dir_cos1[1]/norm1,reflected_dir_cos1[2]/norm1])
                angel_inci_reflect = acos(np.dot(reflected_dir_cos,normal_dir_cos))
                pts = np.vstack((pts,ep))
                dir_cos = np.vstack((dir_cos,reflected_dir_cos))
                ii = ii+1
                I = I - I*eta
    
            else: #refraction
                n = nta_air/nta
                cphi = np.dot(-normal_dir_cos,dir_cos[ii])
                phi = acos(np.dot(-normal_dir_cos,dir_cos[ii]))
                phi2 = asin(n*sin(phi))
                c1 = np.dot(-normal_dir_cos,dir_cos[ii])
                c2 = 1 - ((n**2)*(1-cphi**2))
                # refracted direction cosines
                refracted_dir_cos1 = np.array([(n)*(dir_cos[ii][0]) +(n*cphi-sqrt(c2))*normal_dir_cos[0], (n)*(dir_cos[ii][1]) +(n*cphi-sqrt(c2))*normal_dir_cos[1], (n)*(dir_cos[ii][2]) +(n*cphi-sqrt(c2))*normal_dir_cos[2]])
                norm1 = LA.norm(refracted_dir_cos1)
                refracted_dir_cos = np.array([refracted_dir_cos1[0]/norm1,refracted_dir_cos1[1]/norm1,refracted_dir_cos1[2]/norm1])
                 # find the final point after the ray has refracted and travelled through the sphere
                p1 = ep
                p2 = np.array([p1[0]+refracted_dir_cos[0]*1,p1[1]+refracted_dir_cos[1]*1,p1[2]+refracted_dir_cos[2]*1])
                u = (p2[0] - p1[0])**2+(p2[1] - p1[1])**2+(p2[2] - p1[2])**2
                v = -2*((p2[0] - p1[0])*(c[0]-ep[0])) + -2*((p2[1] - p1[1])*(c[1]-ep[1]))+ -2*((p2[2] - p1[2])*(c[2]-ep[2]))
                w = (c[0]-ep[0])**2+(c[1]-ep[1])**2+(c[2]-ep[2])**2 - (D/2)**2
                t1 = (-v + sqrt(v**2-(4*u*w)))/(2*u)
                ep2 = np.array([ep[0]+ t1*refracted_dir_cos[0], ep[1]+ t1*refracted_dir_cos[1], ep[2]+ t1*refracted_dir_cos[2]])
                s = sqrt((ep2[0]-ep[0])**2 + (ep2[1]-ep[1])**2+ (ep2[2]-ep[2])**2)
                pts = np.vstack((pts,ep2))
                dir_cos = np.vstack((dir_cos,refracted_dir_cos))
                ii = ii+1
                I = I - I*eta # need to decide how much enerfy gets lost 
        elif ep[2]< 0:
            I = -I
            E2.append(I)
            #print("end2")
            break
        else:
            E1.append(I)
            #print("end1")
            break
    
    
    
    
    
    
    