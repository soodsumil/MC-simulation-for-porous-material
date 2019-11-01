# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 20:14:50 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from math import sqrt, asin, pi, sin, cos

d = 0.00476
y = np.array([d/2,d/2,d/2])
x = np.array([0,0,0])
x = np.vstack((x,y))
c = np.array([])
total = 0

V = (4/24)*pi*d*d*d
Dir_cos = np.array([0,0,0])
z =np.array([1/sqrt(3),1/sqrt(3),1/sqrt(3)])
Dir_cos= np.vstack((Dir_cos,z))

p = (1 - V)/1
round = 1
while p>0.9:
    #print(1)
    b = []
    dis_list = []
    i1 = np.array([1,0,0])
    j1 = np.array([0,1,0])
    k1 = np.array([0,0,1])
    for i in range(3):
        L = x.size
        ii = 1
        #print(L+5)
        while x.size <= L:
            #print(2)
            sphere_considered = x[round]
            R_gamma = np.random.random(1)[0]
            R_theta = np.random.random(1)[0]
            gamma = asin(1 - 2*R_gamma)
            theta = 2*pi*R_theta
            k2 = [Dir_cos[round]]
            
            norm = np.cross(k2,k1)
            i2 = (np.cross(k2,k1))/LA.norm(norm)
            j2 = np.cross(k2,i2)
            rot = np.array([[i2],[j2],[k2]])
            rot = rot.transpose()

            dir = np.array([sin(gamma)*cos(theta),sin(gamma)*sin(theta),cos(gamma)])
            dir_grd = np.matmul(rot,dir)
            dir_grd = np.reshape(dir_grd, (1, 3))
            #print(dir_grd)
            
            X = np.array([x[round][0]+(dir_grd[0][0]*d),x[round][1]+(dir_grd[0][1]*d),x[round][2]+(dir_grd[0][2]*d)])
            #print(X)
            #print(L)
            for a in x:
                dis = (X[0] - a[0])**2 + (X[1] - a[1])**2 + (X[2] - a[2])**2       
                dis_list.append(dis)
            
            for a in dis_list:
                if a<(d**2):
                    b.append(1)
                else:
                    b.append(0)
    
            if X[0] < 1 and X[0]> 0 and X[1] > 0 and  X[1] < 1 and X[2] < 1 and X[2] > 0 and sum(b)< 0.5:  
                x = np.vstack((x,X))
                Dir_cos = np.vstack((Dir_cos,dir_grd))
                #print(x)
            ii = ii +1
            if ii > 5000:
                #print(1)
                break    
       
    i = len(x)
    p = (1-i*V)
    round = round+1


for a in x:
    m = []
    for b in x:
        
        distance = (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2
        if distance <= (0.00476**2):
            m.append(1)
    c = np.append(c,m)
    
        
print(total/len(x))
print(len(x))
print(len(c))
