# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 23:02:48 2019

@author: Dell
"""

import numpy as np
from numpy import linalg as LA
from Class_intersection import Intersection


class Overlap:
    def sphere_sphere(self,center1,center2,diameter1, diameter2=-1):
        if(diameter2==-1):
            diameter2 = diameter1
        dis = LA.norm((np.asarray(center1)- np.asarray(center2)))
        if (dis < (diameter1/2+diameter2/2)) :
            return True
        else:
            return False
        
    def line_sphere(self,position,direction,center,diameter):
        x = Intersection(position,direction)
        y = x.line_sphere_intersection(diameter,center)
        if y != False:
            return True
        else:
            return False
        
        