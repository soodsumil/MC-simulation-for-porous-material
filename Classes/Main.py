# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:51:13 2019

@author: Dell
"""

from Monte_carlo_using_classes import Monte_carlo

f = Monte_carlo()
E1 = f.monte_carlo_specular(5000, 0.476,0.4,0.4,2,'my_data',[50,50,2])
       
f = Monte_carlo()
E2 = f.monte_carlo_diffuse(5000, 0.476,0.4,0.4,2,'my_data',[50,50,2])