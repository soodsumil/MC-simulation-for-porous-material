# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 21:26:32 2019

@author: Dell
"""

import csv
import numpy as np

class Import_export_data:
        
    def import_data(self, name):
        f = open(name + '.csv','r')
        reader = csv.reader(f)
        data = []
        for row in reader:
            data.append(row)
        data = np.array(data,float)
        data = data.tolist()
        return data
    
    def export_data_row(self, name, data):
        with open(name + '.csv', "w", newline = '') as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerow(data)
            
    def export_data_rows(self, name, data):
        with open(name + '.csv', "w", newline = '') as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerows(data)
    
    def slicing_data(self, data, initial_dim, final_dim):
        sliced_data = []
        for a in data:
            if a[0]<final_dim[0] and a[0]>initial_dim[0] and  a[1]<final_dim[1] and a[1]>initial_dim[1] and a[2]<final_dim[2] and  a[2]>initial_dim[2]:
                sliced_data.append(a)
        return sliced_data






