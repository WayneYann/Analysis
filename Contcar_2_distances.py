# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:58:40 2016

@author: Josh
"""
from __future__ import division
import os
import numpy as np
import Text_Parse

directory=os.path.expanduser('~/Box Sync/Synced_Files\Coding\Research\Analysis\VASP_Files\O_Hollow')
full_file_paths = Text_Parse.get_filepaths(directory)
CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
Surface_Names = Text_Parse.between_values(CONTCAR_FILES,"\\O_Hollow\\","\\CONTCAR")

#labels1212 = np.genfromtxt('./Lansford_Stats_HW3_data/12.12.csv', delimiter=',', dtype=str)[0,:]
#data1212 = np.genfromtxt('C:\Users\Josh\Desktop\XSD files', delimiter=',')[1:,:]
Surface_Distances = []
num_Files = len(CONTCAR_FILES)
for i in range(0,num_Files):
    Multiplier = Text_Parse.file2listfixed(CONTCAR_FILES[i],2,4,0,None)
    Multiplier = Text_Parse.list_split(Multiplier," ")
    Multiplier = Text_Parse.string2float(Multiplier)
    Multiplier = np.array(Multiplier)
    Coordinates = Text_Parse.file2listflex(CONTCAR_FILES[i],71,['T','F'],0,59)
    Coordinates = Text_Parse.list_split(Coordinates," ")
    Coordinates = Text_Parse.string2float(Coordinates)
    Coordinates = np.array(Coordinates)
    Coordinates = np.dot(Coordinates,Multiplier)
    numAtoms = len(Coordinates)
    Oatom = Coordinates[numAtoms-1]
    mindistance = 100
    for j in range(0,numAtoms-1):
        distance = sum((Oatom-Coordinates[j])**2)**(0.5)
        if distance <= mindistance:
            mindistance = distance
            horizontal = ((Oatom[0]-Coordinates[j][0])**2+(Oatom[1]-Coordinates[j][1])**2)**(0.5)
            vertical = abs(Oatom[2]-Coordinates[j][2])
            #index = j
    #surface_distances = [Surface_Names[i], vertical/mindistance, horizontal/mindistance, vertical, horizontal, mindistance, index]
    surface_distances = [Surface_Names[i], vertical/mindistance, horizontal/mindistance]    
    Surface_Distances.append(surface_distances)