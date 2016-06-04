# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:58:40 2016

@author: Josh
"""
from __future__ import division
import os
import Text_Parse
from ase.io import read
import csv



directory=os.path.expanduser('~/Box Sync/Synced_Files\Coding\Research\Analysis\VASP_Files\N_Hollow')
full_file_paths = Text_Parse.get_filepaths(directory)
CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
CONTCAR_FILES = [i for i in CONTCAR_FILES if not ('freq' in i)]
Surface_Names = Text_Parse.between_values(CONTCAR_FILES,"\\N_Hollow\\","\\CONTCAR")

#Repeating surface
a=2
b=2

Distances = [['Surface', 'Distance', 'Veritical', 'Horizontal','cos(alpha)','sin(alpha)', 'NN2','NN2/NN','NN2cos(alpha)','NN2sin(alpha)']]
num_Files = len(CONTCAR_FILES)
for i in range(0,num_Files):
    ASE_CONTCAR = read(CONTCAR_FILES[i])
    ASE_CONTCAR.set_constraint()
    ASE_CONTCAR = ASE_CONTCAR.repeat((a,b,1))
    numAtoms = len(ASE_CONTCAR)
    Coordinates = ASE_CONTCAR.positions
    Oatom = ASE_CONTCAR.positions[numAtoms-1,:]
    mindistance = 1000
    NN2 = 1000
    for j in range(0,numAtoms-1):
        distance = sum((Oatom-Coordinates[j])**2)**(0.5)
        if distance <= mindistance:
            mindistance = distance
            horizontal = ((Oatom[0]-Coordinates[j][0])**2+(Oatom[1]-Coordinates[j][1])**2)**(0.5)
            vertical = abs(Oatom[2]-Coordinates[j][2])
    for j in range(0,numAtoms-1):
        distance = sum((Oatom-Coordinates[j])**2)**(0.5)
        if distance < NN2 and vertical*1.2 < abs(Oatom[2]-Coordinates[j][2]):
            NN2 = distance
            NN2vertical = abs(Oatom[2]-Coordinates[j][2])
            NN2horizontal = ((Oatom[0]-Coordinates[j][0])**2+(Oatom[1]-Coordinates[j][1])**2)**(0.5)
            #index = j
    #surface_distances = [Surface_Names[i], vertical/mindistance, horizontal/mindistance, vertical, horizontal, mindistance, index]
    surface_distances = [Surface_Names[i], mindistance, NN2/mindistances,vertical/mindistance,horizontal/mindistance,NN2,NN2vertical,NN2horizontal,NN2vertical/NN2,NN2horizontal/NN2]    
    Distances.append(surface_distances)
        
            
myfile = open('N_Distances.csv', 'wb')
wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
for i in range(0,len(Distances)):
    wr.writerow(Distances[i])
myfile.close()