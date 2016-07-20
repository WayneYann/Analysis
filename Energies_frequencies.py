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
import numpy as np



directory=os.path.expanduser('~/Box Sync/Synced_Files\Coding\Research\VASP_Files\VASP_Ad_v02')
full_file_paths = Text_Parse.get_filepaths(directory)
CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
OSZICAR_FILES = Text_Parse.list_contains(full_file_paths,"OSZICAR")
freq_FILES = Text_Parse.list_contains(full_file_paths,"freq")
CONTCAR_FILES = [i for i in CONTCAR_FILES if not ('Free_Surfaces' in i)]
OSZICAR_FILES = [i for i in OSZICAR_FILES if not ('Free_Surfaces' in i)]
OSZICAR_FILES = [i for i in OSZICAR_FILES if not ('Free_Surfaces' in i)]
Surface_Names = Text_Parse.between_values(CONTCAR_FILES,"\\VASP_Ad_v02\\","_CONTCAR")
Free_Surfaces = Text_Parse.list_contains(full_file_paths,"Free_Surfaces")
Free_Surfaces = [i for i in Free_Surfaces if ('OSZICAR' in i)]
#Repeating surface
repeatX=2
repeatY=2

Data = [['Surface', 'Substrate','Substrate Mass', 'NN','NN2','Adsorbate','Adsorbate mass','cos(alpha)','sin(alpha)', 'Mr','Zfrequency','frequency_corrected','Energy','Surface_Energy','Eads']]
num_Files = len(Surface_Names)
for i in range(0,num_Files):
    ASE_CONTCAR = read(CONTCAR_FILES[i])
    numAtoms = len(ASE_CONTCAR)
    numAdsorbates = numAtoms-64
    if CONTCAR_FILES[i].find('111') > 0:
        Surface = '111'
    elif CONTCAR_FILES[i].find('110') > 0:
        Surface = '110'
    elif CONTCAR_FILES[i].find('100') > 0:
        Surface = '100'
    Substrate = ASE_CONTCAR[0].symbol
    Substrate_Mass = ASE_CONTCAR[0].mass
    Adsorbate = ''
    Adsorbate_Mass = 0
    for j in range(0,numAdsorbates):
        Adsorbate = (Adsorbate+ASE_CONTCAR[-1*numAdsorbates+j].symbol)
        Adsorbate_Mass = (Adsorbate_Mass+ASE_CONTCAR[-1*numAdsorbates+j].mass)
    ASE_CONTCAR.set_constraint()
    ASE_CONTCAR2 = ASE_CONTCAR.repeat((repeatX,repeatY,1))  
    Coordinates = ASE_CONTCAR2.positions
    Oatom = ASE_CONTCAR2[repeatX*repeatY*numAtoms-numAdsorbates].position
    mindistance = 1000
    NN2min = 1000
    for j in range(0,repeatX*repeatY*numAtoms-numAdsorbates):
        distance = sum((Oatom-Coordinates[j])**2)**(0.5)
        if distance <= mindistance and np.abs(Oatom[2]-Coordinates[j][2]) > 0.05:
            mindistance = distance
            horizontal = ((Oatom[0]-Coordinates[j][0])**2+(Oatom[1]-Coordinates[j][1])**2)**(0.5)
            vertical = abs(Oatom[2]-Coordinates[j][2])
    for j in range(0,repeatX*repeatY*numAtoms-numAdsorbates):
        distance = sum((Oatom-Coordinates[j])**2)**(0.5)
        NN2vertTest = abs(Oatom[2]-Coordinates[j][2])
        NN2horTest = ((Oatom[0]-Coordinates[j][0])**2+(Oatom[1]-Coordinates[j][1])**2)**(0.5)
        if distance < NN2min and abs(vertical-NN2vertTest) >1:
            NN2min = distance
            NN2vertical = abs(Oatom[2]-Coordinates[j][2])
            NN2horizontal = ((Oatom[0]-Coordinates[j][0])**2+(Oatom[1]-Coordinates[j][1])**2)**(0.5)

    cosa = vertical/mindistance
    sina = horizontal/mindistance
    NN2cosa = NN2vertical/NN2min
    NN2sina = NN2horizontal/NN2min
    if Surface == '100' and CONTCAR_FILES[i].find('Bridge') > 0:
        NN = 2
        NN2 = 2
    elif Surface == '100':
        NN = 4
        NN2=1
    elif Surface == '110' and Substrate in ['Cr','Fe','Mo','V','W']:
        NN = 3
        NN2 = 1
    elif Surface == '110':
        NN = 2
        NN2 = 2
    elif Surface =='111':
        NN = 3
        NN2 = 3
    Mr = (1/Adsorbate_Mass+1/(NN*Substrate_Mass*cosa**2+(mindistance/NN2min)**2*NN2*Substrate_Mass*NN2cosa**2))**(-1)
    Frequencies = Text_Parse.file2listflex(freq_FILES[i],[58,61],'cm-1',46,57)
    Frequencies = Text_Parse.string2float(Frequencies)
    Eigenvector = Text_Parse.file2listfixed(freq_FILES[i],0,10**6,0,None)
    Eigenvector = Text_Parse.list_split(Eigenvector,' ')
    Eigenvector = Text_Parse.string2float(Eigenvector)
    Eigenvector = Eigenvector[0:(2+numAtoms)*3*numAdsorbates]
    newEigen = []
    for j in range(0,len(Eigenvector)):
        if len(Eigenvector[j]) == 6:
            newEigen.append(Eigenvector[j])
    Eigenvector=newEigen
    newEigen = []
    for j in range(0,len(Eigenvector)):
        if abs(Eigenvector[j][3]) + abs(Eigenvector[j][4]) +abs(Eigenvector[j][5]) != 0:
            newEigen.append(Eigenvector[j])
    Eigenvector=newEigen
    maxval = 0
    Zfrequency = 0
    if len(Frequencies) == numAdsorbates*3:
        for j in range(0,numAdsorbates*3):
            if abs(Eigenvector[j][5]) - abs(Eigenvector[j][4]) - abs(Eigenvector[j][3]) > maxval:
                Zfrequency = Frequencies[j]
                maxval = abs(Eigenvector[j][5]) - abs(Eigenvector[j][4]) - abs(Eigenvector[j][3])

    Energy = Text_Parse.file2listflex(OSZICAR_FILES[i],23,['E'],27,42)
    Energy = Text_Parse.list_split(Energy," ")
    Energy = Text_Parse.string2float(Energy)
    Energy = np.array(Energy)
    if len(Energy) > 0:
        Energy=Energy[-1]
        Energy = Energy[0]
    Metal_Surf = Substrate + Surface
    Surface_Energy = [j for j in Free_Surfaces if (Metal_Surf in j)]
    Surface_Energy = Text_Parse.file2listflex(Surface_Energy[0],23,['E'],27,42)
    Surface_Energy = Text_Parse.list_split(Surface_Energy," ")
    Surface_Energy = Text_Parse.string2float(Surface_Energy)
    Surface_Energy = np.array(Surface_Energy)
    if len(Surface_Energy) > 0:
        Surface_Energy=Surface_Energy[-1]
        Surface_Energy = Surface_Energy[0]
    frequency_corrected = Zfrequency*(Adsorbate_Mass/Mr)**0.5
    data = [Surface, Substrate, Substrate_Mass, NN, NN2, Adsorbate, Adsorbate_Mass,cosa,sina, Mr, Zfrequency,frequency_corrected,Energy, Surface_Energy,Energy-Surface_Energy]    
    Data.append(data)
        
myfile = open('Energy_Frequency_Data.csv', 'wb')
wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
for i in range(0,len(Data)):
    wr.writerow(Data[i])
myfile.close()