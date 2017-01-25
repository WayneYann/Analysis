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
import statmech



directory=os.path.expanduser('~/Box Sync/Synced_Files\Coding\Research\VASP_Files\\551K500eV_compressed')
full_file_paths = Text_Parse.get_filepaths(directory)

CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
OSZICAR_FILES = Text_Parse.list_contains(full_file_paths,"OSZICAR")
OUTCAR_FILES = Text_Parse.list_contains(full_file_paths,"OUTCAR")
freq_FILES = Text_Parse.list_contains(full_file_paths,"freq")
freq_FILES = [i for i in freq_FILES if not any(word in i for word in ['Free_Surfaces','Gases'])]
CONTCAR_FILES = [i for i in CONTCAR_FILES if not any(word in i for word in ['Free_Surfaces','Gases'])]
OSZICAR_FILES = [i for i in OSZICAR_FILES if not any(word in i for word in ['Free_Surfaces','Gases'])]
OSZICAR_FILES = [i for i in OSZICAR_FILES if not any(word in i for word in ['Free_Surfaces','Gases'])]
OUTCAR_FILES = [i for i in OUTCAR_FILES if not any(word in i for word in ['Free_Surfaces','Gases'])]
Surface_Names = Text_Parse.between_values(CONTCAR_FILES,"\\551K500eV_compressed\\","_CONTCAR")
Free_Surfaces = Text_Parse.list_contains(full_file_paths,"Free_Surfaces")
Free_OUTCAR = [i for i in Free_Surfaces if ('OUTCAR' in i)]
Free_Surfaces = [i for i in Free_Surfaces if ('OSZICAR' in i)]
Data = [['Surface', 'Substrate','Adsorbate','Substrate Mass', 'NN','mindistance','NN2','NN2min','Adsorbate mass','cos(alpha)','sin(alpha)', 'NN2cosa', 'NN2sina', 'Mr','Energy','Surface_Energy','Egas','Eads','Edisp','SurfaceDisp','EminusDisp','Gvib100C','ZPE','vibEntropy','vibEnthalpy']]
num_Files = len(Surface_Names)
c = 29979245800 #cm/s
JeV = 1.60217653*10**(-19) #eV to Joules
NA = 6.0221409*10**(23)
freqoffset = 0
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
    if Adsorbate == 'C':
        Egas = -1.3953603
    elif Adsorbate == 'CH':
        Egas = -6.2283107
    elif Adsorbate == 'CHH':
        Egas = -12.137645
    elif Adsorbate == 'CHHH':
        Egas = -18.230307
    elif Adsorbate == 'O':
        Egas = -1.9555224
    elif Adsorbate == 'OH':
        Egas = -7.7312472
    elif Adsorbate == 'OHH':
        Egas = -14.160252
    elif Adsorbate == 'N':
        Egas = -3.1718078
    elif Adsorbate == 'NH':
        Egas = -8.1413821
    elif Adsorbate == 'NHH':
        Egas = -13.541416
    elif Adsorbate == 'H':
        Egas =-1.2065081
    elif Adsorbate =='CO':
        Egas = -14.454577
    elif Adsorbate =='OO':
        Egas = -9.5884297
    elif Adsorbate =='OOH': 
        Egas = -10.625912
    elif Adsorbate == 'NN':
        Egas = -15.476786
    elif Adsorbate =='CCHHHHH':
        Egas = -34.751018
    elif Adsorbate == 'CCHHHH':
        Egas = -31.743830
    ASE_CONTCAR.set_constraint()
    Adsorbate_index = numAtoms-numAdsorbates
    vectDistance = ASE_CONTCAR.get_distances(Adsorbate_index,range(0,Adsorbate_index),mic=True,vector=True)
    distances = ASE_CONTCAR.get_distances(Adsorbate_index,range(0,Adsorbate_index),mic=True,vector=False)
    Oatom = ASE_CONTCAR[numAtoms-numAdsorbates].position
    mindistance = 1000
    NN2min = 1000
    for j in range(0,Adsorbate_index):
        if distances[j] <= mindistance and np.abs(vectDistance[j][2]) > 0.05:
            mindistance = distances[j]
            horizontal = ((vectDistance[j][0])**2+(vectDistance[j][1])**2)**(0.5)
            vertical = abs(vectDistance[j][2])
    for j in range(0,Adsorbate_index):
        NN2vertTest = abs(vectDistance[j][2])
        NN2horTest = ((vectDistance[j][0])**2+(vectDistance[1])**2)**(0.5)
        if distances[j] < NN2min and abs(vertical-NN2vertTest) >1:
            NN2min = distances[j]
            NN2vertical = abs(vectDistance[j][2])
            NN2horizontal = ((vectDistance[j][0])**2+(vectDistance[j][1])**2)**(0.5)

    cosa = vertical/mindistance
    sina = horizontal/mindistance
    NN2cosa = NN2vertical/NN2min
    NN2sina = NN2horizontal/NN2min
    if CONTCAR_FILES[i].find('atop') > 0:
        NN = 1
        NN2 = 1
    elif CONTCAR_FILES[i].find('Bridge') > 0 and CONTCAR_FILES[i].find('Hollow') < 0:
        NN = 2
        NN2 = 4
    elif CONTCAR_FILES[i].find('Bridge') > 0:
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
    Mr = (1/Adsorbate_Mass+1/(NN*Substrate_Mass*cosa**2))**(-1)
    if CONTCAR_FILES[i][:-8] ==  freq_FILES[i-freqoffset][:-5]:
        Frequencies = Text_Parse.file2listflex(freq_FILES[i-freqoffset],[58,61],'cm-1',46,57)
        RorI = Text_Parse.file2listflex(freq_FILES[i-freqoffset],[58,61],'cm-1',5,7)
        RorI = [x.strip(' ') for x in RorI]
        Frequencies = Text_Parse.string2float(Frequencies)
        vibGibbs = statmech.vibgibbs(Frequencies,373.15)
        ZPE = statmech.vibenergy(Frequencies,0)
        vibEntropy = statmech.vibentropy(Frequencies,373.15)
        vibEnthalpy = statmech.vibenergy(Frequencies,373.15)
        if 'f/i' in RorI:
            vibGibbs = ''
            ZPE = ''
            vibEntropy = ''
            vibEnthalpy = ''               
    else:
        freqoffset = freqoffset+1
        vibGibbs = ''
        vibEntropy = ''
        vibEnthalpy = ''
        ZPE = ''
    
    if vibGibbs == 0:
        vibGibbs = ''
        vibEntropy = ''
        vibEnthalpy = ''
        ZPE = ''
    numAdsorbates = numAtoms-64
    Energy = Text_Parse.file2listflex(OSZICAR_FILES[i],23,['E'],27,42)
    Energy = Text_Parse.list_split(Energy," ")
    Energy = Text_Parse.string2float(Energy)
    Energy = np.array(Energy)
    if len(Energy) > 0:
        Energy=Energy[-1]
        Energy = Energy[0]
    else:
        Energy = float('NaN')
        
    Edisp = Text_Parse.file2listflex(OUTCAR_FILES[i],[1,10],['Edisp (eV)'],11,23)
    if len(Edisp) > 0:
        Edisp = Text_Parse.string2float(Edisp)
        Edisp = Edisp[0]
    else:
        Edisp = float('NaN')
    
    Metal_Surf = Substrate + Surface
    Surface_Energy = [j for j in Free_Surfaces if (Metal_Surf in j)]
    Surface_Energy = Text_Parse.file2listflex(Surface_Energy[0],23,['E'],27,42)
    Surface_Energy = Text_Parse.list_split(Surface_Energy," ")
    Surface_Energy = Text_Parse.string2float(Surface_Energy)
    Surface_Energy = np.array(Surface_Energy)
    if len(Surface_Energy) > 0:
        Surface_Energy=Surface_Energy[-1]
        Surface_Energy = Surface_Energy[0]
    De = -1*(Energy-Surface_Energy-Egas)
    
    Surface_Disp = [j for j in Free_OUTCAR if (Metal_Surf in j)]
    Surface_Disp = Text_Parse.file2listflex(Surface_Disp[0],[1,10],['Edisp (eV)'],11,23)
    Surface_Disp = Text_Parse.string2float(Surface_Disp)
    SurfaceDisp = Surface_Disp[0]    
    
    """substracts ratio of 12-6 potential to estimate dispersion contribution with molecule below it."""
    Mr = (1/Adsorbate_Mass+1/(NN*Substrate_Mass*cosa**2+((mindistance/NN2min)**6-(mindistance/NN2min)**12)*NN2*Substrate_Mass*NN2cosa**2))**(-1)
    data = [Surface, Substrate, Adsorbate, Substrate_Mass, NN, mindistance,NN2, NN2min, Adsorbate_Mass,cosa,sina, NN2cosa, NN2sina, Mr, Energy, Surface_Energy,Egas,Energy-Surface_Energy-Egas,Edisp,SurfaceDisp,float(Energy-Edisp-Surface_Energy+SurfaceDisp-Egas),vibGibbs,ZPE,vibEntropy,vibEnthalpy]    
    Data.append(data)

    
myfile = open('OOH_Data.csv', 'wb')
wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
for j in range(0,len(Data)):
    wr.writerow(Data[j])
myfile.close()