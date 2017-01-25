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



directory=os.path.expanduser('~/Box Sync/Synced_Files\Coding\Research\VASP_Files\VASP_Ad_v02')
full_file_paths = Text_Parse.get_filepaths(directory)

CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
OSZICAR_FILES = Text_Parse.list_contains(full_file_paths,"OSZICAR")
OUTCAR_FILES = Text_Parse.list_contains(full_file_paths,"OUTCAR")
freq_FILES = Text_Parse.list_contains(full_file_paths,"freq")
CONTCAR_FILES = [i for i in CONTCAR_FILES if not ('Free_Surfaces' in i)]
OSZICAR_FILES = [i for i in OSZICAR_FILES if not ('Free_Surfaces' in i)]
OSZICAR_FILES = [i for i in OSZICAR_FILES if not ('Free_Surfaces' in i)]
OUTCAR_FILES = [i for i in OUTCAR_FILES if not ('Free_Surfaces' in i)]
Surface_Names = Text_Parse.between_values(CONTCAR_FILES,"\\VASP_Ad_v02\\","_CONTCAR")
Free_Surfaces = Text_Parse.list_contains(full_file_paths,"Free_Surfaces")
Free_OUTCAR = [i for i in Free_Surfaces if ('OUTCAR' in i)]
Free_Surfaces = [i for i in Free_Surfaces if ('OSZICAR' in i)]
Data = [['Surface', 'Substrate','Adsorbate','Substrate Mass', 'NN','mindistance','NN2','NN2min','Adsorbate mass','cos(alpha)','sin(alpha)', 'NN2cosa', 'NN2sina', 'Mr','Zfrequency','frequency_corrected','Energy','Surface_Energy','Egas','Eads','Edisp','SurfaceDisp','EminusDisp','vibGibbs','ZPE','vibEntropy','vibEnthalpy','Gvibxyz','GvibxyzCF','Xfrequency','Yfrequency','HorC','Rotation','High_Frequency','RorIx','RorIy']]
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
        vibGibbs = statmech.vibgibbs(Frequencies,298)
        ZPE = statmech.vibenergy(Frequencies,0)
        vibEntropy = statmech.vibentropy(Frequencies,298)
        vibEnthalpy = statmech.vibenergy(Frequencies,298)
        Eigenvector = Text_Parse.file2listfixed(freq_FILES[i-freqoffset],0,10**6,0,None)
        Eigenvector = Text_Parse.list_split(Eigenvector,' ')
        Eigenvector = Text_Parse.string2float(Eigenvector)
        newEigen = []
        if 'f/i' in RorI:
            vibGibbs = ''
            ZPE = ''
            vibEntropy = ''
            vibEnthalpy = ''
        
        for j in range(0,len(Eigenvector)):
            if len(Eigenvector[j]) == 6:
                newEigen.append(Eigenvector[j])
        Eigenvector=newEigen
        newEigen = []
        for j in range(0,len(Eigenvector)):
            if abs(Eigenvector[j][3]) + abs(Eigenvector[j][4]) +abs(Eigenvector[j][5]) <> 0:
                newEigen.append(Eigenvector[j])
        Eigenvector=newEigen
        maxZval = 0
        maxZ2val = 0
        maxXval = 0
        maxYval = 0
        Zfrequency = 0
        Xfrequency = 0
        Yfrequency = 0
        Rotation = 0
        High_Frequency = 0
        RorIx = 0
        RorIy = 0
        if len(Frequencies) > 0:
            High_Frequency = max(Frequencies)
        else:
            High_Frequency = float('NaN')
        if len(Frequencies) == numAdsorbates*3:
            for j in range(0,numAdsorbates*3):
                val = 0
                val2 = 0
                Zsign=0
                Ysign = 0
                Xsign = 0
                #New CODE
                for k in range(0,numAdsorbates):
                    val = np.array(Eigenvector[(numAdsorbates*j+k)]) + val
                    val2 = abs(np.array(Eigenvector[(numAdsorbates*j+k)])) + val2
                    #val = val2
                    #checking that z movment is in same direction
                    Zsign = np.sign(Eigenvector[(numAdsorbates*j+k)][5]) + Zsign
                    Ysign = np.sign(Eigenvector[(numAdsorbates*j+k)][4]) + Ysign
                    Xsign = np.sign(Eigenvector[(numAdsorbates*j+k)][3]) + Xsign
                if (abs(val[5]) - abs(val[4]) - abs(val[3])) > maxZval and abs(Zsign) == numAdsorbates and abs(Eigenvector[j*numAdsorbates][5]) == max(abs(np.array([item[5] for item in Eigenvector[j*numAdsorbates:((j+1)*numAdsorbates)]]))): 
                    Zfrequency = Frequencies[j]
                    maxZval = abs(val[5]) - abs(val[4]) - abs(val[3])
                if (abs(val[4]) - abs(val[5]) - abs(val[3])) > maxYval and abs(Ysign) == numAdsorbates and abs(Eigenvector[j*numAdsorbates][4]) == max(abs(np.array([item[4] for item in Eigenvector[j*numAdsorbates:((j+1)*numAdsorbates)]]))):
                    Yfrequency = Frequencies[j]
                    RorIy = RorI[j]
                    maxYval = abs(val[4]) - abs(val[5]) - abs(val[3])
                if (abs(val[3]) - abs(val[5]) - abs(val[4])) > maxXval and abs(Xsign) == numAdsorbates and abs(Eigenvector[j*numAdsorbates][3]) == max(abs(np.array([item[3] for item in Eigenvector[j*numAdsorbates:((j+1)*numAdsorbates)]]))):
                    Xfrequency = Frequencies[j]
                    RorIx = RorI[j]
                    maxXval = abs(val[3]) - abs(val[5]) - abs(val[4])
                if (val2[5] - val2[4] - val2[3]) > maxZ2val and abs(Zsign) < numAdsorbates: 
                    Rotation = Frequencies[j]
                    maxZ2val = abs(val2[5]) - abs(val2[4]) - abs(val2[3])
            for j in range(0,numAdsorbates*3):
                val = 0
                Zsign=0
                Ysign = 0
                Xsign = 0
                for k in range(0,numAdsorbates):
                    val = np.array(Eigenvector[(numAdsorbates*j+k)]) + val
                    #val = abs(np.array(Eigenvector[(numAdsorbates*j+k)])) + val
                    #checking that movment is in same direction
                    Zsign = np.sign(Eigenvector[(numAdsorbates*j+k)][5]) + Zsign
                    Ysign = np.sign(Eigenvector[(numAdsorbates*j+k)][4]) + Ysign
                    Xsign = np.sign(Eigenvector[(numAdsorbates*j+k)][3]) + Xsign
                if (abs(val[5]) - abs(val[4]) - abs(val[3])) > maxZval and abs(Zsign) == numAdsorbates and Zfrequency == 0: 
                    Zfrequency = Frequencies[j]
                    maxZval = abs(val[5]) - abs(val[4]) - abs(val[3])
                if (abs(val[4]) - abs(val[5]) - abs(val[3])) > maxYval and abs(Ysign) == numAdsorbates and Yfrequency == 0:
                    Yfrequency = Frequencies[j]
                    RorIy = RorI[j]
                    maxYval = abs(val[4]) - abs(val[5]) - abs(val[3])
                if (abs(val[3]) - abs(val[5]) - abs(val[4])) > maxXval and abs(Xsign) == numAdsorbates and Xfrequency == 0:
                    Xfrequency = Frequencies[j]
                    RorIx = RorI[j]
                    maxXval = abs(val[3]) - abs(val[5]) - abs(val[4])

            if Adsorbate =='OO' and NN ==2:
                maxZval=0
                Zfrequency ==0
                for j in range(0,numAdsorbates*3):
                    val = 0
                    Zsign=0
                    for k in range(0,numAdsorbates):
                        val = np.array(Eigenvector[(numAdsorbates*j+k)]) + val
                        #checking that z movment is in same direction
                        Zsign = np.sign(Eigenvector[(numAdsorbates*j+k)][5]) + Zsign
                    if abs(val[5]) > maxZval and abs(Zsign) == numAdsorbates: 
                        Zfrequency = Frequencies[j]
                        maxZval = abs(val[5])
                for j in range(0,numAdsorbates*3):
                    val = 0
                    Zsign=0
                    for k in range(0,numAdsorbates):
                        val = np.array(Eigenvector[(numAdsorbates*j+k)]) + val
                        #checking that movment is in same direction
                        Zsign = np.sign(Eigenvector[(numAdsorbates*j+k)][5]) + Zsign
                    if (abs(val[5]) - abs(val[4]) - abs(val[3])) > maxZval and abs(Zsign) == numAdsorbates and Zfrequency == 0: 
                        Zfrequency = Frequencies[j]
                        maxZval = abs(val[5]) - abs(val[4]) - abs(val[3])

            if Adsorbate =='CCHHHHH':
                Zfrequency = 0
                Yfrequency = 0
                Xfrequency= 0
                maxZval=0
                for j in range(0,numAdsorbates*3):
                    val = 0
                    for k in range(0,numAdsorbates):
                        val = np.array(Eigenvector[(numAdsorbates*j+k)]) + val
                        #val = abs(np.array(Eigenvector[(numAdsorbates*j+k)])) + val
                    if abs(Eigenvector[j*numAdsorbates][5]) > maxZval:
                        Zfrequency = Frequencies[j]
                        maxZval = abs(Eigenvector[j*numAdsorbates][5])
                    #if (abs(val[5]) - abs(val[4]) - abs(val[3])) > maxZval and Zfrequency == 0 and abs(Eigenvector[j*numAdsorbates][5]) == max(abs(np.array([item[5] for item in Eigenvector[j*numAdsorbates:((j+1)*numAdsorbates)]]))): # 
                        #Zfrequency = Frequencies[j]
                        #maxZval = abs(val[5]) - abs(val[4]) - abs(val[3])
                    if (abs(val[4]) - abs(val[5]) - abs(val[3])) > maxYval and Yfrequency == 0 and abs(Eigenvector[j*numAdsorbates][4]) == max(abs(np.array([item[4] for item in Eigenvector[j*numAdsorbates:((j+1)*numAdsorbates)]]))):
                        Yfrequency = Frequencies[j]
                        RorIy = RorI[j]
                        maxYval = abs(val[4]) - abs(val[5]) - abs(val[3])
                    if (abs(val[3]) - abs(val[5]) - abs(val[4])) > maxXval and Xfrequency == 0 and abs(Eigenvector[j*numAdsorbates][3]) == max(abs(np.array([item[3] for item in Eigenvector[j*numAdsorbates:((j+1)*numAdsorbates)]]))):
                        Xfrequency = Frequencies[j]
                        RorIx = RorI[j]
                        maxXval = abs(val[3]) - abs(val[5]) - abs(val[4])
                
            if Adsorbate =='CCHHHHH':
                maxZval=0               
                for k in range(0,numAdsorbates):
                    val = np.array(Eigenvector[(numAdsorbates*j+k)]) + val
                    val = abs(np.array(Eigenvector[(numAdsorbates*j+k)])) + val
                if (abs(val[5]) - abs(val[4]) - abs(val[3])) > maxZval and Zfrequency == 0: # 
                    Zfrequency = Frequencies[j]
                    maxZval = abs(val[5]) - abs(val[4]) - abs(val[3])
                if (abs(val[4]) - abs(val[5]) - abs(val[3])) > maxYval and Yfrequency == 0:
                    Yfrequency = Frequencies[j]
                    RorIy = RorI[j]
                    maxYval = abs(val[4]) - abs(val[5]) - abs(val[3])
                if (abs(val[3]) - abs(val[5]) - abs(val[4])) > maxXval and Xfrequency == 0:
                    Xfrequency = Frequencies[j]
                    RorIx = RorI[j]
                    maxXval = abs(val[3]) - abs(val[5]) - abs(val[4])
                
    else:
        freqoffset = freqoffset+1
        vibGibbs = ''
        vibEntropy = ''
        vibEnthalpy = ''
        Zfrequency = 0
        Xfrequency = ''
        Yfrequency = ''
        
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
    a = (2*np.pi)*Zfrequency/c*((10**(-3)*Adsorbate_Mass/NA)/(2*De*JeV))**0.5*10**(10)
    
    Surface_Disp = [j for j in Free_OUTCAR if (Metal_Surf in j)]
    Surface_Disp = Text_Parse.file2listflex(Surface_Disp[0],[1,10],['Edisp (eV)'],11,23)
    Surface_Disp = Text_Parse.string2float(Surface_Disp)
    SurfaceDisp = Surface_Disp[0]    
    
    """substracts ratio of 12-6 potential to estimate dispersion contribution with molecule below it."""
    Hor_CF = Zfrequency*((1+NN*(Substrate_Mass/Adsorbate_Mass)*cosa**2)/(1+NN/2*(Substrate_Mass/Adsorbate_Mass)*sina**2))**0.5
    Mr = (1/Adsorbate_Mass+1/(NN*Substrate_Mass*cosa**2+((mindistance/NN2min)**6-(mindistance/NN2min)**12)*NN2*Substrate_Mass*NN2cosa**2))**(-1)
    frequency_corrected = Zfrequency*(Adsorbate_Mass/Mr)**0.5
    Hor_CF = frequency_corrected*(sina*(1/Adsorbate_Mass+1/(NN/2*Substrate_Mass*sina**2))/(2*cosa*(1/Adsorbate_Mass+1/(NN*Substrate_Mass*cosa**2))))**0.5
    if Zfrequency ==0:
        Zfrequency = ''
    if frequency_corrected ==0:
        frequency_corrected = ''       
    if Zfrequency <> '' and frequency_corrected <> '' and Xfrequency >0 and Yfrequency >0:
        Gvibxyz = statmech.vibgibbs([Zfrequency,Xfrequency,Yfrequency],298)
        GvibxyzCF = statmech.vibgibbs([frequency_corrected,Hor_CF,Hor_CF],298)
    else:
        Gvibxyz = ''
        GvibxyzCF = ''
    data = [Surface, Substrate, Adsorbate, Substrate_Mass, NN, mindistance,NN2, NN2min, Adsorbate_Mass,cosa,sina, NN2cosa, NN2sina, Mr, Zfrequency,frequency_corrected,Energy, Surface_Energy,Egas,Energy-Surface_Energy-Egas,Edisp,SurfaceDisp,float(Energy-Edisp-Surface_Energy+SurfaceDisp-Egas),vibGibbs,ZPE,vibEntropy,vibEnthalpy,Gvibxyz,GvibxyzCF,Xfrequency,Yfrequency,Hor_CF,Rotation,High_Frequency,RorIx,RorIy]    
    Data.append(data)

    
myfile = open('Energy_Frequency_Data_v02.csv', 'wb')
wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
for j in range(0,len(Data)):
    wr.writerow(Data[j])
myfile.close()