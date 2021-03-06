# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:58:40 2016

@author: Josh
"""
from __future__ import division
import os
import numpy as np
import Text_Parse
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mat
from ase.io import read

Folder = 'Pt111'
mAdsorbate = 16
mMetal = 195
NN=3

mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 10
m1 = mAdsorbate*1.6737236*10**(-27) #kg
m2 = 1000000*1.6737236*10**(-27) #kg
mu = m1*m2/(m1+m2) #kg
directory=os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research//VASP_Files/O_Hollow_noD/%s' % Folder)
full_file_paths = Text_Parse.get_filepaths(directory)
CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
CONTCAR_FILES = [i for i in CONTCAR_FILES if not ('freq' in i)]
OSZICAR_FILES = Text_Parse.list_contains(full_file_paths,"OSZICAR")
File_Names = Text_Parse.between_values(CONTCAR_FILES,"%s\\" % Folder,"\\CONTCAR")

#labels1212 = np.genfromtxt('./Lansford_Stats_HW3_data/12.12.csv', delimiter=',', dtype=str)[0,:]
#data1212 = np.genfromtxt('C:\Users\Josh\Desktop\XSD files', delimiter=',')[1:,:]
Distances = []
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
    distance = abs((Coordinates[64,2]-Coordinates[63,2]))
    distance = distance.flatten()    
    Distances = np.hstack((Distances,distance))
Energies = []
for i in range(0,num_Files):
    Energy = Text_Parse.file2listflex(OSZICAR_FILES[i],23,['E'],27,42)
    Energy = Text_Parse.list_split(Energy," ")
    Energy = Text_Parse.string2float(Energy)
    Energy = np.array(Energy)
    Energy=Energy[-1]
    Energies = np.hstack((Energies,Energy))
    
#re = Distances[minIndex]
#Distances = Distances - re
D0 = min(Energies)
minIndex = (np.ndarray.tolist(Energies)).index(D0)
Energies = Energies - D0

def func(x,De,a,re):
    if a >0 and re > 0:
        return De*(1-np.exp(-1*a*(x-re)))**2
    else:
        return 1e6

c = 2.9979245800*10**(8) #m/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules

#Full fit
plt.figure(figsize=(7,5),dpi=500)

#idx = (Distances>.9)*(Distances<1.36)
#Distances = Distances[idx]
#Energies=Energies[idx]

ppot, pcov = curve_fit(func,Distances,Energies,p0=np.array([3.2,4,Distances[minIndex]]),maxfev = 5000)
De = ppot[0]
a = ppot[1]
re=ppot[2]
Distance_sorted = np.linspace(min(Distances),max(Distances),50)
Energy_fitted = func(Distance_sorted,De,a,re)
plt.plot(Distances,Energies,'o',Distance_sorted,Energy_fitted)
plt.title('Pt Frequency D2',size=14, fontweight='bold')
plt.ylabel('Potential Energy (eV)',size=14, fontweight='bold')
plt.xlabel('Distance from Equilibrium (A)',size=14, fontweight='bold')
plt.xticks(size=14)
plt.yticks(size=14)
plt.show()


#Harmonic Approximation
p1 = np.polyfit(Distances, Energies, 2)

#spring constant (J/m^2)
k = p1[0]*2*JeV*10**(20)
#Reduced mass (kg)
vH = 1/(c*2*np.pi)*(k/mu)**0.5/100 #Harmonic Frequency
print(Folder)
print('De')
print(De)
print('a')
print(a)
print('Harmonic')
print(vH)
frequency = a*10**(8)/np.pi*(De*JeV/(2*mu))**(0.5)/c #cm^-1 Morse frequency
print('Morse')
print(frequency)
print('bond length')
print(re)
'Variance and Covariance'
print(pcov)



# Reduced Mass
#Repeating surface
a=2
b=2
Geometry = [['Surface', 'Distance', 'cos(alpha)','sin(alpha)', 'NN2','NN2/NN','NN2cos(alpha)','NN2sin(alpha)']]
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
            
    
    surface_distances = [File_Names[i], mindistance, vertical/mindistance,horizontal/mindistance,NN2,NN/mindistance,NN2vertical/NN2,NN2horizontal/NN2]    
    Geometry.append(surface_distances)
        
            
MrAdsorbate = (1/mAdsorbate+1/(NN*mMetal*(np.array(np.array(Geometry,ndmin=0)[1:,2],dtype='d'))**2))**(-0.5)
print('Reduced mass of system^.5')
print(MrAdsorbate)