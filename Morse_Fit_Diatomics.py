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
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 10

directory=os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/VASP_Files/Diatomics/N2_harmonic')
full_file_paths = Text_Parse.get_filepaths(directory)
CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
OSZICAR_FILES = Text_Parse.list_contains(full_file_paths,"OSZICAR")
File_Names = Text_Parse.between_values(CONTCAR_FILES,"N2_harmonic\\","\\CONTCAR")

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
    distance = (sum((Coordinates[0,:]-Coordinates[1,:])**2))**0.5
    distance = distance.flatten()    
    Distances = np.hstack((Distances,distance))
Energies = []
for i in range(0,num_Files):
    Energy = Text_Parse.file2listflex(OSZICAR_FILES[i],23,['E'],27,42)
    Energy = Text_Parse.list_split(Energy," ")
    Energy = Text_Parse.string2float(Energy)
    Energy = np.array(Energy)
    Energy=Energy[0]
    Energies = np.hstack((Energies,Energy))
    
D0 = min(Energies)
minIndex = (np.ndarray.tolist(Energies)).index(D0)
Energies = Energies - D0
re = Distances[minIndex]
Distances = Distances - re

def func(x,De,a):
    return De*(1-np.exp(-1*a*x))**2
mH = 14*1.6737236*10**(-27) #kg
mu = mH*mH/(mH+mH) #kg
c = 2.9979245800*10**(8) #m/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules

#Full fit
plt.figure(figsize=(7,5),dpi=500)    
ppot, pcov = curve_fit(func,Distances,Energies)
De = ppot[0]
a = ppot[1]
Distance_sorted = np.linspace(min(Distances),max(Distances),50)
Energy_fitted = func(Distance_sorted,De,a)
plt.plot(Distances,Energies,'o',Distance_sorted,Energy_fitted)
plt.title('Fitted Morse Potential: Far from Minima',size=14, fontweight='bold')
plt.ylabel('Potential Energy (eV)',size=14, fontweight='bold')
plt.xlabel('Distance from Equilibrium (A)',size=14, fontweight='bold')
plt.xticks(size=14)
plt.yticks(size=14)
plt.show()
print(De)
print(a)

#Harmonic Approximation
p1 = np.polyfit(Distances, Energies, 2)

#spring constant (J/m^2)
k = p1[0]*2*JeV*10**(20)
#Reduced mass (kg)
vH = 1/(c*2*np.pi)*(k/mu)**0.5/100 #Harmonic Frequency
print(vH)
frequency = a*10**(8)/np.pi*(De*JeV/(2*mu))**(0.5)/c #cm^-1 Morse frequency
print(frequency)
print(pcov)