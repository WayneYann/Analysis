# -*- coding: utf-8 -*-
"""
Created on Sat October 15 2016

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

#Input parameters
Directory = 'NH_Hollow_Morse'
#Directory = 'Jellium/K_pvfixed'
Eads = 5.3189
Folder = 'Pt111'
#Folder = 'O111'
mAdsorbate = 16
mMetal = 27
NN=3
mEOtoOH = 0.43
mENtoNH = 0.68
''''''''''''''''''''''''''''
No need to change below here
'''''''''''''''''''''''''''''

mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 10
directory=os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/VASP_Files/%s/%s' % (Directory,Folder))
full_file_paths = Text_Parse.get_filepaths(directory)
CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
CONTCAR_FILES = [i for i in CONTCAR_FILES if not ('freq' in i)]
OSZICAR_FILES = Text_Parse.list_contains(full_file_paths,"OSZICAR")
OSZICAR_FILES = [i for i in OSZICAR_FILES if not ('freq' in i)]
File_Names = Text_Parse.between_values(CONTCAR_FILES,"%s" % Directory,"CONTCAR")

#labels1212 = np.genfromtxt('./Lansford_Stats_HW3_data/12.12.csv', delimiter=',', dtype=str)[0,:]
#data1212 = np.genfromtxt('C:\Users\Josh\Desktop\XSD files', delimiter=',')[1:,:]
Distances = []
num_Files = len(CONTCAR_FILES)
for i in range(0,num_Files):
    ASE_CONTCAR = read(CONTCAR_FILES[i])
    numAtoms = len(ASE_CONTCAR)
    numAdsorbates = numAtoms-64
    ASE_CONTCAR.set_constraint()
    Adsorbate_index = numAtoms-numAdsorbates
    vectDistance = ASE_CONTCAR.get_distances(Adsorbate_index,range(0,Adsorbate_index),mic=True,vector=True)
    distances = ASE_CONTCAR.get_distances(Adsorbate_index,range(0,Adsorbate_index),mic=True,vector=False)
    #distances = abs(vectDistance[:,2])
    mindistance = min(distances)
    Distances = np.hstack((Distances,mindistance))
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
mindistance = Distances[minIndex]
def func(x,De,a,re):
    if a > 0 and re > 0:
        return De*(1-np.exp(-1*a*(x-re)))**2
    else:
        return 1e6

def func2(x,a):
    if a >0:
        return Eads*(1-np.exp(-1*a*(x-mindistance)))**2
    else:
        return 1e6
c = 2.9979245800*10**(8) #m/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules
eVperwave = 1.23981*10**(-4)

#Full fit
idx = (Distances-mindistance>=.1); Distances = Distances[idx]; Energies=Energies[idx]

plt.figure(figsize=(7,5),dpi=500)
ppot, pcov = curve_fit(func,Distances,Energies,p0=np.array([3.2,0.8,mindistance]),maxfev = 5000)
De = ppot[0]
a = ppot[1]
re=ppot[2]
ppot2, pcov2 = curve_fit(func2,Distances,Energies,p0=np.array([0.8]),maxfev = 5000)
a2 = ppot2[0]
Distance_sorted = np.linspace(min(Distances),max(Distances),50)
Energy_fitted = func(Distance_sorted,De,a,re)
Energy_fitted2 = func2(Distance_sorted,a2)
plt.plot(Distances,Energies,'o',Distance_sorted,Energy_fitted,'g',Distance_sorted,Energy_fitted2,'r')
plt.title('(Pt-O) Energy on 111-Hollow Site: PBE-D3',size=14, fontweight='bold')
plt.ylabel('Potential Energy (eV)',size=14, fontweight='bold')
plt.xlabel(r'Adsorbate-NN distance ($\AA$)',size=14, fontweight='bold')
plt.legend(['Data','Fit De and re','Fixed De and re'],loc=2)
plt.xticks(size=14)
plt.yticks(size=14)
plt.show()

print(Folder)
print('De')
print(De)
print('a')
print(a)
print('a2')
print(a2)
print('bond length')
print(re)
'Variance and Covariance'
#print(pcov)     
