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

#Input parameters
Directory = 'Reduced_Mass_2'
Folder = 'CO111'
vCO = 493
mAdsorbate = 12+16
NN=1



''''''''''''''''''''''''''''
No need to change below here
'''''''''''''''''''''''''''''

mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 10
mH = 1.6737236*10**(-27) #kg
m1 = mAdsorbate*mH
m2 = 1000000*mH
mu = m1*m2/(m1+m2) #kg
directory=os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/VASP_Files/%s/%s' % (Directory,Folder))
full_file_paths = Text_Parse.get_filepaths(directory)
VASP_FILES = Text_Parse.list_contains(full_file_paths,"vasprun")
File_Names = Text_Parse.between_values(VASP_FILES,"%s" % Directory,"vasprun")

Distances = []
Energies = []
num_Files = len(VASP_FILES)
#for i in range(0,num_Files):

for i in range(0,num_Files):
    ASE_FILE = read(VASP_FILES[i])
    ASE_FILE.set_constraint()
    numatoms = len(ASE_FILE)
    numAdsorbates = numatoms-64
    Mass = 0
    MassDist=0
    for j in range(64,numatoms):
        distances = ASE_FILE.get_distances(j,range(0,64),mic=True)
        mindist = min(distances)
        mass = ASE_FILE.get_masses()[j]
        massdist = mindist*mass
        MassDist = massdist + MassDist
        Mass = Mass + mass
    Distance = MassDist/Mass
    #Distance = min(ASE_FILE.get_distances(64,range(0,64),mic=True))
    Distances.append(Distance)
    Energies.append(ASE_FILE.get_potential_energy())

#re = Distances[minIndex]
#Distances = Distances - re

c = 2.9979245800*10**(8) #m/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules
eVperwave = 1.23981*10**(-4)
mindistance = Distances[0]
def func(x,k,constant):
    if k > 0:
        
        return k*(x-mindistance)**2 + constant
    else:
        return 1e6
ppot, pcov = curve_fit(func,np.array(Distances),np.array(Energies),p0=np.array([26,-300]),maxfev = 5000)
#Harmonic Approximation
p1 = np.polyfit(Distances, Energies, 2)

#spring constant (J/m^2)
k = p1[0]*2*JeV*10**(20)
#Reduced mass (kg)
vH = 1/(c*2*np.pi)*(k/mu)**0.5/100 #Harmonic Frequency
print('Harmonic')
print(vH)

mu = k/(vCO*100*(c*2*np.pi))**2
print('Reduced Mass')
print(mu/mH)
