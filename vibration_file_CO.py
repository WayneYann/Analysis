# -*- coding: utf-8 -*-
"""
Created on Sun Jan 08 12:04:03 2017

@author: Josh
"""
from __future__ import division
import numpy as np
import matplotlib as mat
import statmech
from ase.thermochemistry import IdealGasThermo
from ase import Atoms
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['legend.numpoints'] = 1
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 20
c = 29979245800 #cm/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules
COr = (0.51388868234-0.486111317659)*41.3427602451982
COfreq = 0.260119879
CO = Atoms('CO')
CO[1].z = COr
COThermo = IdealGasThermo([COfreq],'linear',atoms=CO,symmetrynumber=1,spin=0,natoms=2)
COh100 = COThermo.get_enthalpy(100)
COh100 = COThermo.get_enthalpy(145)
COh100 = COThermo.get_enthalpy(130)
COh100 = COThermo.get_enthalpy(420)
COh100 = COThermo.get_enthalpy(350)
COh100 = COThermo.get_enthalpy(450)
COh100 = COThermo.get_enthalpy(340)
COh100 = COThermo.get_enthalpy(500)

H_Ag=np.array([234.060289,22.823011,16.915597,16.786227,3.576371,2.197334])*8.06573
vibAg = statmech.vibenergy(H_Ag,100)
print(vibAg)
H_Cu=np.array([223.437829,31.524364,27.338489,27.112873,13.374554,12.939494])*8.06573
vibCu = statmech.vibenergy(H_Cu,130)
print(vibCu)
H_Ir=np.array([211.367202,44.397834,31.685160,30.837878,18.670770,17.712435])*8.06573
vibIr = statmech.vibenergy(H_Ir,420)
print(vibIr)
H_Ni=np.array([217.116607,43.218982,34.278379,33.953895,17.172639,16.697658])*8.06573
vibNi = statmech.vibenergy(H_Ni,350)
print(vibNi)
H_Pd=np.array([216.982268,43.592890,43.227109,40.159294,20.334542,20.123495])*8.06573
vibPd = statmech.vibenergy(H_Pd,450)
print(vibPd)
H_Pt=np.array([214.733604,43.390842,39.941662,39.565127,20.018194,19.838733])*8.06573
vibPt = statmech.vibenergy(H_Pt,340)
print(vibPt)
H_Rh=np.array([214.533897,42.275580,33.186760,33.035735,18.156254,18.030850])*8.06573
vibRh = statmech.vibenergy(H_Rh,500)
print(vibRh)

A_Ag=np.array([252.208474,21.286667,17.518939,17.512364,2.346201,2.023377])*8.06573
vibAg = statmech.vibenergy(A_Ag,100)
print(vibAg)
A_Au=np.array([253.789165,33.837780,22.946976,22.905899,3.650928,2.801318])*8.06573
vibAu = statmech.vibenergy(A_Au,145)
print(vibAu)  
A_Cu=np.array([247.546678,38.896534,34.556057,34.315866,6.603755,5.506249])*8.06573
vibCu = statmech.vibenergy(A_Cu,130)
print(vibCu)
A_Ir=np.array([247.941017,62.155735,56.711446,56.487965,8.368375,7.347044])*8.06573
vibIr = statmech.vibenergy(A_Ir,420)
print(vibIr)
A_Ni=np.array([246.348131,52.368468,43.062344,43.012748,5.886140,3.446441])*8.06573
vibNi = statmech.vibenergy(A_Ni,350)
print(vibNi)
A_Pd=np.array([249.947125,50.552814,37.486792,37.196734,4.173107,2.577255])*8.06573
vibPd = statmech.vibenergy(A_Pd,450)
print(vibPd)
A_Pt=np.array([253.646662,60.692124,49.159540,48.774832,6.950418,6.466714])*8.06573
vibPt = statmech.vibenergy(A_Pt,340)
print(vibPt)
A_Rh=np.array([245.273834,57.170317,50.020844,49.807444,7.081101,5.850451])*8.06573
vibRh = statmech.vibenergy(A_Rh,500)
print(vibRh)


H_Ag=np.array([234.060289,22.823011,16.915597,16.786227,3.576371,2.197334])*8.06573
vibAg = statmech.vibenergy(H_Ag,298)
print(vibAg)
H_Cu=np.array([223.437829,31.524364,27.338489,27.112873,13.374554,12.939494])*8.06573
vibCu = statmech.vibenergy(H_Cu,298)
print(vibCu)
H_Ir=np.array([211.367202,44.397834,31.685160,30.837878,18.670770,17.712435])*8.06573
vibIr = statmech.vibenergy(H_Ir,298)
print(vibIr)
H_Ni=np.array([217.116607,43.218982,34.278379,33.953895,17.172639,16.697658])*8.06573
vibNi = statmech.vibenergy(H_Ni,298)
print(vibNi)
H_Pd=np.array([216.982268,43.592890,43.227109,40.159294,20.334542,20.123495])*8.06573
vibPd = statmech.vibenergy(H_Pd,298)
print(vibPd)
H_Pt=np.array([214.733604,43.390842,39.941662,39.565127,20.018194,19.838733])*8.06573
vibPt = statmech.vibenergy(H_Pt,298)
print(vibPt)
H_Rh=np.array([214.533897,42.275580,33.186760,33.035735,18.156254,18.030850])*8.06573
vibRh = statmech.vibenergy(H_Rh,298)
print(vibRh)

A_Ag=np.array([250.91862,22.8854,16.33,16.044,1.099,1.099])*8.06573
vibAg = statmech.vibenergy(A_Ag,298)
print(vibAg)
A_Au=np.array([253.789165,33.837780,22.946976,22.905899,3.650928,2.801318])*8.06573
vibAu = statmech.vibenergy(A_Au,298)
print(vibAu)  
A_Cu=np.array([247.546678,38.896534,34.556057,34.315866,6.603755,5.506249])*8.06573
vibCu = statmech.vibenergy(A_Cu,298)
print(vibCu)
A_Ir=np.array([247.941017,62.155735,56.711446,56.487965,8.368375,7.347044])*8.06573
vibIr = statmech.vibenergy(A_Ir,298)
print(vibIr)
A_Ni=np.array([246.348131,52.368468,43.062344,43.012748,5.886140,3.446441])*8.06573
vibNi = statmech.vibenergy(A_Ni,298)
print(vibNi)
A_Pd=np.array([249.947125,50.552814,37.486792,37.196734,4.173107,2.577255])*8.06573
vibPd = statmech.vibenergy(A_Pd,298)
print(vibPd)
A_Pt=np.array([253.646662,60.692124,49.159540,48.774832,6.950418,6.466714])*8.06573
vibPt = statmech.vibenergy(A_Pt,298)
print(vibPt)
A_Rh=np.array([245.273834,57.170317,50.020844,49.807444,7.081101,5.850451])*8.06573
vibRh = statmech.vibenergy(A_Rh,298)
print(vibRh)

