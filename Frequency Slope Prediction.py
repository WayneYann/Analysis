# -*- coding: utf-8 -*-
"""
Created on Sat October 15 2016

@author: Josh
"""
from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import Morse_Func

#Input parameters
Directory1 = '\Jellium\Kfixed'
Folder1 = 'N111'
m1 = 14
Directory2 = '\Jellium\Kfixed'
Folder2 = 'NH111'
m2 = 15
Directory3 = '\Jellium\Kfixed'
Folder3 = 'NH2111'
Directory1 = 'PBE_Morse/O_Hollow_Morse'
Folder1 = 'Rh111'
m2 = 16
mE = 0.25
mindist = 0
maxdist = 1
#Directory1 = 'Jellium/Kfixed'
#Directory2 = 'Jellium/Kfixed'
#Folder1 = 'N111'
#Folder2 = 'NH2111'
''''''''''''''''''''''''''''
No need to change below here
'''''''''''''''''''''''''''''

mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 10
mat.rcParams['legend.numpoints'] = 1

directory=os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/VASP_Files/%s/%s' % (Directory1,Folder1))
ppot, Energies, Distances = Morse_Func.Morse_Fit(directory,mindist,maxdist)
De1 = ppot[0]
a1 = ppot[1]
re1=ppot[2]
print(a1)


directory2=os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/VASP_Files/%s/%s' % (Directory2,Folder2))
ppot2, Energies2, Distances2 = Morse_Func.Morse_Fit(directory2,mindist,maxdist)
De2 = ppot2[0]
a2 = ppot2[1]
re2=ppot2[2]

directory3=os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/VASP_Files/%s/%s' % (Directory3,Folder3))
ppot3, Energies3, Distances3 = Morse_Func.Morse_Fit(directory3,mindist,maxdist)
De3 = ppot3[0]
a3 = ppot3[1]
re3=ppot3[2]

x = np.linspace(-0.5,0.5,50)
Energy_fitted = Morse_Func.MorseEval(x,De1,a1)
Energy_fitted2 =  Morse_Func.MorseEval(x,De2,a2)
Energy_fitted3 =  Morse_Func.MorseEval(x,De3,a3)

#Full fit
plt.figure(figsize=(7,5),dpi=500)
plt.plot(Distances,Energies,'go',Distances2,Energies2,'bs',x,Energy_fitted,'g-',x,Energy_fitted2,'b-')
#plt.title('(Pt-O) Energy on 111-Hollow Site: RPBE-D3',size=14, fontweight='bold')
plt.ylabel(r'(${\Delta}$E$_{ads}$ - ${\Delta}$E$_{o}$) [eV]',size=14, fontweight='bold')
plt.xlabel(r'(r-r$_{o}$) [$\AA$]',size=14, fontweight='bold')
plt.legend(['O','OH'],loc=1)
plt.xticks(size=14)
plt.yticks(size=14)
plt.show()

plt.figure(figsize=(7,5),dpi=500)
plt.plot(Distances,Energies,'go',Distances2,Energies2,'bs',Distances3,Energies3,'r^',x,Energy_fitted,'g-',x,Energy_fitted2,'b-'
            ,x,Energy_fitted3,'r-')
#plt.title('(Pt-O) Energy on 111-Hollow Site: RPBE-D3',size=14, fontweight='bold')
plt.ylabel(r'(${\Delta}$E$_{ads}$ - ${\Delta}$E$_{o}$) [eV]',size=14, fontweight='bold')
plt.xlabel(r'(r-r$_{o}$) [$\AA$]',size=14, fontweight='bold')
plt.legend(['N','NH','NH2'],loc=1)
plt.xticks(size=14)
plt.yticks(size=14)
plt.show()

print(Folder1)
print('De')
print(De1)
print('a')
print(a1)
print('De2')
print(De2)
print('a2')
print(a2)
print('bond length')
print(re1)

"""predicted using bond length"""
mp = (re1/re2)**2*(m1/m2*re2**2/re1**2)**0.5
eD = np.linspace(-6,6,100)
mp = mE*(re1/re2)**2*(m1/m2*(a1**2*De1+eD/re1**2*6*7)/(a2**2*De2+mE*eD/re2**2*8*9))**0.5
mp = (re1/re2)**2*(m1/m2*a1**2*De1/(a2**2*De2))**0.5
#mpvect = mE*(re1/re2)**2*(m1/m2)**0.5*(2*a1**2*De1+eD*6*7/re1**2)**0.5/(2*a2**2*De2+mE*eD*6*7/re2**2)**0.5
#plt.plot(eD,mpvect)
"""predicted slope"""
mp2 = mE*(m1/m2*a2**2/a1**2*De1/De2)**0.5
"""predicted multiplier"""
#mp3 = (a2**2/a1**2*(De1-De1/2)/(De2-mE*De1/2))**0.5
mp3 = (m1/m2*a2**2/a1**2*(De1)/(De2))**0.5
print('multiplier')
print(mp3)
mp = re1**2/re2**2*(a1**2*De1/(a2**2*De2))**0.5
mp2 = re1/re2
print(mp)
print(mp2)