# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 22:06:40 2016

@author: Josh
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import adjustText
import os
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 15
mat.rcParams['legend.numpoints'] = 1
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/O_H_N_Scaling_Darpa.csv')
Metal_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,0:2]
Data_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,2:]
Data = np.genfromtxt(vibration_file, delimiter=',')[1:,2:]
Metal_Info = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[1:,0:2]



#PLotting Eads vs Mr*frequency
'''Eads and frequency'''
EadsO = Data[:,0]
vMO = Data[:,1]
EadsH = Data[:,2]
vMH = Data[:,3]
EadsN = Data[:,4]
vMN = Data[:,5]
'''Reduced masses'''
MrO = 16
MrH = 1
MrN = 14

'''frequency times square root of reduced mass'''
vMOmR = vMO*MrO**0.5
vMNmR = vMN*MrN**0.5
vMHmR = vMH*MrH**0.5

'''Concatenating'''
vMmR = np.concatenate((vMOmR,vMHmR,vMNmR))
Eads = np.concatenate((EadsO,EadsH,EadsN))

'''Linear fit'''
p = np.polyfit(Eads, vMmR, 1) 
m1 = p[0]
b1=p[1]
plt.figure(1)
plt.figure(figsize=(16,8),dpi=500)


plt.plot(EadsN,vMNmR,'go')
plt.plot(EadsO,vMOmR,'bo')
plt.plot(EadsH,vMHmR,'ro')
plt.plot(np.array([min(Eads),max(Eads)]),b1+m1*np.array([min(Eads),max(Eads)]),'k')
plt.title('Atomic Adsorbates on Face Centered Cubic Metal Surfaces: 111 Facet, PBE-D3',size=20, fontweight='bold')
plt.xlabel('${\Delta}$E$_{ads}$ (eV)',size=20, fontweight='bold')
plt.ylabel('Frequency X $\sqrt{Reduced\:Mass}}$ (amu$^{0.5}$/cm)',size=20, fontweight='bold')
plt.legend(['v(M-N)','v(M-O)','v(M-H)'],loc=1,prop={'size':25})
plt.xticks(size=20)
plt.yticks(size=20)

mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,1]))
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,1]))
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,1]))
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(Eads, vMmR, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=30, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()