# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 23:35:06 2016

@author: Josh
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import adjustText
import os
mat.rcParams['mathtext.default'] = 'regular'
molecule_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Diatomic_Molecules_Data.csv')
Molecule_Label = np.genfromtxt(molecule_file, delimiter=',', dtype=str)[0,0]
Data_Labels = np.genfromtxt(molecule_file, delimiter=',', dtype=str)[0,1:]
Data = np.genfromtxt(molecule_file, delimiter=',')[1:,1:]
Molecule_Info = np.genfromtxt(molecule_file, delimiter=',', dtype=str)[1:,0]
DataLabels = []
for i in range(0,len(Data_Labels)):
    DataLabels.append(Data_Labels[i])
    DataLabels.append(i)
    
M1 = Data[:,3]
M2 = Data[:,4]
Mr = M1*M2/(M1+M2)
wMrsquared = Data[:,0]*Mr**(0.5)
De = Data[:,2]

#PLotting v(M-O) vs v(M-O2)
idx = np.isfinite(wMrsquared) & np.isfinite(De)
m,b = np.polyfit(De[idx], wMrsquared[idx], 1) 

mat.rc('text', usetex = True)
plt.figure(1)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(De,wMrsquared,'o',markersize=12)
plt.plot(De, m*De+b,'-b',lw=2)
plt.xlabel('Dissociation Energy (kJ/mol)',size=20, fontweight='bold')
plt.ylabel('$omega*Mr^{0.5}$ (g/mol)^{0.5}/cm',size=20, fontweight='bold')
plt.legend(['Frequency x (square of Mr): %sx+%s' %(round(m,2),round(b,2))],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)


mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(De, wMrsquared, Molecule_Info):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()