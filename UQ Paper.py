# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:05:01 2016

@author: lansford
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import adjustText
import os
import pandas as pd
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['text.latex.unicode'] = 'False'
mat.rcParams['legend.numpoints'] = 1
mat.rcParams['lines.linewidth'] = 4
mat.rcParams['lines.markersize'] = 16
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Data for UQ.csv')

Data = pd.read_csv(vibration_file,sep=',',header=0)

#PLotting v(M-O) vs v(M-OH)
p = np.polyfit(Data.H2,Data.O2, 1) 
m1 = p[0]
b1=p[1]
plt.figure(1)
plt.figure(figsize=(14,10))
plt.plot(np.array([np.min(Data.H2),np.max(Data.H2)]), m1*np.array([np.min(Data.H2),np.max(Data.H2)])+b1,'-b')
plt.plot(Data.H2, Data.O2,'ro')
plt.xlabel('${\Delta}$E$_{ads,H_{2}}$ (eV)',size=24)
plt.ylabel('${\Delta}$E$_{ads,O_{2}}$ (eV)',size=24)
plt.legend(['${\Delta}$E$_{ads,O_{2}}$=%s${\Delta}$E$_{ads,H_{2}}$ - %s' %(round(m1,2),abs(round(b1,2)))],loc=2,prop={'size':20},frameon=False)
plt.xticks(size=20)
plt.yticks(size=20)
plt.xlim([-2.1,0.6])
plt.ylim([-7.1,0])

mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(Data.Metal)):
    Marker.append(''.join(Data.Metal[i]))
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(Data.H2, Data.O2, Data.Metal):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, )
plt.show()
