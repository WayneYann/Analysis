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
mat.rcParams['lines.linewidth'] = 8
mat.rcParams['lines.markersize'] = 24
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Data for UQ.csv')

Data = pd.read_csv(vibration_file,sep=',',header=0)

#PLotting v(M-H) vs v(M-O)
EH2 = -(Data.H2-4.477)/2
EO2 = -(Data.O2-5.1158)/2
p = np.polyfit(EH2,EO2, 1) 
m1 = p[0]
b1=p[1]
R2 = 1 - sum(((m1*EH2+b1)-EO2)**2)/sum((EO2-np.mean(EO2))**2) 
plt.figure(1)
plt.figure(figsize=(14,10))
plt.plot(np.array([np.min(EH2),np.max(EH2)]), m1*np.array([np.min(EH2),np.max(EH2)])+b1,'-b')
plt.plot(EH2, EO2,'ro')
plt.xlabel('${\Delta}$E$_{H}$ (eV)',size=32)
plt.ylabel('${\Delta}$E$_{O}$ (eV)',size=32)
plt.legend(['${\Delta}$E$_{O}$=%s${\Delta}$E$_{H}$ + %s eV' %(round(m1,2),(round(b1,2)))],loc=2,prop={'size':32},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([1.97,3.2])
plt.ylim([2.4,6])
plt.gcf().subplots_adjust(bottom=0.15)
mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(Data.Metal)):
    Marker.append(''.join(Data.Metal[i]))
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(EH2, EO2, Data.Metal):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, )
plt.show()
print(R2)
