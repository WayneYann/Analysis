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

#PLotting v(M-O) vs v(M-OH)
p = np.polyfit(Data.H2/2,Data.O2/2, 1) 
m1 = p[0]
b1=p[1]
plt.figure(1)
plt.figure(figsize=(14,10))
plt.plot(np.array([np.min(Data.H2/2),np.max(Data.H2/2)]), m1*np.array([np.min(Data.H2/2),np.max(Data.H2/2)])+b1,'-b')
plt.plot(Data.H2/2, Data.O2/2,'ro')
plt.xlabel('${\Delta}$E$_{ads,H}$ (eV)',size=32)
plt.ylabel('${\Delta}$E$_{ads,O}$ (eV)',size=32)
plt.legend(['${\Delta}$E$_{ads,O}$=%s${\Delta}$E$_{ads,H}$ - %s eV' %(round(m1,2),abs(round(b1,2)))],loc=2,prop={'size':32},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([-1.1,0.3])
plt.ylim([-3.6,0])
plt.gcf().subplots_adjust(bottom=0.15)
mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(Data.Metal)):
    Marker.append(''.join(Data.Metal[i]))
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(Data.H2/2, Data.O2/2, Data.Metal):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, )
plt.show()

vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Data.csv')
Data = pd.read_csv(vibration_file,sep=',',header=0)

vO = [];GvibO = [];SvibO = [];vH = [];GvibH = [];SvibH = [];
MetalLabels=[]
DataPoints = len(Data.frequency_corrected)
for i in range(0,DataPoints):
    if (Data.Surface[i] == 111 and Data.Substrate[i]<>'Al'):
        if Data.Adsorbate[i] =='O':
            MetalLabels.append(Data.Substrate[i])
            vO.append(Data.frequency_corrected[i])
            GvibO.append(Data.vibGibbs[i])
            SvibO.append(Data.vibEntropy[i])
        elif Data.Adsorbate[i] =='H':
            vH.append(Data.frequency_corrected[i])
            GvibH.append(Data.vibGibbs[i])
            SvibH.append(Data.vibEntropy[i])
        

p = np.polyfit(vH,vO, 1) 
m1 = p[0]
b1=p[1]
plt.figure(1)
plt.figure(figsize=(14,10))
plt.plot(np.array([np.min(vH),np.max(vH)]), m1*np.array([np.min(vH),np.max(vH)])+b1,'-b')
plt.plot(vH, vO,'ro')
plt.xlabel(r'$\mathbf{\nu}_{\perp (M-H)}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp (M-O)}$ [$cm^{-1}$]',size=32)
plt.legend([r'$\mathbf{\nu}_{\perp (M-O)}$=%s$\mathbf{\nu}_{\perp (M-H)}$ - %s $cm^{-1}$' %(round(m1,2),abs(round(b1,2)))],loc=2,prop={'size':32},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([780,1150])
plt.ylim([330,550])
plt.gcf().subplots_adjust(bottom=0.15)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(MetalLabels)):
    Marker.append(MetalLabels)
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(vH, vO, MetalLabels):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, )
plt.show()

p = np.polyfit(GvibH,GvibO, 1) 
m1 = p[0]
b1=p[1]
plt.figure(1)
plt.figure(figsize=(14,10))
plt.plot(np.array([np.min(GvibH),np.max(GvibH)]), m1*np.array([np.min(GvibH),np.max(GvibH)])+b1,'-b')
plt.plot(GvibH, GvibO,'ro')
plt.xlabel('$G_{vib,H}$ at 298K (eV)',size=32)
plt.ylabel('$G_{vib,O}$ at 298K (eV)',size=32)
plt.legend(['$G_{vib,O}$=%s$G_{vib,H}$ - %s eV' %(round(m1,2),abs(round(b1,2)))],loc=2,prop={'size':32},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([.115,.17])
plt.ylim([0.032,0.065])
plt.gcf().subplots_adjust(bottom=0.15)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(MetalLabels)):
    Marker.append(MetalLabels)
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(GvibH, GvibO, MetalLabels):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, )
plt.show()

p = np.polyfit(SvibH,SvibO, 1) 
m1 = p[0]
b1=p[1]
plt.figure(1)
plt.figure(figsize=(14,10))
plt.plot(np.array([np.min(SvibH),np.max(SvibH)]), m1*np.array([np.min(SvibH),np.max(SvibH)])+b1,'-b')
plt.plot(SvibH, SvibO,'ro')
plt.xlabel('$S_{vib,H}$ at 298K (10$^{-5}$eV/K)',size=32)
plt.ylabel('$S_{vib,O}$ at 298K (10$^{-4}$eV/K)',size=32)
plt.legend(['$S_{vib,O}$=%s$S_{vib,H}$ + %.2e eV/K' %(round(m1,2),b1)],loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.ticklabel_format(style='scientific')
plt.ticklabel_format(scilimits=(1,2))
plt.ticklabel_format(useOffset=False)
plt.xlim([1.7*10**(-5),4.5*10**(-5)])
plt.ylim([1.21*10**(-4),1.9*10**(-4)])
plt.gcf().subplots_adjust(bottom=0.15)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(MetalLabels)):
    Marker.append(MetalLabels)
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(SvibH, SvibO, MetalLabels):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, )
plt.show()