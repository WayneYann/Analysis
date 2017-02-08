# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:20:17 2016
　
@author: lansford
　
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import adjustText
import os
import pandas as pd
from ase.thermochemistry import IdealGasThermo
from ase import Atoms
from scipy.optimize import curve_fit
import pylab
import matplotlib.gridspec as gridspec
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['legend.numpoints'] = 1
mat.rcParams['lines.linewidth'] = 5
mat.rcParams['lines.markersize'] = 20

#Results
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Data_v02.csv')
Data = pd.read_csv(vibration_file,sep=',',header=0)
cols = Data.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
Data.columns = cols
c = 29979245800 #cm/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules

NN = 1; Metals = ['Au','Ag','Cu','Pt','Ni','Rh','Ir','Pd']
freqs = Data.Zfrequency; masses = Data.Adsorbate_mass 
DataPoints = len(freqs)

idx100 = Data.index[(Data.Surface == 100) & Data.Substrate.isin(Metals) & (Data.NN ==1)]
MetalLabels = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='O')]

O = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='O')])**2
OH = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='OH')])**2
N = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='N')])**2
NH = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='NH')])**2
NH2 = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='NHH')])**2
C = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='C')])**2
CH = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CH')])**2
CH2 = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CHH')])**2
CH3 = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CHHH')])**2
CH2CH3 = (Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CCHHHHH')])**2
C2H4 = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CCHHHH')]

O111 = (Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate=='O')])**2
OH111 = (Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate=='OH')])**2
O110 = (Data.Zfrequency[(Data.Surface == 110) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate=='O')])**2
OH110 = (Data.Zfrequency[(Data.Surface == 110) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate=='OH')])**2

Metals2 = ['Au','Ag','Cu','Pt','Rh','Ir','Pd']
CH3110 = (Data.Zfrequency[(Data.Surface == 110) & Data.Substrate.isin(Metals2) & (Data.NN ==1) & (Data.Adsorbate=='CHHH')])**2
CH2CH3110 = (Data.Zfrequency[(Data.Surface == 110) & Data.Substrate.isin(Metals2) & (Data.NN ==1) & (Data.Adsorbate=='CCHHHHH')])**2

CH3111 = (Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate=='CHHH')])**2
CH2CH3111 = (Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate=='CCHHHHH')])**2

C2H4Eads = Data.Eads[Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')]
C2H4vib = Data.Zfrequency[Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')]
C2H4Metal = list(Data.Substrate[Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')])
C2H4Surface = list(Data.Surface[Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')])

C2H4Eads100 = Data.Eads[(Data.Surface == 100) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')]
C2H4vib100 = Data.Zfrequency[(Data.Surface == 100) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')]

C2H4Eads110 = Data.Eads[(Data.Surface == 110) & (Data.Substrate <> 'Ir') & Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')]
C2H4vib110 = Data.Zfrequency[(Data.Surface == 110) & (Data.Substrate <> 'Ir') & Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')]

C2H4Eads111 = Data.Eads[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')]
C2H4vib111 = Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='CCHHHH')]

C2H4Labels = []
for i in range(0,len(C2H4Metal)):
    C2H4Labels.append(C2H4Metal[i]+'('+str(C2H4Surface[i])+')')

plt.figure(1)
plt.figure(figsize=(12,10))
plt.plot(C2H4Eads111, C2H4vib111,'rs')
plt.plot(C2H4Eads110, C2H4vib110,'g^')
plt.plot(C2H4Eads100, C2H4vib100,'bo')
plt.xlabel('${\Delta}$E$_{ads,C_{2}H_{4}}$ [eV]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,C_{2}H_{4}} [cm^{-1}]$',size=32)
plt.legend(['111','110','100']
,loc=1,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([-1.9,-0.2])
plt.ylim([45,460])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
for x, y, s in zip(C2H4Eads, C2H4vib, C2H4Metal):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#plt.figure(2)
#plt.figure(figsize=(12,10))
mat.rcParams['lines.markersize'] = 10
G = gridspec.GridSpec(2, 2)
#plt.figure(figsize=(15,10.5))
plt.figure(figsize=(15,15))
ax1 = plt.subplot(G[0,1])
pNH = np.polyfit(N, NH, 1)
pNH2 = np.polyfit(N, NH2, 1)
plt.plot(N, NH,'g^')
plt.plot(N, NH2,'bo')
plt.xlabel(r'$\mathbf{\nu}_{\perp,N}$ [$cm^{-1}$]',size=16)
plt.ylabel(r'$\mathbf{\nu}_{\perp,NHx}$ [$cm^{-1}$]',size=16)
plt.plot(np.array([np.min(N),np.max(N)]), pNH[0]*np.array([np.min(N),np.max(N)])+pNH[1],'-g')
plt.plot(np.array([np.min(N),np.max(N)]), pNH2[0]*np.array([np.min(N),np.max(N)])+pNH2[1],'-b')
plt.legend([r'$\mathbf{\nu}_{\perp,NH}$=%.2f$\mathbf{\nu}_{\perp,N}$ + %.0f $cm^{-1}$' %(pNH[0],pNH[1])
,r'$\mathbf{\nu}_{\perp,NH_{2}}$=%.2f$\mathbf{\nu}_{\perp,N}$ + %.0f $cm^{-1}$' %(pNH2[0],pNH2[1])]
,loc=2,prop={'size':12},frameon=False)
plt.xticks(size=14)
plt.yticks(size=14)
#plt.xlim([480**2,1020**2])
#plt.ylim([380**2,860**2])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
M2 = np.concatenate((MetalLabels,MetalLabels))
v11 = np.concatenate((N,N))
v22 = np.concatenate((NH,NH2))
for x, y, s in zip(v11, v22, M2):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=16, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'})
texts.append(ax1.text(0.01,0.95,'(b)',size=12,name ='Calibri',fontweight='bold',transform=ax1.transAxes))

#plt.figure(3)
#plt.figure(figsize=(12,10))
mat.rcParams['lines.markersize'] = 10
ax2 = plt.subplot(G[0,0])
pCH = np.polyfit(C, CH, 1)
pCH2 = np.polyfit(C, CH2, 1)
pCH3 = np.polyfit(C, CH3, 1)
plt.plot(C, CH,'g^')
plt.plot(C, CH2,'bo')
plt.plot(C, CH3,'rs')
plt.xlabel(r'$\mathbf{\nu}_{\perp,C}$ [$cm^{-1}$]',size=16)
plt.ylabel(r'$\mathbf{\nu}_{\perp,CHx}$ [$cm^{-1}$]',size=16)
plt.plot(np.array([np.min(C),np.max(C)]), pCH[0]*np.array([np.min(C),np.max(C)])+pCH[1],'-g')
plt.plot(np.array([np.min(C),np.max(C)]), pCH2[0]*np.array([np.min(C),np.max(C)])+pCH2[1],'-b')
plt.plot(np.array([np.min(C),np.max(C)]), pCH3[0]*np.array([np.min(C),np.max(C)])+pCH3[1],'-r')
plt.legend([r'$\mathbf{\nu}_{\perp,CH}$=%.2f$\mathbf{\nu}_{\perp,C}$ + %.0f $cm^{-1}$' %(pCH[0],pCH[1])
,r'$\mathbf{\nu}_{\perp,CH_{2}}$=%.2f$\mathbf{\nu}_{\perp,C}$ + %.0f $cm^{-1}$' %(pCH2[0],pCH2[1])
,r'$\mathbf{\nu}_{\perp,CH_{3}}$=%.2f$\mathbf{\nu}_{\perp,C}$ + %.0f $cm^{-1}$' %(pCH3[0],pCH3[1])]
,loc=2,prop={'size':12},frameon=False)
plt.xticks(size=14)
plt.yticks(size=14)
#plt.xlim([530,1020])
#plt.ylim([280,1020])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
M2 = np.concatenate((MetalLabels,MetalLabels,MetalLabels))
v11 = np.concatenate((C,C,C))
v22 = np.concatenate((CH,CH2,CH3))
for x, y, s in zip(v11, v22, M2):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=16, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'})
texts.append(ax2.text(0.01,0.95,'(a)',size=12,name ='Calibri',fontweight='bold',transform=ax2.transAxes))
#plt.show()

#plt.figure(4)
#plt.figure(figsize=(33,12))
mat.rcParams['lines.markersize'] = 16
ax3 = plt.subplot(G[1,:])
pOH = np.polyfit(O, OH, 1)
pOH111 = np.polyfit(O111, OH111, 1)
pOH110 = np.polyfit(O110, OH110, 1)
plt.plot(O, OH,'g^')
plt.plot(O111, OH111,'bo')
plt.plot(O110, OH110,'rs')
plt.xlabel(r'$\mathbf{\nu}_{\perp,O}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,OH}$ [$cm^{-1}$]',size=32)
plt.plot(np.array([np.min(O),np.max(O)]), pOH[0]*np.array([np.min(O),np.max(O)])+pOH[1],'-g')
plt.plot(np.array([np.min(O111),np.max(O111)]), pOH111[0]*np.array([np.min(O111),np.max(O111)])+pOH111[1],'-b')
plt.plot(np.array([np.min(O110),np.max(O110)]), pOH110[0]*np.array([np.min(O110),np.max(O110)])+pOH110[1],'-r')
plt.legend([r'$\mathbf{\nu}_{\perp,OH-100}$=%.2f$\mathbf{\nu}_{\perp,O-100}$ + %.0f $cm^{-1}$' %(pOH[0],pOH[1])
,r'$\mathbf{\nu}_{\perp,OH-111}$=%.2f$\mathbf{\nu}_{\perp,O-111}$ + %.0f $cm^{-1}$' %(pOH111[0],pOH111[1])
,r'$\mathbf{\nu}_{\perp,OH-110}$=%.2f$\mathbf{\nu}_{\perp,O-110}$ + %.0f $cm^{-1}$' %(pOH110[0],pOH110[1])]
,loc=2,prop={'size':18},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
#plt.xlim([480,810])
#plt.ylim([340,600])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
M2 = np.concatenate((MetalLabels,MetalLabels,MetalLabels))
v11 = np.concatenate((O,O111,O110))
v22 = np.concatenate((OH,OH111,OH110))
for x, y, s in zip(v11, v22, M2):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
texts.append(ax3.text(0.01,0.92,'(c)',size=18,name ='Calibri',fontweight='bold',transform=ax3.transAxes))
plt.show() 

plt.figure(5)
plt.figure(figsize=(20,10))
pCH2CH3 = np.polyfit(CH3, CH2CH3, 1)
pCH2CH3111 = np.polyfit(CH3111, CH2CH3111, 1)
pCH2CH3110 = np.polyfit(CH3110, CH2CH3110, 1)
plt.plot(CH3, CH2CH3,'g^')
plt.plot(CH3111, CH2CH3111,'bo')
plt.plot(CH3110, CH2CH3110,'rs')
plt.xlabel(r'$\mathbf{\nu}_{\perp,CH_{3}}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,CH_{2}CH_{3}}$ [$cm^{-1}$]',size=32)
plt.plot(np.array([np.min(CH3),np.max(CH3)]), pCH2CH3[0]*np.array([np.min(CH3),np.max(CH3)])+pCH2CH3[1],'-g')
plt.plot(np.array([np.min(CH3111),np.max(CH3111)]), pCH2CH3111[0]*np.array([np.min(CH3111),np.max(CH3111)])+pCH2CH3111[1],'-b')
plt.plot(np.array([np.min(CH3110),np.max(CH3110)]), pCH2CH3110[0]*np.array([np.min(CH3110),np.max(CH3110)])+pCH2CH3110[1],'-r')
plt.legend([r'$\mathbf{\nu}_{\perp,CH_{2}CH_{3}-100}$=%.2f$\mathbf{\nu}_{\perp,CH_{3}-100}$ + %.0f $cm^{-1}$' %(pCH2CH3[0],pCH2CH3[1])
,r'$\mathbf{\nu}_{\perp,CH_{2}CH_{3}-111}$=%.2f$\mathbf{\nu}_{\perp,CH_{3}-111}$ + %.0f $cm^{-1}$' %(pCH2CH3111[0],pCH2CH3111[1])
,r'$\mathbf{\nu}_{\perp,CH_{2}CH_{3}-110}$=%.2f$\mathbf{\nu}_{\perp,CH_{3}-110}$ + %.0f $cm^{-1}$' %(pCH2CH3110[0],pCH2CH3110[1])]
,loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([350,520])
plt.ylim([340,520])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
M2 = np.concatenate((MetalLabels,MetalLabels,MetalLabels))
v11 = np.concatenate((CH3,CH3111,CH3110))
v22 = np.concatenate((CH2CH3,CH2CH3111,CH2CH3110))
for x, y, s in zip(v11, v22, M2):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show() 


Metals = ['Au','Ag','Cu','Pt','Ni','Rh','Ir']
def func(E,a,b):
    if a < 0:
        return (a*E+b)**0.5
    else:
        return 1e15
plt.figure(6)
plt.figure(figsize=(20,10))
Eads = Data.Eads[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
vAHx = Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
Mass = Data.Adsorbate_mass[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
ppAtop, pcov = curve_fit(func,Eads,vAHx,p0=np.array([-200000,-80000]),maxfev = 50000)
Eatop = np.linspace(np.min(Eads),np.max(Eads),num=100)
vatop = (ppAtop[0]*Eatop+ppAtop[1])**0.5
plt.plot(Eads,vAHx,'ob')
Eads = Data.Eads[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['CH','CHH','CHHH']))]
vAHx = Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['CH','CHH','CHHH']))]
finite=np.isfinite(vAHx)
ppHollow, pcov = curve_fit(func,Eads[finite],vAHx[finite],p0=np.array([-200000,-80000]),maxfev = 50000)
Ehollow = np.linspace(np.min(Eads),np.max(Eads),num=100)
vhollow = (ppHollow[0]*Eatop+ppHollow[1])**0.5
Mass = Data.Adsorbate_mass[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['CH','CHH','CHHH']))]
plt.plot(Eads,vAHx,'sg')
plt.plot(Eatop,vatop,'-b')
plt.plot(Ehollow,vhollow,'--g')
plt.xlabel(r'$\mathbf{\Delta}E_{ads,CHx}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,CHx}$ [$cm^{-1}$]',size=32)
plt.xticks(size=28)
plt.yticks(size=28)
plt.legend([r'$\mathbf{\nu}_{atop}$=$\sqrt{%.0f\mathbf{\Delta}E_{ads} - %.0f\: cm^{-2}}$' %(ppAtop[0],abs(ppAtop[1]))
,r'$\mathbf{\nu}_{hollow}$=$\sqrt{%.0f\mathbf{\Delta}E_{ads} - %.0f\: cm^{-2}}$' %(ppHollow[0],abs(ppHollow[1]))]
,loc=2,prop={'size':24},frameon=False)
plt.xlim([-7,-0.5])
plt.ylim([200,1400])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.show()



Metals = ['Au','Ag','Cu','Pt','Ni','Rh','Ir','Pd']
plt.figure(7)
plt.figure(figsize=(12,10))
Eads = Data.Eads[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
vAHx = Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
Mass = Data.Adsorbate_mass[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
plt.plot(Eads,vAHx,'sb')

Eads = Data.Eads[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
vAHx = Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
Mass = Data.Adsorbate_mass[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
plt.plot(Eads,vAHx,'sg')
plt.xlabel(r'$\mathbf{\Delta}E_{ads,OHx}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,OHx}$ [$cm^{-1}$]',size=32)
plt.xticks(size=28)
plt.yticks(size=28)
plt.legend(['atop','fcc'],loc=2,prop={'size':24},frameon=False)
plt.xlim([-5.3,-1.5])
#plt.ylim([230,800])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.show()

plt.figure(8)
plt.figure(figsize=(12,10))
Eads = Data.Eads[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
vAHx = Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
Mass = Data.Adsorbate_mass[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
plt.plot(Eads,vAHx,'sb')

Eads = Data.Eads[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
vAHx = Data.Zfrequency[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
Mass = Data.Adsorbate_mass[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
plt.plot(Eads,vAHx,'sg')
plt.xlabel(r'$\mathbf{\Delta}E_{ads,NHx}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,NHx}$ [$cm^{-1}$]',size=32)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([-5.5,0])
plt.ylim([290,1100])
plt.legend(['atop','fcc'],loc=2,prop={'size':24},frameon=False)
plt.xlim([-5.3,0])
plt.ylim([290,1100])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.show()


Metals = ['Au','Ag','Cu','Pt','Ni','Rh','Ir']
plt.figure(9)
plt.figure(figsize=(12,10))

Eads = Data.Eads[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
vAHx = Data.Xfrequency[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
Mass = Data.Adsorbate_mass[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
plt.plot(Eads,vAHx,'ob')

Eads = Data.Eads[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
vAHx = Data.Yfrequency[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
Mass = Data.Adsorbate_mass[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
plt.plot(Eads,vAHx,'^b')

Eads = Data.Eads[(Data.RorIx == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
vAHx = Data.Xfrequency[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
Mass = Data.Adsorbate_mass[(Data.RorIx == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
plt.plot(Eads,vAHx,'og')

Eads = Data.Eads[(Data.RorIy == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
vAHx = Data.Yfrequency[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
Mass = Data.Adsorbate_mass[(Data.RorIy == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['C','CH','CHH','CHHH']))]
plt.plot(Eads,vAHx,'^g')
plt.xlabel(r'$\mathbf{\Delta}E_{ads,CHx}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\parallel,CHx}$ [$cm^{-1}$]',size=32)
plt.xticks(size=28)
plt.yticks(size=28)
plt.legend(['atop-X','atop-Y','fcc-X','fcc-y'],loc=1,prop={'size':24},frameon=False)
plt.xlim([-7.2,-1])
#plt.ylim([200,1100])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.show()

Metals = ['Au','Ag','Cu','Pt','Ni','Rh','Ir','Pd']
plt.figure(10)
plt.figure(figsize=(12,10))

Eads = Data.Eads[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
vAHx = Data.Xfrequency[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
Mass = Data.Adsorbate_mass[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
plt.plot(Eads,vAHx,'ob')

Eads = Data.Eads[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
vAHx = Data.Yfrequency[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
Mass = Data.Adsorbate_mass[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['O','OH']))]
plt.plot(Eads,vAHx,'^b')

Eads = Data.Eads[(Data.RorIx == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
vAHx = Data.Xfrequency[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
Mass = Data.Adsorbate_mass[(Data.RorIx == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
plt.plot(Eads,vAHx,'og')

Eads = Data.Eads[(Data.RorIy == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
vAHx = Data.Yfrequency[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
Mass = Data.Adsorbate_mass[(Data.RorIy == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['O','OH']))]
plt.plot(Eads,vAHx,'^g')

plt.xlabel(r'$\mathbf{\Delta}E_{ads,OHx}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\parallel,OHx}$ [$cm^{-1}$]',size=32)
plt.xticks(size=28)
plt.yticks(size=28)
plt.legend(['atop-X','atop-Y','fcc-X','fcc-y'],loc=1,prop={'size':24},frameon=False)
plt.xlim([-5.3,-1.5])
#plt.ylim([230,800])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.show()

plt.figure(11)
plt.figure(figsize=(12,10))
Eads = Data.Eads[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
vAHx = Data.Xfrequency[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
Mass = Data.Adsorbate_mass[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
plt.plot(Eads,vAHx,'ob')

Eads = Data.Eads[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
vAHx = Data.Yfrequency[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
Mass = Data.Adsorbate_mass[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==1) & (Data.Adsorbate.isin(['N','NH','NHH']))]
plt.plot(Eads,vAHx,'^b')

Eads = Data.Eads[(Data.RorIx == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
vAHx = Data.Xfrequency[(Data.RorIx == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
Mass = Data.Adsorbate_mass[(Data.RorIx == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
plt.plot(Eads,vAHx,'og')

Eads = Data.Eads[(Data.RorIy == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
vAHx = Data.Yfrequency[(Data.RorIy == 'f') & (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
Mass = Data.Adsorbate_mass[(Data.RorIy == 'f') &  (Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3) & (Data.Adsorbate.isin(['N','NH','NHH']))]
plt.plot(Eads,vAHx,'^g')

plt.xlabel(r'$\mathbf{\Delta}E_{ads,NHx}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\parallel,NHx}$ [$cm^{-1}$]',size=32)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([-5.3,-0.8])
plt.ylim([0,520])
plt.legend(['atop-X','atop-Y','fcc-X','fcc-y'],loc=1,prop={'size':24},frameon=False)
#plt.xlim([-5.5,0])
#plt.ylim([290,1100])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.show()


#O2 vs O
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Data_v02.csv')
Data = pd.read_csv(vibration_file,sep=',',header=0)
cols = Data.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
Data.columns = cols
c = 29979245800 #cm/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules

NN = 1; Metals = ['Au','Ag','Cu','Pt','Ni','Rh','Ir','Pd']
freqs = Data.Zfrequency; masses = Data.Adsorbate_mass 
DataPoints = len(freqs)

idx100 = Data.index[(Data.Surface == 100) & Data.Substrate.isin(Metals) & (Data.NN ==1)]
MetalLabels = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='O')]
OHollow = Data.Zfrequency[(Data.Surface == 111) & (Data.NN ==3) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='O')]
O2 = Data.Zfrequency[(Data.Surface == 111) & (Data.NN ==2) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='OO')]
O2Rotation = Data.Rotation[(Data.Surface == 111) & (Data.NN ==2) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='OO')]
O2Stretch = Data.High_Frequency[(Data.Surface == 111) & (Data.NN ==2) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='OO')]

OHollow = (Data.Zfrequency[Data.Surface.isin([111,100]) & Data.NN.isin([3,2,4]) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='O')])**2
O2 = (Data.Zfrequency[Data.Surface.isin([111,100]) & (Data.NN ==2) & (Data.NN2 ==4) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='OO')])**2
O2Rotation = (Data.Rotation[Data.Surface.isin([111,100]) & (Data.NN ==2) & (Data.NN2 ==4) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='OO')])**2
O2Stretch = (Data.High_Frequency[Data.Surface.isin([111,100]) & (Data.NN ==2) & (Data.NN2 ==4) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='OO')])**2
MetalLabels = Data.Substrate[Data.Surface.isin([111,100]) & (Data.NN ==2) & (Data.NN2 ==4) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='OO')]
Surface = Data.Surface[Data.Surface.isin([111,100]) & (Data.NN ==2) & (Data.NN2 ==4) & Data.Substrate.isin(Metals) & (Data.Adsorbate=='OO')]
Marker = []
MetalLabels = np.array(MetalLabels)
Surface = np.array(Surface)
for i in range(0,len(MetalLabels)):
    Marker.append(MetalLabels[i] + '(' + str(Surface[i]) + ')')

mat.rcParams['lines.markersize'] = 15
plt.figure(12,dpi=1200)
#plt.figure(figsize=(15,10.5))
plt.figure(figsize=(10,15))
ax1 = plt.subplot(2,1,2)
pO2 = np.polyfit(OHollow, O2, 1)
pO2Rotation = np.polyfit(OHollow, O2Rotation, 1)
pO2Stretch = np.polyfit(OHollow, O2Stretch, 1)
plt.plot(OHollow, O2,'bo')
plt.plot(OHollow, O2Rotation,'g^')
plt.plot(OHollow, O2Stretch,'rs')
plt.xlabel(r'$\mathbf{\nu}_{\perp,O}$ [$cm^{-1}$]',size=16)
plt.ylabel(r'$\mathbf{\nu}_{O_{2}}$/$\mathbf{\delta}_{O_{2}}$ [$cm^{-1}$]',size=16)
plt.plot(np.array([np.min(OHollow),np.max(OHollow)]), pO2[0]*np.array([np.min(OHollow),np.max(OHollow)])+pO2[1],'-b')
plt.plot(np.array([np.min(OHollow),np.max(OHollow)]), pO2Rotation[0]*np.array([np.min(OHollow),np.max(OHollow)])+pO2Rotation[1],'-g')
plt.plot(np.array([np.min(OHollow),np.max(OHollow)]), pO2Stretch[0]*np.array([np.min(OHollow),np.max(OHollow)])+pO2Stretch[1],'-r')
plt.legend([r'$\mathbf{\nu}_{O_{2},\perp}$=%.2f$\mathbf{\nu}_{\perp,O}$ + %.0f $cm^{-1}$' %(pO2[0],pO2[1])
,r'$\mathbf{\delta}_{O_{2},\perp}$=%.2f$\mathbf{\nu}_{\perp,O}$ + %.0f $cm^{-1}$' %(pO2Rotation[0],pO2Rotation[1])
,r'$\mathbf{\nu}_{O_{2},O-O}$=%.2f$\mathbf{\nu}_{\perp,O}$ + %.0f $cm^{-1}$' %(pO2Stretch[0],pO2Stretch[1])]
,loc=2,prop={'size':16},frameon=False)
plt.xticks(size=14)
plt.yticks(size=14)
#plt.xlim([190,550])
#plt.ylim([100,1650])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
#plt.annotate('(b)',xy=(1,1),textcoords='figure fraction')
#M2 = np.concatenate((Marker,Marker,Marker))
#v11 = np.concatenate((OHollow,OHollow,OHollow))
#v22 = np.concatenate((O2,O2Rotation,O2Stretch))
for x, y, s in zip(OHollow, O2, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=16, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
texts.append(ax1.text(0.9,0.9,'(b)',size=16,name ='Calibri',fontweight='bold',transform=ax1.transAxes))
#plt.show()

"""Experimental Data"""
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/O_H_N_Quals_Data.csv')
Metal_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,0]
Data_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,1:]
Data = np.genfromtxt(vibration_file, delimiter=',')[1:,1:]
Metal_Info = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[1:,0]

#Calculating reduced mass
MassH = 1.00794
MassO = 16
MassN = 14
vMO = Data[:,0]**2
vMO2 = Data[:,3]**2
vOOsite1 = Data[:,6]**2
vOOsite2 = Data[:,9]**2
vMN = Data[:,12]**2
DataLabels = []
vMH = Data[:,15]**2
vMHhorlit = Data[:,16]**2
for i in range(0,len(Data_Labels)):
    DataLabels.append(Data_Labels[i])
    DataLabels.append(i)

#Plotting v(M-O) vs v(M-O2) vibrations
OHollow = vMO
O2 = vMO2
O2Rotation = vOOsite2
O2Stretch = vOOsite1
idx = np.isfinite(vMO) & np.isfinite(vMO2)
pO2 = np.polyfit(vMO[idx], vMO2[idx], 1)
idx = np.isfinite(vMO) & np.isfinite(vOOsite2)
pO2Rotation = np.polyfit(vMO[idx], vOOsite2[idx], 1)
idx = np.isfinite(vMO) & np.isfinite(vOOsite1)
pO2Stretch = np.polyfit(vMO[idx], vOOsite1[idx], 1)

#plt.figure(13)
#plt.figure(figsize=(15,10.5))
ax2 = plt.subplot(2,1,1)
plt.plot(OHollow, O2,'bo')
plt.plot(OHollow, O2Rotation,'g^')
plt.plot(OHollow, O2Stretch,'rs')
plt.xlabel(r'$\mathbf{\nu}_{\perp,O}$ [$cm^{-1}$]',size=16)
plt.ylabel(r'$\mathbf{\nu}_{O_{2}}$/$\mathbf{\delta}_{O_{2}}$ [$cm^{-1}$]',size=16)
plt.plot(np.array([np.min(OHollow),np.max(OHollow)]), pO2[0]*np.array([np.min(OHollow),np.max(OHollow)])+pO2[1],'-b')
plt.plot(np.array([np.min(OHollow),np.max(OHollow)]), pO2Rotation[0]*np.array([np.min(OHollow),np.max(OHollow)])+pO2Rotation[1],'-g')
plt.plot(np.array([np.min(OHollow),np.max(OHollow)]), pO2Stretch[0]*np.array([np.min(OHollow),np.max(OHollow)])+pO2Stretch[1],'-r')
plt.legend([r'$\mathbf{\nu}_{O_{2},\perp}$=%.2f$\mathbf{\nu}_{\perp,O}$ + %.0f $cm^{-1}$' %(pO2[0],pO2[1])
,r'$\mathbf{\delta}_{O_{2},\perp}$=%.2f$\mathbf{\nu}_{\perp,O}$ + %.0f $cm^{-1}$' %(pO2Rotation[0],pO2Rotation[1])
,r'$\mathbf{\nu}_{O_{2},O-O}$=%.2f$\mathbf{\nu}_{\perp,O}$ + %.0f $cm^{-1}$' %(pO2Stretch[0],pO2Stretch[1])]
,loc=2,prop={'size':16},frameon=False)
plt.xticks(size=14)
plt.yticks(size=14)
#plt.xlim([240,700])
#plt.ylim([190,1400])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
plt.annotate('(a)',xy=(1,1),textcoords='figure fraction')
M2 = np.concatenate((MetalLabels,MetalLabels,MetalLabels))
v11 = np.concatenate((OHollow,OHollow,OHollow))
v22 = np.concatenate((O2,O2Rotation,O2Stretch))

Metal_Info = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[1:,0]
Marker = []
for i in [i for i, x in enumerate(vMO2) if x]:
    if Metal_Info[i] == 'Ag(100)' or Metal_Info[i] == 'Ag(110)':
        Marker.append(Metal_Info[i])
    else:
        Marker.append(Metal_Info[i])
for i in [i for i, x in enumerate(vOOsite1) if x]:
    if Metal_Info[i] == 'Ag(100)' or Metal_Info[i] == 'Ag(110)':        
        Marker.append(Metal_Info[i])
    else:
        Marker.append(Metal_Info[i])
for i in [i for i, x in enumerate(vOOsite2) if x]:
    if Metal_Info[i] == 'Ag(100)' or Metal_Info[i] == 'Ag(110)':
        Marker.append(Metal_Info[i])
    else:
        Marker.append(Metal_Info[i])

for x, y, s in zip(v11, v22, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=16, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
texts.append(ax2.text(0.9,0.9,'(a)',size=16,name='Calibri',fontweight='bold',transform=ax2.transAxes))
plt.show()

vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/O_H_N_Quals_Data.csv')
Data = pd.read_csv(vibration_file,sep=',',header=0)
vMO = Data['v(M-O)']**2
vMH = Data['v(M-H)']**2
vMN = Data['v(M-N)']**2
vMH2 = Data['parallel (M-H)']**2
idx1 = np.isfinite(vMO) & np.isfinite(vMH)
p = np.polyfit(vMO[idx1],vMH[idx1], 1) 
m1 = p[0]
b1=p[1]
idx2 = np.isfinite(vMO) & np.isfinite(vMN)
p = np.polyfit(vMO[idx2],vMN[idx2], 1) 
m2 = p[0]
b2=p[1]
idx3 = np.isfinite(vMO) & np.isfinite(vMH2)
p = np.polyfit(vMO[idx3],vMH2[idx3], 1) 
m3 = p[0]
b3=p[1]

plt.figure(14)
plt.figure(figsize=(20,10))
plt.plot(vMO, vMH,'rs')
plt.plot(vMO, vMH2,'bo')
plt.plot(vMO, vMN,'^g')
plt.plot(np.array([np.min(vMO[idx1]),np.max(vMO[idx1])]), m1*np.array([np.min(vMO[idx1]),np.max(vMO[idx1])])+b1,'-r')
plt.plot(np.array([np.min(vMO[idx3]),np.max(vMO[idx3])]), m3*np.array([np.min(vMO[idx3]),np.max(vMO[idx3])])+b3,'-b')
plt.plot(np.array([np.min(vMO[idx2]),np.max(vMO[idx2])]), m2*np.array([np.min(vMO[idx2]),np.max(vMO[idx2])])+b2,'-g')
plt.xlabel(r'$\mathbf{\nu}_{\perp,O}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp/\parallel,N/H}$ [$cm^{-1}$]',size=32)
plt.legend([r'$\mathbf{\nu}_{\perp,H}$=%s$\mathbf{\nu}_{\perp,O}$ + %.0f eV' %(round(m1,2),b1)
,r'$\mathbf{\nu}_{\parallel,H}$=%s$\mathbf{\nu}_{\perp,O}$ + %.0f eV' %(round(m3,2),b3)
,r'$\mathbf{\nu}_{\perp,N}$=%s$\mathbf{\nu}_{\perp,O}$ + %.0f eV' %(round(m2,2),b2)]
,loc=2,prop={'size':28},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
#plt.xlim([270,650])
#plt.ylim([200,1450])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
M3 = np.concatenate((Data['Surface'][idx1],Data.Surface[idx3],Data.Surface[idx2]))
vOO = np.concatenate((vMO[idx1],vMO[idx3],vMO[idx2]))
vO2N = np.concatenate((vMH[idx1],vMH2[idx3],vMN[idx2]))
for x, y, s in zip(vOO, vO2N, M3):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=28, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

"""DFT Frequency Slopes"""

vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Slopes_v03.csv')
Data = pd.read_csv(vibration_file,sep=',',header=0)
cols = Data.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
Data.columns = cols


nS=6;nP=8
rA = Data.rA; rB = Data.rB
s_fracA = Data.s_frac_A; p_fracA = 1-s_fracA
s_fracB = Data.s_frac_B; p_fracB = 1-s_fracB
RnA = (s_fracA*nS*(nS+1)*rA**nP+p_fracA*nP*(nP+1)*rA**nS)/(rA**2*(s_fracA*rA**nP+p_fracA*rA**nS))
RnB = (s_fracB*nS*(nS+1)*rB**nP+p_fracB*nP*(nP+1)*rB**nS)/(rB**2*(s_fracB*rB**nP+p_fracB*rB**nS))
mVpred = Data.mE_no_Disp*Data.MA/Data.MB*(RnB/RnA)
mV = Data.mv
Marker = np.array(Data.Adsorbates)
for i in range(0,len(Marker)):
    Marker[i] = Marker[i].replace("2","$_{2}$")
    Marker[i] = Marker[i].replace("3","$_{3}$")
idx100 = Data.index[(Data.Surface == 100)]
idx110 = Data.index[(Data.Surface == 110)]
idx111 = Data.index[(Data.Surface == 111)]

mV100 = mV[idx100]
mVpred100 = mVpred[idx100]
mV110 = mV[idx110]
mVpred110 = mVpred[idx110]
mV111 = mV[idx111]
mVpred111 = mVpred[idx111]
Marker111 = Marker[idx111]

Marker = np.array(Data.Adsorbates)
for i in range(0,len(Marker)):
    Marker[i] = Marker[i].replace("2","$_{2}$")
    Marker[i] = Marker[i].replace("3","$_{3}$")
Cidx100 = Data.index[(Data.Surface == 100) & (Data['Adsorbates'].str.contains('C'))]
Oidx100 = Data.index[(Data.Surface == 100) & (Data['Adsorbates'].str.contains('O'))]
Nidx100 = Data.index[(Data.Surface == 100) & (Data['Adsorbates'].str.contains('N'))]
Cidx110 = Data.index[(Data.Surface == 110) & (Data['Adsorbates'].str.contains('C'))]
Oidx110 = Data.index[(Data.Surface == 110) & (Data['Adsorbates'].str.contains('O'))]
Nidx110 = Data.index[(Data.Surface == 110) & (Data['Adsorbates'].str.contains('N'))]
Cidx111 = Data.index[(Data.Surface == 111) & (Data['Adsorbates'].str.contains('C'))]
Oidx111 = Data.index[(Data.Surface == 111) & (Data['Adsorbates'].str.contains('O'))]
Nidx111 = Data.index[(Data.Surface == 111) & (Data['Adsorbates'].str.contains('N'))]

CmV100 = mV[Cidx100]
CmVpred100 = mVpred[Cidx100]
CmV110 = mV[Cidx110]
CmVpred110 = mVpred[Cidx110]
CmV111 = mV[Cidx111]
CmVpred111 = mVpred[Cidx111]

NmV100 = mV[Nidx100]
NmVpred100 = mVpred[Nidx100]
NmV110 = mV[Nidx110]
NmVpred110 = mVpred[Nidx110]
NmV111 = mV[Nidx111]
NmVpred111 = mVpred[Nidx111]

OmV100 = mV[Oidx100]
OmVpred100 = mVpred[Oidx100]
OmV110 = mV[Oidx110]
OmVpred110 = mVpred[Oidx110]
OmV111 = mV[Oidx111]
OmVpred111 = mVpred[Oidx111]

plt.figure(15)
plt.figure(figsize=(20,10))
plt.plot([0,1],[0,1],'k-')
plt.plot(CmV100,CmVpred100,'bo')
plt.plot(CmV110,CmVpred110,'b^')
plt.plot(CmV111,CmVpred111,'bs')
plt.plot(NmV100,NmVpred100,'go')
plt.plot(NmV110,NmVpred110,'g^')
plt.plot(NmV111,NmVpred111,'gs')
plt.plot(OmV100,OmVpred100,'ro')
plt.plot(OmV110,OmVpred110,'r^')
plt.plot(OmV111,OmVpred111,'rs')
plt.xlabel(r'Slope Fit to DFT Data',size=28)
plt.ylabel(r'Predicted Slope',size=28)
plt.legend(['Parity']
,loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.xlim([0.08,0.9])
plt.ylim([0.080,0.9])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []

for x, y, s in zip(mV110, mVpred110, Marker111):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'})
plt.show()

#vibrational Gibbs effects
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Data_v02.csv')
Data = pd.read_csv(vibration_file,sep=',',header=0)
cols = Data.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
Data.columns = cols

NN = 3; Metals = ['Au','Ag','Cu','Pt','Ni','Rh','Ir','Pd']
freqs = Data.Zfrequency; masses = Data.Adsorbate_mass 
DataPoints = len(freqs)

idx111 = Data.index[(Data.Surface == 111) & Data.Substrate.isin(Metals) & (Data.NN ==3)]
MetalLabels = list(Data.Substrate[Data.index.isin(idx111) & (Data.Adsorbate=='O')])
EO = np.array(Data.Eads[Data.index.isin(idx111) & (Data.Adsorbate=='O')])
ECO = np.array(Data.Eads[Data.index.isin(idx111) & (Data.Adsorbate=='CO')])
GO = np.array(Data.vibGibbs[Data.index.isin(idx111) & (Data.Adsorbate=='O')])
GCO = np.array(Data.vibGibbs[Data.index.isin(idx111) & (Data.Adsorbate=='CO')])

O2r = (0.5149910055-0.4850089945)*41.3427602451982
O2freq = 0.1904905 #eV
COr = (0.51388868234-0.486111317659)*41.3427602451982
COfreq = 0.260119879
O2 = Atoms('O2')
O2[1].z = O2r
CO = Atoms('CO')
CO[1].z = COr
O2Thermo = IdealGasThermo([O2freq],'linear',atoms=O2,symmetrynumber=2,spin=1,natoms=2)
O2gibbs = O2Thermo.get_gibbs_energy(298,101325);
COThermo = IdealGasThermo([COfreq],'linear',atoms=CO,symmetrynumber=1,spin=0,natoms=2)
COgibbs = COThermo.get_gibbs_energy(298,101325);

PtI = MetalLabels.index('Pt')
AgI = MetalLabels.index('Ag')
CuI = MetalLabels.index('Cu')
PdI = MetalLabels.index('Pd')

#redefining O2 eenergy change from molecular oxygen
EO = EO+9.5884297/2-1.9555224
delGCOPt = GCO[PtI]+ECO[PtI]-COgibbs

delGCOAg = GCO[AgI]+ECO[AgI]-COgibbs
delGCOPtAg = GCO[PtI]+ECO[AgI]-COgibbs
delGOPt = GO[PtI]+EO[PtI]-O2gibbs/2
delGOAg = GO[AgI]+EO[AgI]-O2gibbs/2
delGOPtAg = GO[PtI]+EO[AgI]-O2gibbs/2


KcoPt = np.exp(-delGCOPt/(kB*298/JeV))
KcoAg = np.exp(-delGCOAg/(kB*298/JeV))
KcoPtAg = np.exp(-delGCOPtAg/(kB*298/JeV))
KoPt = np.exp(-delGOPt/(kB*298/JeV))
KoAg = np.exp(-delGOAg/(kB*298/JeV))
KoPtAg = np.exp(-delGOPtAg/(kB*298/JeV))

idx1 = np.isfinite(EO) & np.isfinite(GO)
idx2 = np.isfinite(ECO) & np.isfinite(GCO)
p = np.polyfit(EO[idx1],GO[idx1], 1) 
mO = p[0]
bO=p[1]
p = np.polyfit(ECO[idx2],GCO[idx2], 1) 
mCO = p[0]
bCO=p[1]

mat.rcParams['lines.markersize'] = 10
pCO = 10**(-8)

xxmin = -2.4
xxmax = -1.4
yymin = -2
yymax = -1.5
xx = pylab.linspace(xxmin, xxmax, 200)
yy = pylab.linspace(yymin, yymax, 200)
zz = pylab.zeros([len(yy), len(xx)])
def thetaCO(EadsO, EadsCO):
    delGCO = EadsCO-COgibbs +mCO*EadsCO+bCO
    Kco = np.exp(-1*delGCO/(kB*298/JeV))
    delGO = EadsO-O2gibbs/2 +mO*EadsO+bO
    Ko = np.exp(-1*delGO/(kB*298/JeV))
    theta = Kco*pCO/(1+Kco*pCO+(Ko*(1-pCO))**0.5)
    return theta

for i in xrange(len(xx)):
    for j in xrange(len(yy)):
        zz[j, i] = thetaCO(xx[i], yy[j])
pylab.figure(16)
pylab.pcolor(xx, yy, zz)
pylab.plot(EO,ECO,'wo')
pylab.ylabel('CO Eads [eV]')
pylab.xlabel('O Eads [eV]')
pylab.colorbar()
pylab.ylim([yymin,yymax])
pylab.xlim([xxmin,xxmax])
texts = []
for x, y, s in zip(EO+0.01, ECO+0.02, MetalLabels):
    texts.append(pylab.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri',color='w'))
pylab.show()

def thetaCO2(EadsO, EadsCO):
    delGCO = EadsCO-COgibbs + GCO[AgI]
    Kco = np.exp(-1*delGCO/(kB*298/JeV))
    delGO = EadsO-O2gibbs/2 +GO[AgI]
    Ko = np.exp(-1*delGO/(kB*298/JeV))
    theta = Kco*pCO/(1+Kco*pCO+(Ko*(1-pCO))**0.5)
    return theta

for i in xrange(len(xx)):
    for j in xrange(len(yy)):
        zz[j, i] = thetaCO2(xx[i], yy[j])
pylab.figure(17)
pylab.pcolor(xx, yy, zz)
pylab.plot(EO,ECO,'wo')
pylab.ylim([yymin,yymax])
pylab.xlim([xxmin,xxmax])
pylab.ylabel('CO Eads [eV]')
pylab.xlabel('O Eads [eV]')
pylab.colorbar()
texts = []
for x, y, s in zip(EO+0.01, ECO+0.02, MetalLabels):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri',color='w'))
pylab.show()

for i in xrange(len(xx)):
    for j in xrange(len(yy)):
        zz[j, i] = thetaCO2(xx[i], yy[j]) - thetaCO(xx[i], yy[j])
pylab.figure(18)
pylab.pcolor(xx, yy, zz)
pylab.plot(EO,ECO,'wo')
pylab.ylabel('CO Eads [eV]')
pylab.xlabel('O Eads [eV]')
pylab.colorbar()
pylab.ylim([yymin,yymax])
pylab.xlim([xxmin,xxmax])
texts = []
for x, y, s in zip(EO+0.01, ECO+0.02, MetalLabels):
    texts.append(pylab.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri',color='w'))
pylab.show()