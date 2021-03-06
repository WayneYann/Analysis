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
mat.rcParams['legend.numpoints'] = 1
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 20
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Data.csv')
#with open(vibration_file,'rb') as csvfile:
#    newfile = csv.reader(csvfile, delimiter=',')
Data = pd.read_csv(vibration_file,sep=',',header=0)

vibration_file2 = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/OHx_and_NHx_condensed_data_v05.csv')
Data2 = pd.read_csv(vibration_file2,sep=',',header=0)    
#PLotting Eads (M-O) vs frequency
vMODFT = [Data.frequency_corrected[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O' and Data.Substrate[i] <> 'Al')]
EadsODFT = [np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O' and Data.Substrate[i] <> 'Al')]
MetalLabels = [Data.Substrate[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O' and Data.Substrate[i] <> 'Al')]
vMOExp = [Data2['v(M-O)'][i] for i in range(0, len(Data2.Surface)) if (Data2.Cut[i] == 111 and np.isnan(Data2['v(M-O)'][i]) == False and Data2.Surface[i] <> 'Ni(111)*' and Data2.Metal[i] <> 'Al' and Data2.Surface[i] <> 'Cr(111)' and Data2.Metal[i] in Data.Substrate.values)]
Substrates = [Data2.Metal[i] for i in range(0, len(Data2.Surface)) if (Data2.Cut[i] == 111 and np.isnan(Data2['v(M-O)'][i]) == False and Data2.Surface[i] <> 'Ni(111)*' and Data2.Metal[i] <> 'Al' and Data2.Surface[i] <> 'Cr(111)' and Data2.Metal[i] in Data.Substrate.values)]

EadsforExp = []
for i in range(0,len(Substrates)):
    for j in range(0,len(Data.Substrate)):
        if Data.Substrate[j] == Substrates[i] and Data.Surface[j] ==111 and Data.Adsorbate[j] =='O':
            EadsforExp.append(np.double(Data.Eads[j]))
vMODFT = np.array(vMODFT)
EadsODFT = np.array(EadsODFT)
vMOExp = np.array(vMOExp)
EadsforExp = np.array(EadsforExp)

p = np.polyfit(EadsODFT, vMODFT, 1) 
m1 = p[0]
b1=p[1]
RvDFT = 1 - sum(((m1*EadsODFT+b1)-vMODFT)**2)/sum(((m1*EadsODFT+b1)-np.mean(vMODFT))**2)
MAEDFT = np.mean(abs((m1*EadsODFT+b1)-vMODFT))
p = np.polyfit(EadsforExp, vMOExp, 1) 
m2 = p[0]
b2=p[1]
RvExp = 1 - sum(((m1*EadsforExp+b1)-vMOExp)**2)/sum(((m1*EadsforExp+b1)-np.mean(vMOExp))**2)
#f1 = plt.figure(1)
#f1.set_figheight(8)
#f1.set_figwidth(16)
#f1.set_dpi(500)
plt.figure(1)
plt.figure(figsize=(16,8))
plt.plot(np.array([np.min(EadsODFT),np.max(EadsODFT)]), m1*np.array([np.min(EadsODFT),np.max(EadsODFT)])+b1,'--b')
#plt.plot(np.array([np.min(EadsforExp),np.max(EadsforExp)]), m2*np.array([np.min(EadsforExp),np.max(EadsforExp)])+b2,'--g')
plt.plot(EadsODFT, vMODFT,'bo')
#plt.plot(EadsforExp, vMOExp,'g^')
plt.title('Frequency Trends with Binding Energy: 111 FCC metals',size=20, fontweight='bold')
plt.xlabel('${\Delta}$E$_{ads}$ atomic O [eV]',size=24)
plt.ylabel(r'$\mathbf{\nu}_{\perp}(M-O)$ [cm$^{-1}$]',size=24)
plt.legend(['RPBE-D3: y=%sx+%s \n $R^{2}$=%s, MAE = %s cm$^{-1}$' %(round(m1,2),round(b1,2),round(RvDFT,2),round(MAEDFT,2))],loc=1,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)
plt.xlim(-5.1, -2.8)
plt.ylim(280, 570)

mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(MetalLabels)):
    Marker.append(''.join(MetalLabels[i]))
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(EadsODFT, vMODFT, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#PLotting v(M-O) vs v(M-OH)
vMODFT = [Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O')]
EadsODFT = [np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O')]
vMOHDFT = [Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='OH')]
EadsOHDFT = [np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='OH')]
MetalLabels = [Data.Substrate[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O')]

vMODFT = np.array(vMODFT)
EadsODFT = np.array(EadsODFT)
vMOHDFT = np.array(vMOHDFT)
EadsOHDFT = np.array(EadsOHDFT)

p = np.polyfit(vMODFT, vMOHDFT, 1) 
m1 = p[0]
b1=p[1]
Rvfreq = 1 - sum(((m1*vMODFT+b1)-vMOHDFT)**2)/sum(((m1*vMODFT+b1)-np.mean(vMOHDFT))**2)
MAEfreq = np.mean(abs((m1*vMODFT+b1)-vMOHDFT))
p = np.polyfit(EadsODFT, EadsOHDFT, 1) 
m2 = p[0]
b2=p[1]
RvEads = 1 - sum(((m2*EadsODFT+b2)-EadsOHDFT)**2)/sum(((m2*EadsODFT+b2)-np.mean(EadsOHDFT))**2)
MAEEads = np.mean(abs((m2*EadsODFT+b2)-EadsOHDFT))
plt.figure(2)
plt.figure(figsize=(10,10))
plt.plot(vMODFT, vMOHDFT,'b^')
plt.plot(np.array([np.min(vMODFT),np.max(vMODFT)]), m1*np.array([np.min(vMODFT),np.max(vMODFT)])+b1,'--b')
#plt.plot(EadsforExp, vMOExp,'g^')
plt.title(r'$\mathbf{\nu}_{\perp}(M-OHx)$ on 111-fcc sites',size=20)
plt.xlabel(r'$\mathbf{\nu}_{\perp}(M-O)$ [cm$^{-1}$]',size=24)
plt.ylabel(r'$\mathbf{\nu}_{\perp}(M-OH)$ [cm$^{-1}$]',size=24)
plt.legend(['RPBE-D3: y=%sx+%s \n $R^{2}$=%s, MAE=%s cm$^{-1}$' %(round(m1,2),round(b1,2),round(Rvfreq,2),round(MAEfreq,2))],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)

plt.xlim(280, 500)
plt.ylim(220, 400)
mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(MetalLabels)):
    Marker.append(''.join(MetalLabels[i]))
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(vMODFT, vMOHDFT, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

plt.figure(3)
plt.figure(figsize=(10,10))
plt.plot(EadsODFT, EadsOHDFT,'bo')
plt.plot(np.array([np.min(EadsODFT),np.max(EadsODFT)]), m2*np.array([np.min(EadsODFT),np.max(EadsODFT)])+b2,'--b')
#plt.plot(EadsforExp, vMOExp,'g^')
plt.title(r'${\Delta}$E$_{ads}(M-OHx)$ on 111-fcc sites',size=20)
plt.xlabel('${\Delta}$E$_{ads}$ atomic O [eV]',size=24)
plt.ylabel('${\Delta}$E$_{ads}$ OH [eV]',size=24)
plt.legend(['RPBE-D3: y=%sx+%s \n $R^{2}$=%s, MAE=%s eV' %(round(m2,2),round(b2,2),round(RvEads,2),round(MAEEads,2))],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)
plt.xlim(-8.2, -2)
plt.ylim(-4.2, -1.5)
mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(MetalLabels)):
    Marker.append(''.join(MetalLabels[i]))
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(EadsODFT, EadsOHDFT, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
Anharmonic_Data = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Anharmonic_Data.csv')
Data3 = pd.read_csv(Anharmonic_Data,sep=',',header=0)

#PLotting v(M-O) vs v(M-OH)

weinfO = np.array([Data3.weinf[i] for i in range(0,len(Data3.weinf)) if (Data3.Adsorbate[i] =='O')])
wexeO = np.array([Data3.wexeinf[i] for i in range(0,len(Data3.weinf)) if (Data3.Adsorbate[i] =='O')])
MetalO = [Data3.Metal[i] for i in range(0,len(Data3.weinf)) if (Data3.Adsorbate[i] =='O')]
weinfO2 = np.array([Data3.weinf[i] for i in range(0,len(Data3.weinf)) if (Data3.Adsorbate[i] =='O2')])
wexeO2 = np.array([Data3.wexeinf[i] for i in range(0,len(Data3.weinf)) if (Data3.Adsorbate[i] =='O2')])
MetalO2 = [Data3.Metal[i] for i in range(0,len(Data3.weinf)) if (Data3.Adsorbate[i] =='O2')]

p = np.polyfit(weinfO, wexeO, 1) 
m1 = p[0]
b1=p[1]
RvO = 1 - sum(((m1*weinfO+b1)-wexeO)**2)/sum(((m1*weinfO+b1)-np.mean(wexeO))**2)
MAEO = np.mean(abs((m1*weinfO+b1)-wexeO))
p = np.polyfit(weinfO2, wexeO2, 1) 
m2 = p[0]
b2=p[1]
RvO2 = 1 - sum(((m2*weinfO2+b2)-wexeO2)**2)/sum(((m2*weinfO2+b2)-np.mean(wexeO2))**2)

plt.figure(4)
plt.figure(figsize=(14,10))
plt.plot(weinfO, wexeO,'bo')
plt.plot(np.array([np.min(weinfO),np.max(weinfO)]), m1*np.array([np.min(weinfO),np.max(weinfO)])+b1,'--b')
plt.title(r'Scaling of $\mathbf{\omega}_{e}\mathbf{\chi}_{e}$ with $\mathbf{\omega}_{e}$ for Atomic O on 111-Hollow Sites',size=20)
plt.ylabel(r'$\mathbf{\omega}_{e}\mathbf{\chi}_{e}$ [cm$^{-1}$]',size=24)
plt.xlabel(r'$\mathbf{\omega}_{e}$ [cm$^{-1}$]',size=24)
plt.legend([r'$\mathbf{\omega}_{e}\mathbf{\chi}_{e}$=%s$\mathbf{\omega}_{e}$-%s cm$^{-1}$' %(round(m1,4),round(np.sign(b1)*b1,2))+'\n $R^{2}$=%s, MAE=%s cm$^{-1}$' %(round(RvO,2),round(MAEO,2))],loc=2,prop={'size':24})
plt.xticks(size=20)
plt.yticks(size=20)
plt.xlim(300, 500)
plt.ylim(0.15, 0.9)
mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(MetalO)):
    Marker.append(MetalO[i])
for i in range(0,len(MetalO2)):
    Marker.append(MetalO2[i])
mat.rc('text', usetex = False)
weinf = np.concatenate((weinfO,weinfO2))
wexe = np.concatenate((wexeO,wexeO2))
texts = []
for x, y, s in zip(weinf, wexe, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'})
plt.show()