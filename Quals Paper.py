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
from scipy.interpolate import interp1d
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['legend.numpoints'] = 1
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/O_H_N_Quals_Data.csv')
Metal_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,0]
Data_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,1:]
Data = np.genfromtxt(vibration_file, delimiter=',')[1:,1:]
Metal_Info = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[1:,0]

#Calculating reduced mass
MassH = 1.00794
MassO = 16
MassN = 14
vMO = Data[:,0]
vMO2 = Data[:,3]
vOOsite1 = Data[:,6]
vOOsite2 = Data[:,9]
vMN = Data[:,12]
DataLabels = []
vMH = Data[:,15]
vMHhorlit = Data[:,16]
for i in range(0,len(Data_Labels)):
    DataLabels.append(Data_Labels[i])
    DataLabels.append(i)

#Plotting v(M-O) vs v(M-O2) vibrations

idx = np.isfinite(vMO) & np.isfinite(vMO2)
m1,b1 = np.polyfit(vMO[idx], vMO2[idx], 1)
RvO2 = 1 - sum(((m1*vMO[idx]+b1)-vMO2[idx])**2)/sum(((m1*vMO[idx]+b1)-np.mean(vMO2[idx]))**2)
idx = np.isfinite(vMO) & np.isfinite(vOOsite1)
m2,b2 = np.polyfit(vMO[idx], vOOsite1[idx], 1)
ROO1 = 1 - sum(((m2*vMO[idx]+b2)-vOOsite1[idx])**2)/sum(((m2*vMO[idx]+b2)-np.mean(vOOsite1[idx]))**2) 
idx = np.isfinite(vMO) & np.isfinite(vOOsite2)
m3,b3 = np.polyfit(vMO[idx], vOOsite2[idx], 1)
ROO2 = 1 - sum(((m3*vMO[idx]+b3)-vOOsite2[idx])**2)/sum(((m3*vMO[idx]+b3)-np.mean(vOOsite2[idx]))**2)
mat.rc('text', usetex = False)
plt.figure(2)
#plt.figure(figsize=(14,8),dpi=2000)
plt.figure(figsize=(14,8))
plt.plot(vMO,vMO2,'ob',markersize=16)
plt.plot(vMO,vOOsite1,'^g',markersize=16)
plt.plot(vMO,vOOsite2,'^r',markersize=16)
plt.plot(np.array([min(vMO),max(vMO)]), m1*np.array([min(vMO),max(vMO)])+b1,'--b',lw=4)
plt.plot(np.array([min(vMO),max(vMO)]), m2*np.array([min(vMO),max(vMO)])+b2,'--g',lw=4) 
plt.plot(np.array([min(vMO),max(vMO)]), m3*np.array([min(vMO),max(vMO)])+b3,'--r',lw=4)
#plt.title('Experimental frequencies in adsorbed O2 vs O',size=12, fontweight='bold')
#mat.rc('text', usetex = True)
#plt.xlabel(r'\textbf{time}',size=28)
mat.rc('text', usetex=False)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
#plt.xlabel(r'$\mathbf{\nu}_{\perp}(M-O)$ [$cm^{-1}$]',size=28)
mat.rc('text', usetex = False)
plt.ylabel(r'$\mathbf{\nu}_{\perp}(M-O_{2})$ and (O-O) [$cm^{-1}$]',size=28)
plt.legend([r'$\mathbf{\nu}_{\perp}(M-O_{2}$): %sx+%s' %(round(m1,2),int(round(b1,2))), \
'(O-O) site 1: %sx+%s' %(round(m2,2),int(round(b2,2))) ,\
'(O-O) site 2: %sx+%s' %(round(m3,2),int(round(b3,2)))],loc=2,prop={'size':20})
plt.xticks(size=24)
plt.yticks(size=24)

mat.rc('text', usetex = False)
Marker = []


idxvMO2 = np.isfinite(vMO2)
idxvOOsite1 = np.isfinite(vOOsite1)
idxvOOsite2 = np.isfinite(vOOsite2)
#vMOs = np.concatenate((vMO[idxvMO2],vMO[idxvOOsite1],vMO[idxvOOsite2]))
#vMO2s = np.concatenate((vMO2[idxvMO2],vOOsite1[idxvOOsite1],vOOsite2[idxvOOsite2]))
vMOs = np.concatenate((vMO,vMO,vMO))
vMO2s = np.concatenate((vMO2,vOOsite1,vOOsite2))

for i in [i for i, x in enumerate(vMO2) if x]:
    if Metal_Info[i] == 'Ag(001)' or Metal_Info[i] == 'Ag(110)':
        Marker.append(''.join(Metal_Info[i]+'${^{'+str(Data[i,2]).replace(".0","")+'}}$'+'${^{,'+str(Data[i,5]).replace(".0","")+'}}$'))
    else:
        Marker.append(''.join(Metal_Info[i]+'${^{'+str(round(Data[i,2])).replace(".0","")+'}}$'))
for i in [i for i, x in enumerate(vOOsite1) if x]:
    if Metal_Info[i] == 'Ag(001)' or Metal_Info[i] == 'Ag(110)':        
        Marker.append(''.join(Metal_Info[i]+'${^{'+str(Data[i,2]).replace(".0","")+'}}$'+'${^{,'+str(Data[i,8]).replace(".0","")+'}}$'))
    else:
        Marker.append(''.join(Metal_Info[i]+'${^{'+str(Data[i,2]).replace(".0","")+'}}$'))
for i in [i for i, x in enumerate(vOOsite2) if x]:
    if Metal_Info[i] == 'Ag(001)' or Metal_Info[i] == 'Ag(110)':
        Marker.append(''.join(Metal_Info[i]+'${^{'+str(Data[i,2]).replace(".0","")+'}}$'+'${^{,'+str(Data[i,11]).replace(".0","")+'}}$'))
    else:
        Marker.append(''.join(Metal_Info[i]+'${^{'+str(Data[i,2]).replace(".0","")+'}}$'))
mat.rc('text', usetex = False)

texts = []

idx = np.isfinite(vMO2s)

for x, y, s in zip(vMOs[idx], vMO2s[idx], np.array(Marker)[idx]):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=22, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#Plotting v(M-O) vs v(M-N), vMH, and vMH parallel vibrations

idx = np.isfinite(vMO) & np.isfinite(vMN)
m1,b1 = np.polyfit(vMO[idx], vMN[idx], 1)
RvN = 1 - sum(((m1*vMO[idx]+b1)-vMN[idx])**2)/sum(((m1*vMO[idx]+b1)-np.mean(vMN[idx]))**2)
idx = np.isfinite(vMO) & np.isfinite(vMH)
m2,b2 = np.polyfit(vMO[idx], vMH[idx], 1)
RvH = 1 - sum(((m2*vMO[idx]+b2)-vMH[idx])**2)/sum(((m2*vMO[idx]+b2)-np.mean(vMH[idx]))**2) 
idx = np.isfinite(vMO) & np.isfinite(vMHhorlit)
m3,b3 = np.polyfit(vMO[idx], vMHhorlit[idx], 1)
RvHpara = 1 - sum(((m3*vMO[idx]+b3)-vMHhorlit[idx])**2)/sum(((m3*vMO[idx]+b3)-np.mean(vMHhorlit[idx]))**2)
mat.rc('text', usetex = False)
plt.figure(2)
plt.figure(figsize=(14,8),dpi=2000)
#plt.figure(figsize=(14,8))
plt.plot(vMO,vMN,'ob',markersize=16)
plt.plot(vMO,vMH,'^g',markersize=16)
plt.plot(vMO,vMHhorlit,'^r',markersize=16)
plt.plot(np.array([min(vMO),max(vMO)]), m1*np.array([min(vMO),max(vMO)])+b1,'--b',lw=4)
plt.plot(np.array([min(vMO),max(vMO)]), m2*np.array([min(vMO),max(vMO)])+b2,'--g',lw=4) 
plt.plot(np.array([min(vMO),max(vMO)]), m3*np.array([min(vMO),max(vMO)])+b3,'--r',lw=4)
#plt.title('Experimental frequencies in adsorbed O2 vs O',size=12, fontweight='bold')
#mat.rc('text', usetex = True)
#plt.xlabel(r'\textbf{time}',size=28)
mat.rc('text', usetex=False)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
plt.xlabel(r'$\mathbf{\nu}_{\perp}(M-O)$ [$cm^{-1}$]',size=24)
mat.rc('text', usetex = False)
plt.ylabel(r'$\mathbf{\nu}_{\perp}(M-N)$, $\mathbf{\nu}_{\perp}(M-H)$, $\mathbf{\nu}_{\parallel}(M-H)$ [$cm^{-1}$]',size=24)
plt.legend([r'$\mathbf{\nu}_{\perp}(M-N}$): %sx+%s' %(round(m1,2),int(round(b1,2))), \
r'$\mathbf{\nu}_{\perp}(M-H}$): %sx+%s' %(round(m2,2),int(round(b2,2))) ,\
r'$\mathbf{\nu}_{\parallel}(M-H}$): %sx+%s' %(round(m3,2),int(round(b3,2)))],loc=2,prop={'size':20})
plt.xticks(size=24)
plt.yticks(size=24)

mat.rc('text', usetex = False)
Marker = []


idxvMN = np.isfinite(vMN)
idxvMH = np.isfinite(vMH)
idxvMHhorlit = np.isfinite(vMHhorlit)
#vMOs = np.concatenate((vMO[idxvMO2],vMO[idxvOOsite1],vMO[idxvOOsite2]))
#vMO2s = np.concatenate((vMO2[idxvMO2],vOOsite1[idxvOOsite1],vOOsite2[idxvOOsite2]))
vMOs = np.concatenate((vMO[idxvMN],vMO[idxvMH],vMO[idxvMHhorlit]))
vMs = np.concatenate((vMN[idxvMN],vMH[idxvMH],vMHhorlit[idxvMHhorlit]))


for i in [i for i, x in enumerate(idxvMN) if x]:
    Marker.append(''.join(Metal_Info[i]+'${^{'+str(Data[i,2]).replace(".0","")+'}}$'+'${^{,'+str(Data[i,14]).replace(".0","")+'}}$'))

for i in [i for i, x in enumerate(idxvMH) if x]:
    Marker.append(''.join(Metal_Info[i]+'${^{'+str(Data[i,2]).replace(".0","")+'}}$'+'${^{,'+str(Data[i,19]).replace(".0","")+'}}$'))

for i in [i for i, x in enumerate(idxvMHhorlit) if x]:
    Marker.append(''.join(Metal_Info[i]+'${^{'+str(Data[i,2]).replace(".0","")+'}}$'+'${^{,'+str(Data[i,19]).replace(".0","")+'}}$'))

mat.rc('text', usetex = False)

texts = []

idx = np.isfinite(vMs)

for x, y, s in zip(vMOs, vMs, np.array(Marker)):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=22, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

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
vMODFTforExp = []
vMOuncorrected = []
for i in range(0,len(Substrates)):
    for j in range(0,len(Data.Substrate)):
        if Data.Substrate[j] == Substrates[i] and Data.Surface[j] ==111 and Data.Adsorbate[j] =='O':
            EadsforExp.append(np.double(Data.Eads[j]))
            vMODFTforExp.append(Data.frequency_corrected[j])
            vMOuncorrected.append(Data.Zfrequency[j])
vMODFT = np.array(vMODFT)
EadsODFT = np.array(EadsODFT)
vMOExp = np.array(vMOExp)
EadsforExp = np.array(EadsforExp)

p = np.polyfit(EadsODFT, vMODFT, 1) 
m1 = p[0]
b1=p[1]
RvDFT = 1 - sum(((m1*EadsODFT+b1)-vMODFT)**2)/sum(((m1*EadsODFT+b1)-np.mean(vMODFT))**2)
p = np.polyfit(EadsforExp, vMOExp, 1) 
m2 = p[0]
b2=p[1]
RvExp = 1 - sum(((m1*EadsforExp+b1)-vMOExp)**2)/sum(((m1*EadsforExp+b1)-np.mean(vMOExp))**2)
#f1 = plt.figure(1)
#f1.set_figheight(8)
#f1.set_figwidth(16)
#f1.set_dpi(500)
plt.figure(1)
plt.figure(figsize=(12,10))
plt.plot(EadsODFT, vMODFT,'bo')
plt.plot(np.array([np.min(EadsODFT),np.max(EadsODFT)]), m1*np.array([np.min(EadsODFT),np.max(EadsODFT)])+b1,'--b')
#plt.plot(np.array([np.min(EadsforExp),np.max(EadsforExp)]), m2*np.array([np.min(EadsforExp),np.max(EadsforExp)])+b2,'--g')
#plt.plot(EadsforExp, vMOExp,'g^')
#plt.title('Frequency Trends with Binding Energy: 111 FCC metals',size=20, fontweight='bold')
plt.xlabel('${\Delta}$E$_{ads}$ atomic O [eV]',size=24)
plt.ylabel(r'$\mathbf{\nu}_{\perp}(M-O)$ [cm$^{-1}$]',size=24)
plt.legend([r'$\mathbf{\nu}_{\perp,O}$=%s${\Delta}$E$_{ads,O}$+%s $cm^{-1}$' %(round(m1,2),round(b1,2))+'\n $R^{2}$=%s' %(round(RvDFT,2))],loc=1,prop={'size':20})
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
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'})
        #,arrowprops=dict(arrowstyle="-", color='k', lw=2))
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
p = np.polyfit(EadsODFT, EadsOHDFT, 1) 
m2 = p[0]
b2=p[1]
RvEads = 1 - sum(((m2*EadsODFT+b2)-EadsOHDFT)**2)/sum(((m2*EadsODFT+b2)-np.mean(EadsOHDFT))**2)

plt.figure(2)
plt.figure(figsize=(12,10))
plt.plot(vMODFT, vMOHDFT,'b^')
plt.plot(np.array([np.min(vMODFT),np.max(vMODFT)]), m1*np.array([np.min(vMODFT),np.max(vMODFT)])+b1,'--b')
#plt.plot(EadsforExp, vMOExp,'g^')
#plt.title(r'$\mathbf{\nu}_{\perp}(M-OHx)$ on 111-fcc sites',size=20)
plt.xlabel(r'$\mathbf{\nu}_{\perp}(M-O)$ [cm$^{-1}$]',size=24)
plt.ylabel(r'$\mathbf{\nu}_{\perp}(M-OH)$ [cm$^{-1}$]',size=24)
plt.legend([r'$\mathbf{\nu}_{\perp,OH}$=%s$\mathbf{\nu}_{\perp,O}$+%s $cm^{-1}$' %(round(m1,2),round(b1,2))+'\n $R^{2}$=%s' %(round(Rvfreq,2))],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)
plt.xlim(290, 500)
plt.ylim(230, 400)

mat.rc('text', usetex = True)

Marker = []
for i in range(0,len(MetalLabels)):
    Marker.append(''.join(MetalLabels[i]))
mat.rc('text', usetex = False)

texts = []
for x, y, s in zip(vMODFT, vMOHDFT, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'})
plt.show()

plt.figure(3)
plt.figure(figsize=(12,10))
plt.plot(EadsODFT, EadsOHDFT,'bo')
plt.plot(np.array([np.min(EadsODFT),np.max(EadsODFT)]), m2*np.array([np.min(EadsODFT),np.max(EadsODFT)])+b2,'--b')
#plt.plot(EadsforExp, vMOExp,'g^')
plt.title(r'${\Delta}$E$_{ads}(M-OHx)$ on 111-fcc sites',size=20)
plt.xlabel('${\Delta}$E$_{ads}$ atomic O [eV]',size=24)
plt.ylabel('${\Delta}$E$_{ads}$ OH [eV]',size=24)
plt.legend(['RPBE-D3: y=%sx+%s \n $R^{2}$=%s' %(round(m2,2),round(b2,2),round(RvEads,2))],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)

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

#predicted and actual slope values
#PLotting v(M-O) vs v(M-OH)
vMODFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O' and Data.Substrate[i]<>'Ni')])
EadsODFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O' and Data.Substrate[i]<>'Ni')])
LO = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O' and Data.Substrate[i]<>'Ni')])
vMOHDFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='OH' and Data.Substrate[i]<>'Ni')])
EadsOHDFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='OH' and Data.Substrate[i]<>'Ni')])
LOH = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='OH' and Data.Substrate[i]<>'Ni')])
vMNDFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='N' and Data.Substrate[i]<>'Ni')])
EadsNDFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='N' and Data.Substrate[i]<>'Ni')])
LN = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='N' and Data.Substrate[i]<>'Ni')])
vMNHDFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NH' and Data.Substrate[i]<>'Ni')])
EadsNHDFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NH' and Data.Substrate[i]<>'Ni')])
LNH = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NH' and Data.Substrate[i]<>'Ni')])
vMCDFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='C' and Data.Substrate[i]<>'Ni')])
EadsCDFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='C' and Data.Substrate[i]<>'Ni')])
LC = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='C' and Data.Substrate[i]<>'Ni')])
vMCHDFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='CH' and Data.Substrate[i]<>'Ni')])
EadsCHDFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='CH' and Data.Substrate[i]<>'Ni')])
LCH = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='CH' and Data.Substrate[i]<>'Ni')])

vMNHHDFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NHH' and Data.Substrate[i]<>'Ni')])
EadsNHHDFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NHH' and Data.Substrate[i]<>'Ni')])
LNHH = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NHH' and Data.Substrate[i]<>'Ni')])
vMNHDFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NH' and Data.Substrate[i]<>'Ni')])
EadsNHDFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NH' and Data.Substrate[i]<>'Ni')])
LNH = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='NH' and Data.Substrate[i]<>'Ni')])

vMOHHDFT = np.array([Data.Zfrequency[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='OHH' and Data.Substrate[i]<>'Ni')])
EadsOHHDFT = np.array([np.double(Data.Eads[i]) for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='OHH' and Data.Substrate[i]<>'Ni')])
LOHH = np.array([Data.mindistance[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='OHH' and Data.Substrate[i]<>'Ni')])

MetalLabels = [Data.Substrate[i] for i in range(0,len(Data.frequency_corrected)) if (Data.Surface[i] == 111 and  Data.Adsorbate[i] =='O' and Data.Substrate[i]<>'Ni')]

vMODFT = np.array(vMODFT)
EadsODFT = np.array(EadsODFT)
vMOHDFT = np.array(vMOHDFT)
EadsOHDFT = np.array(EadsOHDFT)

