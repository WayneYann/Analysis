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
from ase import Atoms
from scipy.interpolate import interp1d
import scipy as stats
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['legend.numpoints'] = 1
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 20
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Data.csv')
#with open(vibration_file,'rb') as csvfile:
# newfile = csv.reader(csvfile, delimiter=',')
Data = pd.read_csv(vibration_file,sep=',',header=0)
cols = Data.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
Data.columns = cols
vibration_file2 = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/OHx_and_NHx_condensed_data_v05.csv')
Data2 = pd.read_csv(vibration_file2,sep=',',header=0)
#Available Adsorbates: OHx, NHx, CHx, O2, OOH, N2, CO

Ad1 = 'N'; Ad2 = 'NH'; Surf = 110; NN = 2; BadMetals = ['Cr','V']
freqs = Data.Zfrequency; masses = Data.Adsorbate_mass 
#freqs = Data.frequency_corrected; masses = Data.Mr 

"""Dont change below this line"""
#predicted and actual slope values
v1 = [];E1 = [];L1 = [];v2 = [];E2 = [];L2 = []; M1 = []; M2 = []; G1 = []; G2 = []
MetalLabels=[]; MetalLabels2=[]
DataPoints = len(Data.frequency_corrected)
for i in range(0,DataPoints):
    if (Data.Surface[i] == Surf and Data.Substrate[i] not in BadMetals and Data.NN[i] == NN):
        if Data.Adsorbate[i] ==Ad1:
            MetalLabels.append(Data.Substrate[i])
            v1.append(freqs[i])
            E1.append(Data.Eads[i])
            L1.append(Data.mindistance[i])
            M1.append(masses[i])
            G1.append(Data.vibGibbs[i])
        elif Data.Adsorbate[i] ==Ad2:
            MetalLabels2.append(Data.Substrate[i])
            v2.append(freqs[i])
            E2.append(Data.Eads[i])
            L2.append(Data.mindistance[i])
            M2.append(masses[i])
            G2.append(Data.vibGibbs[i])
for i in range(0,len(MetalLabels)):
    if MetalLabels[i] <> MetalLabels2[i]:
        MetalLabels2 = MetalLabels2[:i] + ['nan'] + MetalLabels2[i:]
        E2 = E2[:i]  + [float('nan')] + E2[i:]
        v2 = v2[:i] + [float('nan')] + v2[i:]
        L2 = L2[:i] + [float('nan')] + L2[i:]
        M2 = M2[:i] + [float('nan')] + M2[i:]
        G2 = G2[:i] + [float('nan')] + G2[i:]


v1 = np.array(v1); v2=np.array(v2); E1 = np.array(E1); E2 = np.array(E2); L1 = np.array(L1); L2 = np.array(L2); M1 = np.array(M1); M2 = np.array(M2); G1 = np.array(G1); G2 = np.array(G2)
idx = np.isfinite(v1) & np.isfinite(v2)
lAB = np.array(L1)[idx]/np.array(L2)[idx]
rAB = np.mean(lAB)
mAB = np.mean(M1[idx]/M2[idx])
v1 = v1[idx]; v2 = v2[idx]; MetalLabels = np.array(MetalLabels)[idx]
E1 = E1[idx]; E2 = E2[idx]; G1 = G1[idx]; G2 = G2[idx]
L1 = L1[idx]; L2 = L2[idx]
PtI = list(MetalLabels).index('Pt')
mat.rc('text', usetex = False)

p = np.polyfit(E1,(1/L2**4),1)
LE = 1/(p[1]**(1/4))
if p[1] <=0:
   LE = float('inf')

res = (p[0]*E1+p[1])-1/L2**4
sres = np.std(res,ddof=2)
sx = np.std(E1,ddof=1)
SE = sres/(sx*(len(E1)-1)**0.5)
SLEinv4 = sres*(np.mean(E1**2)/sum((E1-np.mean(E1))**2))**0.5
sLE = np.abs(LE**5/4*SLEinv4)

plt.figure(1)
plt.figure(figsize=(6,4))
p = np.polyfit(E1, E2, 1)
m1 = p[0]
b1=p[1]
R2 = 1 - sum(((m1*E1+b1)-E2)**2)/sum((E2-np.mean(E2))**2)
#plt.figure(figsize=(12,10))
plt.plot(E1, E2,'bo')
plt.plot(np.array([np.min(E1),np.max(E1)]), m1*np.array([np.min(E1),np.max(E1)])+b1,'--b')
#plt.plot(EadsforExp, vMOExp,'g^')
plt.title(r'${\Delta}$E$_{ads} slope$ on 111-fcc sites',size=20)
plt.xlabel('${\Delta}$E$_{ads}$ atomic %s [eV]' %(Ad1),size=24)
plt.ylabel('${\Delta}$E$_{ads}$ %s [eV]' %(Ad2),size=24)
#plt.legend(['RPBE-D3: y=%sx+%s \n $R^{2}$=%s' %(round(m1,2),round(b1,2),round(R2,2))],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)

texts = []
for a, b, s in zip(E1,E2,MetalLabels):
    texts.append(plt.text(a, b, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'},
    arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
print('Energy slope')
print(m1)
print('Energy intercept')
print(b1)
print('Energy on Pt')
print(E1[PtI])
print('LA')
print(L1[PtI])
print('LB')
print(L2[PtI])
print('LE')
print(LE)
print('bE')
print(b1)
print('R^2')
print(R2)
res = (m1*E1+b1)-E2
sres = np.std(res,ddof=2)
sx = np.std(E1,ddof=1)
srAB = np.std(lAB, ddof=1)/(len(E1)-1)**0.5
SE = sres/(sx*(len(E1)-1)**0.5)
SbE = sres*(np.mean(E1**2)/sum((E1-np.mean(E1))**2))**0.5
print('Standard Error of slope')
print(SE)
print('Standard Error intercept')
print(SbE)
print(SbE/np.mean(E1))
print('rAB')
print(rAB)
print('mAB')
print(mAB)
predslope = (m1*rAB**2*mAB)**0.5
print('predicted slope')
print(predslope)
SEpred = (0.25/m1*rAB**2*mAB*SE**2+m1*mAB*srAB**2)**0.5
print('Standard Error of predicted slope')
print(SEpred)



plt.figure(2)
plt.figure(figsize=(6,4))
p = np.polyfit(v1, v2, 1)
m1 = p[0]
b1=p[1]
R2 = 1 - sum(((m1*v1+b1)-v2)**2)/sum((v2-np.mean(v2))**2)
#plt.figure(figsize=(12,10))
plt.plot(v1, v2,'bo')
plt.plot(np.array([np.min(v1),np.max(v1)]), m1*np.array([np.min(v1),np.max(v1)])+b1,'--b')
#plt.plot(EadsforExp, vMOExp,'g^')
plt.title(r'Frequency',size=20)
plt.xlabel('v%s cm-1' %(Ad1),size=24)
plt.ylabel('v%s cm-1' %(Ad2),size=24)
#plt.legend(['RPBE-D3: y=%sx+%s \n $R^{2}$=%s' %(round(m1,2),round(b1,2),round(R2,2))],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)

texts = []
for a, b, s in zip(v1,v2,MetalLabels):
    texts.append(plt.text(a, b, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'},
    arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
print('Frequency slope')
print(m1)
print('Frequency intercept')
print(b1)
print('R^2')
print(R2)
res = (m1*v1+b1)-v2
sres = np.std(res,ddof=2)
sx = np.std(v1,ddof=1)
SE = sres/(sx*(len(v1)-1)**0.5)
print('Standard Error')
print(SE)
#talpha = stats.t.ppf(1-.025,len(v1)-2)
print('Error')
Error = predslope- m1
print(Error)
print('Sum of standard errors')
SumErrors = SEpred+SE
print(SumErrors)
print('95% CI')
#print(talpha*SE)



"""
G1 = E2
#G2 = v2
plt.figure(3)
plt.figure(figsize=(6,4))
p = np.polyfit(G1, G2, 1)
m1 = p[0]
b1=p[1]
R2 = 1 - sum(((m1*G1+b1)-G2)**2)/sum((G2-np.mean(G2))**2)
#plt.figure(figsize=(12,10))
plt.plot(G1, G2,'bo')
plt.plot(np.array([np.min(G1),np.max(G1)]), m1*np.array([np.min(G1),np.max(G1)])+b1,'--b')
#plt.plot(EadsforExp, vMOExp,'g^')
plt.title(r'Gvib on 111-fcc sites',size=20)
plt.xlabel('${\Delta}$E$_{ads}$ %s [eV]' %(Ad2),size=24)
plt.ylabel('Gvib %s [eV]' %(Ad2),size=24)
#plt.legend(['RPBE-D3: y=%sx+%s \n $R^{2}$=%s' %(round(m1,2),round(b1,2),round(R2,2))],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)

texts = []
for a, b, s in zip(G1,G2,MetalLabels):
    texts.append(plt.text(a, b, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'},
    arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
print('Gibbs slope')
print(m1)
print('R^2')
print(R2)
"""