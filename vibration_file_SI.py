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
from statmech import vibgibbs2
import statmech
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['legend.numpoints'] = 1
mat.rcParams['lines.linewidth'] = 3
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

idx100 = Data.index[(Data.Surface == 100) & Data.Substrate.isin(Metals)]
idx110 = Data.index[(Data.Surface == 110) & Data.Substrate.isin(Metals)]
idx111 = Data.index[(Data.Surface == 111) & Data.Substrate.isin(Metals)]
MetalLabels = Data.Substrate[Data.index.isin(idx111) & (Data.Adsorbate=='O')]

OZ111 = Data.Zfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='O') & (Data.NN==3)]
OHZ111 = Data.Zfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='OH') & (Data.NN==3)]
OZ110 = Data.Zfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='O') & (Data.NN==2)]
OHZ110 = Data.Zfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='OH') & (Data.NN==2)]
OZ100 = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='O') & (Data.NN==4)]
OHZ100 = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='OH') & (Data.NN==4)]
OZ100B = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='O') & (Data.NN==2)]
OHZ100B = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='OH') & (Data.NN==2)]

OX111 = Data.Xfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='O') & (Data.NN==3)]
OHX111 = Data.Xfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='OH') & (Data.NN==3)]
OX110 = Data.Xfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='O') & (Data.NN==2)]
OHX110 = Data.Xfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='OH') & (Data.NN==2)]
OY110 = Data.Yfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='O') & (Data.NN==2)]
OHY110 = Data.Yfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='OH') & (Data.NN==2)]
OX100 = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='O') & (Data.NN==4)]
OHX100 = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='OH') & (Data.NN==4)]
OX100B = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='O') & (Data.NN==2)]
OHX100B = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='OH') & (Data.NN==2)]
OY100B = Data.Yfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='O') & (Data.NN==2)]
OHY100B = Data.Yfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='OH') & (Data.NN==2)]

Metals111 = Data.Substrate[Data.index.isin(idx111) & (Data.Adsorbate=='O') & (Data.NN==3)]
Metals110 = Data.Substrate[Data.index.isin(idx110) & (Data.Adsorbate=='O') & (Data.NN==2)]
Metals100 = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='O')]
Metals100H = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='O') & (Data.NN==4)]
Metals100B = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='O') & (Data.NN==4)]
plt.figure(1)
plt.figure(figsize=(12,10))
plt.plot(OZ111, OHZ111,'gs')
plt.plot(OZ110, OHZ110,'bo')
plt.plot(OZ100, OHZ100,'g^')
plt.plot(OZ100B, OHZ100B,'bv')
plt.xlabel(r'$\mathbf{\nu}_{\perp,O}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,OH}$ [$cm^{-1}$]',size=32)
plt.legend(['111-Hollow','110-Bridge','100-Hollow','100-Bridge']
,loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)

plt.xlim([190,550])
plt.ylim([210,450])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
Metals = pd.concat([Metals111,Metals110,Metals100])
OZ = pd.concat([OZ111,OZ110,OZ100,OZ100B])
OHZ = pd.concat([OHZ111,OHZ110,OHZ100,OHZ100B])
for x, y, s in zip(OZ, OHZ, Metals):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
p = np.polyfit(OZ, OHZ, 1)
m1 = p[0]
b1=p[1]
print('slope')
print(m1)
print('intercept')
print(b1)

plt.figure(2)
plt.figure(figsize=(12,10))
plt.plot(OX111, OHX111,'gs')
plt.plot(pd.concat([OX110,OY110]), pd.concat([OHX110,OHY110]),'bo')
plt.plot(OX100, OHX100,'g^')
plt.plot(pd.concat([OX100B,OY100B]), pd.concat([OHX100B,OHY100B]),'bv')
plt.xlabel(r'$\mathbf{\nu}_{\parallel,O}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\parallel,OH}$ [$cm^{-1}$]',size=32)
plt.legend(['111-Hollow','110-Bridge','100-Hollow','100-Bridge']
,loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)

plt.xlim([50,500])
plt.ylim([10,350])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
Metals = pd.concat([Metals111,Metals110,Metals110,Metals100H,Metals100B,Metals100B])
OX = pd.concat([OX111,OX110,OY110,OX100,OX100B,OY100B])
OHX = pd.concat([OHX111,OHX110,OHY110,OHX100,OHX100B,OHY100B])
for x, y, s in zip(OX, OHX, Metals):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
p = np.polyfit(OX, OHX, 1)
m1 = p[0]
b1=p[1]
print('slope')
print(m1)
print('intercept')
print(b1)

MetalLabels = Data.Substrate[Data.index.isin(idx111) & (Data.Adsorbate=='N')]

NZ111 = Data.Zfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='N') & (Data.NN==3)]
NHZ111 = Data.Zfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='NH') & (Data.NN==3)]
NZ110 = Data.Zfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='N') & (Data.NN==2)]
NHZ110 = Data.Zfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='NH') & (Data.NN==2)]
NZ100 = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='N') & (Data.NN==4)]
NHZ100 = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='NH') & (Data.NN==4)]
NZ100B = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='N') & (Data.NN==2)]
NHZ100B = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='NH') & (Data.NN==2)]

NX111 = Data.Xfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='N') & (Data.NN==3)]
NHX111 = Data.Xfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='NH') & (Data.NN==3)]
NX110 = Data.Xfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='N') & (Data.NN==2)]
NHX110 = Data.Xfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='NH') & (Data.NN==2)]
NY110 = Data.Yfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='N') & (Data.NN==2)]
NHY110 = Data.Yfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='NH') & (Data.NN==2)]
NX100 = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='N') & (Data.NN==4)]
NHX100 = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='NH') & (Data.NN==4)]
NX100B = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='N') & (Data.NN==2)]
NHX100B = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='NH') & (Data.NN==2)]
NY100B = Data.Yfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='N') & (Data.NN==2)]
NHY100B = Data.Yfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='NH') & (Data.NN==2)]

Metals111 = Data.Substrate[Data.index.isin(idx111) & (Data.Adsorbate=='N') & (Data.NN==3)]
Metals110 = Data.Substrate[Data.index.isin(idx110) & (Data.Adsorbate=='N') & (Data.NN==2)]
Metals100 = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='N')]
Metals100H = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='N') & (Data.NN==4)]
Metals100B = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='N') & (Data.NN==4)]
plt.figure(3)
plt.figure(figsize=(12,10))
plt.plot(NZ111, NHZ111,'gs')
plt.plot(NZ110, NHZ110,'bo')
plt.plot(NZ100, NHZ100,'g^')
plt.plot(NZ100B, NHZ100B,'bv')
plt.xlabel(r'$\mathbf{\nu}_{\perp,N}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,NH}$ [$cm^{-1}$]',size=32)
plt.legend(['111-Hollow','110-Bridge','100-Hollow','100-Bridge']
,loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)

plt.xlim([150,700])
plt.ylim([245,650])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
Metals = pd.concat([Metals111,Metals110,Metals100])
NZ = pd.concat([NZ111,NZ110,NZ100,NZ100B])
NHZ = pd.concat([NHZ111,NHZ110,NHZ100,NHZ100B])
for x, y, s in zip(NZ, NHZ, Metals):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Nalibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
p = np.polyfit(NZ, NHZ, 1)
m1 = p[0]
b1=p[1]
print('slope')
print(m1)
print('intercept')
print(b1)

plt.figure(4)
plt.figure(figsize=(12,10))
plt.plot(NX111, NHX111,'gs')
plt.plot(pd.concat([NX110,NY110]), pd.concat([NHX110,NHY110]),'bo')
plt.plot(NX100, NHX100,'g^')
plt.plot(pd.concat([NX100B,NY100B]), pd.concat([NHX100B,NHY100B]),'bv')
plt.xlabel(r'$\mathbf{\nu}_{\parallel,N}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\parallel,NH}$ [$cm^{-1}$]',size=32)
plt.legend(['111-Hollow','110-Bridge','100-Hollow','100-Bridge']
,loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)

#plt.xlim([50,500])
#plt.ylim([10,350])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
Metals = pd.concat([Metals111,Metals110,Metals110,Metals100H,Metals100B,Metals100B])
NX = pd.concat([NX111,NX110,NY110,NX100,NX100B,NY100B])
NHX = pd.concat([NHX111,NHX110,NHY110,NHX100,NHX100B,NHY100B])
for x, y, s in zip(NX, NHX, Metals):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Nalibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
p = np.polyfit(NX, NHX, 1)
m1 = p[0]
b1=p[1]
print('slope')
print(m1)
print('intercept')
print(b1)


MetalLabels = Data.Substrate[Data.index.isin(idx111) & (Data.Adsorbate=='C')]

CZ111 = Data.Zfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='C') & (Data.NN==3)]
CHZ111 = Data.Zfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='CH') & (Data.NN==3)]
CZ110 = Data.Zfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='C') & (Data.NN==2)]
CHZ110 = Data.Zfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='CH') & (Data.NN==2)]
CZ100 = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='C') & (Data.NN==4)]
CHZ100 = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CH') & (Data.NN==4)]
CZ100B = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='C') & (Data.NN==2)]
CHZ100B = Data.Zfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CH') & (Data.NN==2)]

CX111 = Data.Xfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='C') & (Data.NN==3)]
CHX111 = Data.Xfrequency[Data.index.isin(idx111) & (Data.Adsorbate=='CH') & (Data.NN==3)]
CX110 = Data.Xfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='C') & (Data.NN==2)]
CHX110 = Data.Xfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='CH') & (Data.NN==2)]
CY110 = Data.Yfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='C') & (Data.NN==2)]
CHY110 = Data.Yfrequency[Data.index.isin(idx110) & (Data.Adsorbate=='CH') & (Data.NN==2)]
CX100 = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='C') & (Data.NN==4)]
CHX100 = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CH') & (Data.NN==4)]
CX100B = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='C') & (Data.NN==2)]
CHX100B = Data.Xfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CH') & (Data.NN==2)]
CY100B = Data.Yfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='C') & (Data.NN==2)]
CHY100B = Data.Yfrequency[Data.index.isin(idx100) & (Data.Adsorbate=='CH') & (Data.NN==2)]

Metals111 = Data.Substrate[Data.index.isin(idx111) & (Data.Adsorbate=='C') & (Data.NN==3)]
Metals110 = Data.Substrate[Data.index.isin(idx110) & (Data.Adsorbate=='C') & (Data.NN==2)]
Metals100 = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='C')]
Metals100H = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='C') & (Data.NN==4)]
Metals100B = Data.Substrate[Data.index.isin(idx100) & (Data.Adsorbate=='C') & (Data.NN==4)]
plt.figure(5)
plt.figure(figsize=(12,10))
plt.plot(CZ111, CHZ111,'gs')
plt.plot(CZ110, CHZ110,'bo')
plt.plot(CZ100, CHZ100,'g^')
plt.plot(CZ100B, CHZ100B,'bv')
plt.xlabel(r'$\mathbf{\nu}_{\perp,C}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp,CH}$ [$cm^{-1}$]',size=32)
plt.legend(['111-Hollow','110-Bridge','100-Hollow','100-Bridge']
,loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)

plt.xlim([225,650])
plt.ylim([295,750])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
Metals = pd.concat([Metals111,Metals110,Metals100])
CZ = pd.concat([CZ111,CZ110,CZ100,CZ100B])
CHZ = pd.concat([CHZ111,CHZ110,CHZ100,CHZ100B])
for x, y, s in zip(CZ, CHZ, Metals):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
p = np.polyfit(CZ, CHZ, 1)
m1 = p[0]
b1=p[1]
print('slope')
print(m1)
print('intercept')
print(b1)

plt.figure(6)
plt.figure(figsize=(12,10))
plt.plot(CX111, CHX111,'gs')
plt.plot(pd.concat([CX110,CY110]), pd.concat([CHX110,CHY110]),'bo')
plt.plot(CX100, CHX100,'g^')
plt.plot(pd.concat([CX100B,CY100B]), pd.concat([CHX100B,CHY100B]),'bv')
plt.xlabel(r'$\mathbf{\nu}_{\parallel,C}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\parallel,CH}$ [$cm^{-1}$]',size=32)
plt.legend(['111-Hollow','110-Bridge','100-Hollow','100-Bridge']
,loc=2,prop={'size':24},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)

#plt.xlim([50,500])
#plt.ylim([10,350])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
Metals = pd.concat([Metals111,Metals110,Metals110,Metals100H,Metals100B,Metals100B])
CX = pd.concat([CX111,CX110,CY110,CX100,CX100B,CY100B])
CHX = pd.concat([CHX111,CHX110,CHY110,CHX100,CHX100B,CHY100B])
for x, y, s in zip(CX, CHX, Metals):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()
p = np.polyfit(CX, CHX, 1)
m1 = p[0]
b1=p[1]
print('slope')
print(m1)
print('intercept')
print(b1)

vCu = np.array([1802.14,254.26,220.5,218.67,107.87,104.36])
vPt = np.array([1731.9,349.97,322.15,319.11,161.45,160.01])
vAg = np.array([1887.8,184.08,136.43,135.39,28.845,17.72])
vRh = np.array([1730.3,340.98,267.67,266.45,146.44,145.42])
vAl = np.array([1591.86,193.18,169.06,168.5,71.08,70.13])
vIr = np.array([1704.79,358.09,255.558,248.724,150.6,142.86])
vNi = np.array([1751.16,348.58,276.47,273.86,138.5,134.67])
vPd = np.array([1750.07,351.6,348.65,323.906,164.009,162.306])

mat.rcParams['lines.markersize'] = 20
Pairs = ['Cu','Pt','Ag','Rh','Al','Ir','Ni','Pd']
vibgibbs = [statmech.vibgibbs(vCu,298),statmech.vibgibbs(vPt,298),statmech.vibgibbs(vAg,298),statmech.vibgibbs(vRh,298),statmech.vibgibbs(vAl,298),statmech.vibgibbs(vIr,298),statmech.vibgibbs(vNi,298),statmech.vibgibbs(vPd,298)]
vibgibbsTaylor = [statmech.vibgibbsTaylor(vCu,298,vCu),statmech.vibgibbsTaylor(vPt,298,vCu),statmech.vibgibbsTaylor(vAg,298,vCu),statmech.vibgibbsTaylor(vRh,298,vCu),statmech.vibgibbsTaylor(vAl,298,vCu),statmech.vibgibbsTaylor(vIr,298,vCu),statmech.vibgibbsTaylor(vNi,298,vCu),statmech.vibgibbsTaylor(vPd,298,vCu)]
plt.figure(18)
plt.figure(figsize=(12,10))
#plt.figure(figsize=(12,10))
plt.plot(vibgibbs, vibgibbsTaylor,'bo')
plt.plot(vibgibbs, vibgibbs,'k-')
#plt.title('Comparing The Predicted Frequency and Energy Slopes',size=20)
plt.xlabel('Gibbs Energy (eV)',size=24)
plt.ylabel('Taylor Expansion (eV)',size=24)
plt.legend(['Metal','Parity Line'],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)
plt.xlim([-0.025,0.15])
plt.ylim([-0.025,0.16])
texts = []
for a, b, s in zip(vibgibbs,vibgibbsTaylor,Pairs):
    texts.append(plt.text(a, b, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'},
    arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

v50 = np.linspace(25,100,10)
v100 = np.linspace(50,200,10)
v500 = np.linspace(250,1000,10)
v1000 = np.linspace(500,2000,10)
v2000 = np.linspace(1000,4000,10)

v50Gibbs50 = vibgibbs2(v50,50)
v50Gibbs150 = vibgibbs2(v50,150)
v50Gibbs300 = vibgibbs2(v50,300)
v50Gibbs500 = vibgibbs2(v50,500)
v50Gibbs1000 = vibgibbs2(v50,1000)

plt.figure(19)
plt.figure(figsize=(12,10))
plt.plot(v50,v50Gibbs50,'b-')
plt.plot(v50,v50Gibbs150,'g--')
plt.plot(v50,v50Gibbs300,'ro')
plt.plot(v50,v50Gibbs500,'c^')
plt.plot(v50,v50Gibbs1000,'ms')
plt.legend(['50K','150K','300K','500K','1000K']
,loc=4,prop={'size':20},frameon=False)
plt.ylabel('G$_{vib}$ [eV]',size=20)
plt.xlabel(r'$\mathbf{\nu}$ [$cm^{-1}$]',size=20)
plt.xticks(size=20)
plt.yticks(size=20)
plt.show

v100Gibbs50 = vibgibbs2(v100,50)
v100Gibbs150 = vibgibbs2(v100,150)
v100Gibbs300 = vibgibbs2(v100,300)
v100Gibbs500 = vibgibbs2(v100,500)
v100Gibbs1000 = vibgibbs2(v100,1000)

plt.figure(20)
plt.figure(figsize=(12,10))
plt.plot(v100,v100Gibbs50,'b-')
plt.plot(v100,v100Gibbs150,'g--')
plt.plot(v100,v100Gibbs300,'ro')
plt.plot(v100,v100Gibbs500,'c^')
plt.plot(v100,v100Gibbs1000,'ms')
plt.legend(['50K','150K','300K','500K','1000K']
,loc=4,prop={'size':20},frameon=False)
plt.ylabel('G$_{vib}$ [eV]',size=20)
plt.xlabel(r'$\mathbf{\nu}$ [$cm^{-1}$]',size=20)
plt.xticks(size=20)
plt.yticks(size=20)
plt.show

v500Gibbs50 = vibgibbs2(v500,50)
v500Gibbs150 = vibgibbs2(v500,150)
v500Gibbs300 = vibgibbs2(v500,300)
v500Gibbs500 = vibgibbs2(v500,500)
v500Gibbs1000 = vibgibbs2(v500,1000)

plt.figure(21)
plt.figure(figsize=(12,10))
plt.plot(v500,v500Gibbs50,'b-')
plt.plot(v500,v500Gibbs150,'g--')
plt.plot(v500,v500Gibbs300,'ro')
plt.plot(v500,v500Gibbs500,'c^')
plt.plot(v500,v500Gibbs1000,'ms')
plt.legend(['50K','150K','300K','500K','1000K']
,loc=4,prop={'size':20},frameon=False)
plt.ylabel('G$_{vib}$ [eV]',size=20)
plt.xlabel(r'$\mathbf{\nu}$ [$cm^{-1}$]',size=20)
plt.xticks(size=20)
plt.yticks(size=20)
plt.show

v1000Gibbs50 = vibgibbs2(v1000,50)
v1000Gibbs150 = vibgibbs2(v1000,150)
v1000Gibbs300 = vibgibbs2(v1000,300)
v1000Gibbs500 = vibgibbs2(v1000,500)
v1000Gibbs1000 = vibgibbs2(v1000,1000)

plt.figure(22)
plt.figure(figsize=(12,10))
plt.plot(v1000,v1000Gibbs50,'b-')
plt.plot(v1000,v1000Gibbs150,'g--')
plt.plot(v1000,v1000Gibbs300,'ro')
plt.plot(v1000,v1000Gibbs500,'c^')
plt.plot(v1000,v1000Gibbs1000,'ms')
plt.legend(['50K','150K','300K','500K','1000K']
,loc=4,prop={'size':20},frameon=False)
plt.ylabel('G$_{vib}$ [eV]',size=20)
plt.xlabel(r'$\mathbf{\nu}$ [$cm^{-1}$]',size=20)
plt.xticks(size=20)
plt.yticks(size=20)
plt.show

v2000Gibbs50 = vibgibbs2(v2000,50)
v2000Gibbs150 = vibgibbs2(v2000,150)
v2000Gibbs300 = vibgibbs2(v2000,300)
v2000Gibbs500 = vibgibbs2(v2000,500)
v2000Gibbs1000 = vibgibbs2(v2000,1000)

plt.figure(23)
plt.plot(v2000,v2000Gibbs50,'b-')
plt.plot(v2000,v2000Gibbs150,'g--')
plt.plot(v2000,v2000Gibbs300,'ro')
plt.plot(v2000,v2000Gibbs500,'c^')
plt.plot(v2000,v2000Gibbs1000,'ms')
plt.legend(['50K','150K','300K','500K','1000K']
,loc=4,prop={'size':20},frameon=False)
plt.ylabel('G$_{vib}$ [eV]',size=20)
plt.xlabel(r'$\mathbf{\nu}$ [$cm^{-1}$]',size=20)
plt.xticks(size=20)
plt.yticks(size=20)
plt.show