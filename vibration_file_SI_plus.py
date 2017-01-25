# -*- coding: utf-8 -*-
"""
Nreated on Thu Jan 05 11:59:18 2017

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

NN = 1; Metals = ['Au','Ag','Nu','Pt','Ni','Rh','Ir','Pd']
freqs = Data.Zfrequency; masses = Data.Adsorbate_mass 
DataPoints = len(freqs)

idx100 = Data.index[(Data.Surface == 100) & Data.Substrate.isin(Metals)]
idx110 = Data.index[(Data.Surface == 110) & Data.Substrate.isin(Metals)]
idx111 = Data.index[(Data.Surface == 111) & Data.Substrate.isin(Metals)]
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
plt.figure(1)
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

#plt.xlim([190,550])
#plt.ylim([210,450])
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

plt.figure(2)
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