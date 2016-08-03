# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 01:14:46 2016

@author: Josh
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import adjustText
import os
from scipy.interpolate import interp1d
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['lines.linewidth'] = 3
mat.rcParams['lines.markersize'] = 20
mat.rcParams['legend.numpoints'] = 1
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/OHx_and_NHx_condensed_data_v05.csv')
Metal_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,0:6]
Data_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,7:]
Data = np.genfromtxt(vibration_file, delimiter=',')[1:,7:]
Metal_Info = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[1:,0:6]

def vibentropy(vibrations,temperature):
   """"Input: vibrations in 1/cm and temperature in Kelvin as list
   Returns entropy in eV/K"""
   T = temperature
   vib = vibrations
   c = 29979245800 #cm/s
   h = 6.6260693*10**(-34) #planks constant
   kB = 1.3806505*10**(-23) #boltzman's constant
   JeV = 1.60217653*10**(-19) #eV to Joules
   n = len(vibrations)
   entropy = 0
   for i in range(0,n):
       entropy = (h/T*c*vib[i]/(np.exp(h/kB/T*c*vib[i])-1)-kB*np.log(1-np.exp(-1*h*c*vib[i]/kB/T)))/JeV + entropy
   return entropy

def vibenergy(vibrations,temperature):
   """"Input: vibrations in 1/cm and temperature in Kelvin as list
   Returns energy in eV"""
   T = temperature
   vib = vibrations
   c = 29979245800 #cm/s
   h = 6.6260693*10**(-34) #planks constant
   kB = 1.3806505*10**(-23) #boltzman's constant
   JeV = 1.60217653*10**(-19) #eV to Joules
   n = len(vibrations)
   energy = 0
   for i in range(0,n):
       energy = h*c*(0.5*vib[i] + vib[i]/(np.exp(h/kB/T*c*vib[i])-1))/JeV + energy
   return energy


#Calculating reduced mass
MassH = 1.00794
MassN = 14*MassH
MassO = 16*MassH
vMO = Data[:,8]
vMO2 = Data[:,29]
vOOatop = Data[:,32]
vOObridge = Data[:,35]
vMN = Data[:,41]
vMNH = Data[:,44]
vMNH2 = Data[:,47]
vMOH = Data[:,14]
vMOH2 = Data[:,11]
Dband = Data[:,53]
DataLabels = []
vMH = Data[:,54]
vMHhorlit = Data[:,55]
for i in range(0,len(Data_Labels)):
    DataLabels.append(Data_Labels[i])
    DataLabels.append(i)
MrO = (1/MassO+1/(Data[:,0]*Data[:,2]*Data[:,3]**2))**(-0.5)
MrH = (1/MassH+1/(Data[:,0]*Data[:,2]*Data[:,63]**2))**(-0.5)
MrN = (1/14+1/(Data[:,0]*Data[:,2]*Data[:,61]**2))**(-0.5)

Eads = Data[:,5]
idx = np.isfinite(Eads) & np.isfinite(vMO)
mEads,bEads = np.polyfit(Eads[idx], vMO[idx], 1)

c = 29979245800 #cm/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules

vMOpara = []
NN = Data[:,0]
Ms = Data[:,2]
cosa = Data[:,3]
sina = Data[:,4]
for i in range(0,len(vMO)):
    if Data[i,0]==3:
        vMOpara.append(vMO[i]*((1+NN[i]*(Ms[i]/MassO)*cosa[i]**2)/((1+NN[i]/2*(Ms[i]/MassO)*sina[i]**2))/2)**0.5)
        #vMOpara.append(vMO[i]*(Mrperp/Mrpara)**(0.5))
    elif Data[i,0] ==4:
        Mr = MassO*NN[i]*Ms[i]/(MassO+NN[i]*Ms[i])
        Mrperp = MassO*NN[i]*cosa[i]*Ms[i]/(MassO+NN[i]*cosa[i]*Ms[i])
        Mrpara = ((Mr**2-Mrperp**2)/2)**(0.5)
        vMOpara.append(vMO[i]*((1+NN[i]*(Ms[i]/MassO)*cosa[i]**2)/((1+NN[i]/2*(Ms[i]/MassO)*sina[i]**2))/2)**0.5)
        #vMOpara.append(vMO[i]*(Mrperp/Mrpara)**(0.5))

    elif Data[i,0]==2:
        vMOpara.append(float('nan'))  
    else:
        vMOpara.append(float('nan'))  
vMOpara = np.array(vMOpara)
idxO = np.isfinite(vMO) & np.isfinite(vMOpara)
pHOR = np.polyfit(vMO[idxO], vMOpara[idxO], 1) 

vMOx = np.arange(250,700,1)
vMOparax = (pHOR[0]*vMOx+pHOR[1])

'''parallel frequencies for M-N'''
vMNpara = []
NN = Data[:,0]
Ms = Data[:,2]
cosa = Data[:,61]
sina = Data[:,62]
for i in range(0,len(vMN)):
    if Data[i,0]==3:
        vMNpara.append(vMN[i]*((1+NN[i]*(Ms[i]/MassN)*cosa[i]**2)/((1+NN[i]/2*(Ms[i]/MassN)*sina[i]**2))/2)**0.5)
        #vMOpara.append(vMO[i]*(Mrperp/Mrpara)**(0.5))
    elif Data[i,0] ==4:
        Mr = MassN*NN[i]*Ms[i]/(MassN+NN[i]*Ms[i])
        Mrperp = MassN*NN[i]*cosa[i]*Ms[i]/(MassN+NN[i]*cosa[i]*Ms[i])
        Mrpara = ((Mr**2-Mrperp**2)/2)**(0.5)
        vMNpara.append(vMN[i]*((1+NN[i]*(Ms[i]/MassN)*cosa[i]**2)/((1+NN[i]/2*(Ms[i]/MassN)*sina[i]**2))/2)**0.5)
        #vMOpara.append(vMO[i]*(Mrperp/Mrpara)**(0.5))

    elif Data[i,0]==2:
        vMNpara.append(float('nan'))  
    else:
        vMNpara.append(float('nan'))  
vMNpara = np.array(vMNpara)
idxN = np.isfinite(vMN) & np.isfinite(vMNpara)
pHOR = np.polyfit(vMN[idxN], vMNpara[idxN], 1) 

vMNx = np.arange(250,700,1)
vMNparax = (pHOR[0]*vMNx+pHOR[1])

'''parallel frequencies for M-H'''

vMHpara = Data[:,55]
idxH = np.isfinite(vMH) & np.isfinite(vMHpara)
pHOR = np.polyfit(vMH[idxH], vMHpara[idxH], 1) 

vMHx = np.arange(800,1300,1)
vMHparax = (pHOR[0]*vMHx+pHOR[1])
GibbsvibO = vibenergy([vMO,vMOpara,vMOpara],298)-298*vibentropy([vMO,vMOpara,vMOpara],298)
GibbsvibN = vibenergy([vMN,vMNpara,vMNpara],298)-298*vibentropy([vMN,vMNpara,vMNpara],298)
GibbsvibH = vibenergy([vMH,vMHpara,vMHpara],298)-298*vibentropy([vMH,vMHpara,vMHpara],298)

GibbsvibOfit = vibenergy([vMOx,vMOparax,vMOparax],298)-298*vibentropy([vMOx,vMOparax,vMOparax],298)
GibbsvibNfit = vibenergy([vMNx,vMNparax,vMNparax],298)-298*vibentropy([vMNx,vMNparax,vMNparax],298)
GibbsvibHfit = vibenergy([vMHx,vMHparax,vMHparax],298)-298*vibentropy([vMHx,vMHparax,vMHparax],298)

"""Getting R2 values"""
fGO = interp1d(vMOx,GibbsvibOfit)
fGN = interp1d(vMNx,GibbsvibNfit)
fGH = interp1d(vMHx,GibbsvibHfit)
R2O = 1 - sum((fGO(vMO[idxO])-GibbsvibO[idxO])**2)/sum((fGO(vMO[idxO])-np.mean(GibbsvibO[idxO]))**2)
R2N = 1 - sum((fGN(vMN[idxN])-GibbsvibN[idxN])**2)/sum((fGN(vMN[idxN])-np.mean(GibbsvibN[idxN]))**2)
R2H = 1 - sum((fGH(vMH[idxH])-GibbsvibH[idxH])**2)/sum((fGH(vMH[idxH])-np.mean(GibbsvibH[idxH]))**2)
"""Plotting Combined Gibbs Energy (Contribution from vibrations)"""

plt.figure(1)
#plt.figure(figsize=(18,10),dpi=500)
plt.figure(figsize=(14,8))
plt.plot(vMOx,GibbsvibOfit,'-b')
plt.plot(vMNx,GibbsvibNfit,'-g')
plt.plot(vMHx,GibbsvibHfit,'-r')


plt.plot(vMO,GibbsvibO,'ob')
plt.plot(vMN,GibbsvibN,'og')
plt.plot(vMH,GibbsvibH,'or')

#plt.title('The Scaling of M-O Entropy',size=12, fontweight='bold')
plt.xlabel('Perependicular Frequency [$cm^{-1}$]',size=24, fontweight='bold')
plt.ylabel('G$_{vib}$ [eV]',size=24, fontweight='bold')
#plt.title('Atomic Oxygen',size=20, fontweight='bold')
plt.legend(['O: R$^{2}$ = %.2f' %(R2O),'N: R$^{2}$ = %.2f' %(R2N),'H: R$^{2}$ = %.2f' %(R2H)],loc=2,prop={'size':24})
plt.xticks(size=20)
plt.yticks(size=20)

Marker = []
for i in range(0,len(vMO)):
    Marker.append(Metal_Info[i,0])
for i in range(0,len(vMN)):
    Marker.append(Metal_Info[i,0])
for i in range(0,len(vMH)):
    Marker.append(Metal_Info[i,0])
mat.rc('text', usetex = False)
texts = []
vM = np.concatenate((vMO,vMN,vMH))
Gibbsvib = np.concatenate((GibbsvibO,GibbsvibN,GibbsvibH))
idx = np.isfinite(vM) & np.isfinite(Gibbsvib)
vM = vM[idx]
Gibbsvib = Gibbsvib[idx]
Marker = np.array(Marker)[idx]
for x, y, s in zip(vM, Gibbsvib, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

"""Plotting Gibbs Energy oxygen (Contribution from vibrations)"""
plt.figure(1)
plt.figure(figsize=(16,10),dpi=500)
plt.plot(vMOx,vibenergy([vMOx,vMOparax,vMOparax],200)-200*vibentropy([vMOx,vMOparax,vMOparax],200),'-g')
plt.plot(vMOx,vibenergy([vMOx,vMOparax,vMOparax],298)-298*vibentropy([vMOx,vMOparax,vMOparax],298),'-r')
plt.plot(vMOx,vibenergy([vMOx,vMOparax,vMOparax],400)-400*vibentropy([vMOx,vMOparax,vMOparax],400),'-b')

plt.plot(vMO,vibenergy([vMO,vMOpara,vMOpara],200)-200*vibentropy([vMO,vMOpara,vMOpara],200),'og')
plt.plot(vMO,vibenergy([vMO,vMOpara,vMOpara],298)-298*vibentropy([vMO,vMOpara,vMOpara],298),'or')
plt.plot(vMO,vibenergy([vMO,vMOpara,vMOpara],400)-400*vibentropy([vMO,vMOpara,vMOpara],400),'ob')

#plt.title('The Scaling of M-O Entropy',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('${\Delta}$G$_{vib}$ (eV)',size=20, fontweight='bold')
plt.title('Atomic Oxygen',size=20, fontweight='bold')
plt.legend(['200 K','298 K','400 K'],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)
Gibbsvib = vibenergy([vMO,vMOpara,vMOpara],298)-298*vibentropy([vMO,vMOpara,vMOpara],298)
Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO, Gibbsvib, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

"""Plotting Gibbs Energy Nitrogen (Contribution from vibration) """
plt.figure(2)
plt.figure(figsize=(16,10),dpi=500)
plt.plot(vMNx,vibenergy([vMNx,vMNparax,vMNparax],200)-200*vibentropy([vMNx,vMNparax,vMNparax],200),'-g')
plt.plot(vMNx,vibenergy([vMNx,vMNparax,vMNparax],298)-298*vibentropy([vMNx,vMNparax,vMNparax],298),'-r')
plt.plot(vMNx,vibenergy([vMNx,vMNparax,vMNparax],400)-400*vibentropy([vMNx,vMNparax,vMNparax],400),'-b')

plt.plot(vMN,vibenergy([vMN,vMNpara,vMNpara],200)-200*vibentropy([vMN,vMNpara,vMNpara],200),'og')
plt.plot(vMN,vibenergy([vMN,vMNpara,vMNpara],298)-298*vibentropy([vMN,vMNpara,vMNpara],298),'or')
plt.plot(vMN,vibenergy([vMN,vMNpara,vMNpara],400)-400*vibentropy([vMN,vMNpara,vMNpara],400),'ob')

#plt.title('The Scaling of M-O Entropy',size=12, fontweight='bold')
plt.xlabel('v(M-N) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('${\Delta}$G$_{vib}$ (eV)',size=20, fontweight='bold')
plt.title('Atomic Nitrogen',size=20, fontweight='bold')
plt.legend(['200 K','298 K','400 K'],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)
Gibbsvib = vibenergy([vMN,vMNpara,vMNpara],298)-298*vibentropy([vMN,vMNpara,vMNpara],298)
Marker = []
for i in range(0,len(vMN)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMN, Gibbsvib, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

"""Plotting Gibbs Energy Hydrogen (Contribution from vibration) """
plt.figure(3)
plt.figure(figsize=(16,10),dpi=500)
plt.plot(vMHx,vibenergy([vMHx,vMHparax,vMHparax],200)-200*vibentropy([vMHx,vMHparax,vMHparax],200),'-g')
plt.plot(vMHx,vibenergy([vMHx,vMHparax,vMHparax],298)-298*vibentropy([vMHx,vMHparax,vMHparax],298),'-r')
plt.plot(vMHx,vibenergy([vMHx,vMHparax,vMHparax],400)-400*vibentropy([vMHx,vMHparax,vMHparax],400),'-b')

plt.plot(vMH,vibenergy([vMH,vMHpara,vMHpara],200)-200*vibentropy([vMH,vMHpara,vMHpara],200),'og')
plt.plot(vMH,vibenergy([vMH,vMHpara,vMHpara],298)-298*vibentropy([vMH,vMHpara,vMHpara],298),'or')
plt.plot(vMH,vibenergy([vMH,vMHpara,vMHpara],400)-400*vibentropy([vMH,vMHpara,vMHpara],400),'ob')

#plt.title('The Scaling of M-O Entropy',size=12, fontweight='bold')
plt.xlabel('v(M-H) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('${\Delta}$G$_{vib}$ (eV)',size=20, fontweight='bold')
plt.title('Atomic Hydrogen',size=20, fontweight='bold')
plt.legend(['200 K','298 K','400 K'],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)
Gibbsvib = vibenergy([vMH,vMHpara,vMHpara],298)-298*vibentropy([vMH,vMHpara,vMHpara],298)
Marker = []
for i in range(0,len(vMH)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMH, Gibbsvib, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()