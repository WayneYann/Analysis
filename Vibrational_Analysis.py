# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 22:11:25 2016

@author: Josh
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import adjustText
import os
mat.rcParams['mathtext.default'] = 'regular'
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/OHx_and_NHx_condensed_data_v04.csv')
Metal_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,0:6]
Data_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,7:]
Data = np.genfromtxt(vibration_file, delimiter=',')[1:,7:]
Metal_Info = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[1:,0:6]

#Calculating reduced mass
MassO = 16
MassH = 1.00794
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

#Plotting Reduced Masses
plt.figure(1)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMO,MrO,'o',markersize=12)
plt.title('Reduced Mass for Oxygen adsorbed on Metals Calculated from DFT',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('$Mr^{0.5}$ $(g/mol)^{0.5}$',size=20, fontweight='bold')
plt.xticks(size=16)
plt.yticks(size=16)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0'))+'}}$')
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO, MrO, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts, vMO,MrO,autoalign=True,va='bottom',
                ha='right',arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#Plotting Reduced Masses
plt.figure(1)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMO,MrH,'or',markersize=12)
plt.title('Reduced Mass for Hydrogen adsorbed on Metals Calculated from DFT',size=16, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('$Mr^{0.5}$ $(g/mol)^{0.5}$',size=20, fontweight='bold')
plt.xticks(size=16)
plt.yticks(size=16)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0'))+'}}$')
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO, MrH, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts, vMO,MrH,autoalign=True,va='bottom',
                ha='right',arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#Plotting Reduced Masses
plt.figure(1)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMO,MrN,'og',markersize=12)
plt.title('Reduced Mass for Nitrogen adsorbed on Metals Calculated from DFT',size=16, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('$Mr^{0.5}$ $(g/mol)^{0.5}$',size=20, fontweight='bold')
plt.xticks(size=16)
plt.yticks(size=16)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0'))+'}}$')
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO, MrN, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts, vMO,MrN,autoalign=True,va='bottom',
                ha='right',arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#Plotting Reduced Masses
plt.figure(1)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMO,MrO,'ob',markersize=12)
plt.plot(vMO,MrN,'og',markersize=12)
plt.title('Reduced Mass for Oxygen and Nitrogen adsorbed on Metals Calculated from DFT',size=16, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('$Mr^{0.5}$ $(g/mol)^{0.5}$',size=20, fontweight='bold')
plt.legend(['Pt-O','Pt-N'],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0'))+'}}$')
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0'))+'}}$')
mat.rc('text', usetex = False)
texts = []

                
Mr2 = np.concatenate((MrO,MrN))
vM2O = np.concatenate((vMO,vMO))

texts = []

for x, y, s in zip(vM2O,Mr2, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True)
plt.show()
#PLotting v(M-O) vs v(M-O2)

idx = np.isfinite(vMO) & np.isfinite(vMO2)
m1,b1 = np.polyfit(vMO[idx], vMO2[idx], 1) 
idx = np.isfinite(vMO) & np.isfinite(vOOatop)
m2,b2 = np.polyfit(vMO[idx], vOOatop[idx], 1) 
idx = np.isfinite(vMO) & np.isfinite(vOObridge)
m3,b3 = np.polyfit(vMO[idx], vOObridge[idx], 1) 

plt.figure(2)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMO,vMO2,'o',markersize=12)
plt.plot(vMO,vOOatop,'o',markersize=12)
plt.plot(vMO,vOObridge,'o',markersize=12)
plt.plot(vMO, m1*vMO+b1,'-k')
plt.plot(vMO, m2*vMO+b2,'-k') 
plt.plot(vMO, m3*vMO+b3,'-k')
#plt.title('Experimental frequencies in adsorbed O2 vs O',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('v(M-O2) and v(O-O) $cm^{-1}$',size=20, fontweight='bold')
plt.legend(['v(M-O2): %sx+%s' %(round(m1,2),round(b1,2)), \
'v(O-O) atop: %sx+%s' %(round(m2,2),round(b2,2)) ,\
'v(O-O) bridge: %sx+%s' %(round(m3,2),round(b3,2))],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)


mat.rc('text', usetex = True)
Marker = []


idxvMO2 = np.isfinite(vMO2)
idxvOOatop = np.isfinite(vOOatop)
idxvOObridge = np.isfinite(vOObridge)
#vMOs = np.concatenate((vMO[idxvMO2],vMO[idxvOOatop],vMO[idxvOObridge]))
#vMO2s = np.concatenate((vMO2[idxvMO2],vOOatop[idxvOOatop],vOObridge[idxvOObridge]))
vMOs = vMO[idxvMO2]
vMO2s = vMO2[idxvMO2]

for i in [i for i, x in enumerate(idxvMO2) if x]:
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,31])).rstrip('.0')+'}}$'))
'''
for i in [i for i, x in enumerate(idxvOOatop) if x]:
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,34])).rstrip('.0')+'}}$'))
for i in [i for i, x in enumerate(idxvOObridge) if x]:
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,37])).rstrip('.0')+'}}$'))
'''
mat.rc('text', usetex = False)

mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMOs, vMO2s, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#Plotting parallel vibrations
vMHhor = []
NN = Data[:,0]
Ms = Data[:,2]
cosa = Data[:,63]
sina = Data[:,64]
ON2overN = Data[:,65]
ONN2cos = Data[:,66]
ONN2sin = Data[:,67]

vMOhor = []
NN = Data[:,0]
Ms = Data[:,2]
cosa = Data[:,3]
sina = Data[:,4]
for i in range(0,len(vMO)):
    if Data[i,0]==3:
        vMOhor.append(vMO[i]*((1+NN[i]*(Ms[i]/MassO)*cosa[i]**2)/((1+NN[i]/2*(Ms[i]/MassO)*sina[i]**2))/2)**0.5)
        #vMOhor.append(vMO[i]*(Mrperp/Mrpara)**(0.5))
    elif Data[i,0] ==4:
        Mr = MassO*NN[i]*Ms[i]/(MassO+NN[i]*Ms[i])
        Mrperp = MassO*NN[i]*cosa[i]*Ms[i]/(MassO+NN[i]*cosa[i]*Ms[i])
        Mrpara = ((Mr**2-Mrperp**2)/2)**(0.5)
        vMOhor.append(vMO[i]*((1+NN[i]*(Ms[i]/MassO)*cosa[i]**2)/((1+NN[i]/2*(Ms[i]/MassO)*sina[i]**2))/2)**0.5)
        #vMOhor.append(vMO[i]*(Mrperp/Mrpara)**(0.5))

    elif Data[i,0]==2:
        vMOhor.append(float('nan'))  
    else:
        vMOhor.append(float('nan'))  
vMOhor = np.array(vMOhor)
idx = np.isfinite(vMO) & np.isfinite(vMOhor)
pHOR = np.polyfit(vMO[idx], vMOhor[idx], 1) 

plt.figure(3)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMO,vMOhor,'o',markersize=12)
plt.plot(vMO, pHOR[0]*vMO+pHOR[1],'-k')
#plt.title('Parallel and Perendicular frequencies for Oxygen',size=12, fontweight='bold')
plt.xlabel('v(M-O) perpendicular $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('v(M-O) parallel $cm^{-1}$',size=20, fontweight='bold')
plt.legend(['v(M-O) parallel: %sx+%s' %(round(pHOR[0],2),round(pHOR[1],1))],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO, vMOhor, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,only_move={'points':'y','text':'y'},autoalign=True, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#Calculating ZPE
c = 29979245800 #cm/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules
T = 298
ZPE = 0.5*c/JeV*h*(vMO +vMOhor*2)

vMOx = np.arange(250,700,1)
ZPEpredicted = 0.5*c/JeV*h*(2*(pHOR[0]*vMOx+pHOR[1])+vMOx)
vMOparax = (pHOR[0]*vMOx+pHOR[1])
pfit = np.polyfit(vMOx, ZPEpredicted, 1)
plt.figure(4)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMOx, ZPEpredicted,'-k',lw=2)
plt.plot(vMO,ZPE,'o',markersize=12)
#plt.title('Zero Point Energy for M-O system',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('ZPE (eV)',size=20, fontweight='bold')
plt.legend(['ZPE predicted: %sx+%s' %(round(pfit[0],5),round(pfit[1],4))],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO, ZPE, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True)
plt.show

#Calculating Entropy
SvRO = h/kB/T*(c*vMO*1/(np.exp(h/kB/T*c*vMO)-1)-np.log(1-np.exp(-1*h*c*vMO/kB/T)) \
+ c*vMOhor*1/(np.exp(h/kB/T*c*vMOhor)-1)-np.log(1-np.exp(-1*h*c*vMOhor/kB/T)) \
+ c*vMOhor*1/(np.exp(h/kB/T*c*vMOhor)-1)-np.log(1-np.exp(-1*h*c*vMOhor/kB/T)) \
)

SvROpredicted = h/kB/T*(c*vMOx*1/(np.exp(h/kB/T*c*vMOx)-1)-np.log(1-np.exp(-1*h*c*vMOx/kB/T)) \
+ c*vMOparax*1/(np.exp(h/kB/T*c*vMOparax)-1)-np.log(1-np.exp(-1*h*c*vMOparax/kB/T)) \
+ c*vMOparax*1/(np.exp(h/kB/T*c*vMOparax)-1)-np.log(1-np.exp(-1*h*c*vMOparax/kB/T)) \
)

SvRO2atop = h/kB/T*(c*vMO2*1/(np.exp(h/kB/T*c*vMO2)-1)-np.log(1-np.exp(-1*h*c*vMO2/kB/T)) \
+ c*vOOatop*1/(np.exp(h/kB/T*c*vOOatop)-1)-np.log(1-np.exp(-1*h*c*vOOatop/kB/T)) \
)

SvRO2bridge = h/kB/T*(c*vMO2*1/(np.exp(h/kB/T*c*vMO2)-1)-np.log(1-np.exp(-1*h*c*vMO2/kB/T)) \
+ c*vOObridge*1/(np.exp(h/kB/T*c*vOObridge)-1)-np.log(1-np.exp(-1*h*c*vOObridge/kB/T)) \
)

idx = np.isfinite(vMO) & np.isfinite(vMO2)
p1 = np.polyfit(vMO[idx], vMO2[idx], 1) 
idx = np.isfinite(vMO) & np.isfinite(vOOatop)
p2 = np.polyfit(vMO[idx], vOOatop[idx], 1) 
idx = np.isfinite(vMO) & np.isfinite(vOObridge)
p3 = np.polyfit(vMO[idx], vOObridge[idx], 1) 
vMO2lin = p1[0]*vMOx+p1[1]
vOOatoplin = p2[0]*vMOx+p2[1]
vOObridgelin = p3[0]*vMOx+p3[1]
SvRO2atoplin = h/kB/T*(c*vMO2lin*1/(np.exp(h/kB/T*c*vMO2lin)-1)-np.log(1-np.exp(-1*h*c*vMO2lin/kB/T)) \
+ c*vOOatoplin*1/(np.exp(h/kB/T*c*vOOatoplin)-1)-np.log(1-np.exp(-1*h*c*vOOatoplin/kB/T)) \
)

SvRO2bridgelin = h/kB/T*(c*vMO2lin*1/(np.exp(h/kB/T*c*vMO2lin)-1)-np.log(1-np.exp(-1*h*c*vMO2lin/kB/T)) \
+ c*vOObridgelin*1/(np.exp(h/kB/T*c*vOObridgelin)-1)-np.log(1-np.exp(-1*h*c*vOObridgelin/kB/T)) \
)

#Plotting vibrations and Oxygen Entropy
T=298
SvRO298 = h/kB/T*(c*vMOx*1/(np.exp(h/kB/T*c*vMOx)-1)-np.log(1-np.exp(-1*h*c*vMOx/kB/T)) \
+ c*vMOparax*1/(np.exp(h/kB/T*c*vMOparax)-1)-np.log(1-np.exp(-1*h*c*vMOparax/kB/T)) \
+ c*vMOparax*1/(np.exp(h/kB/T*c*vMOparax)-1)-np.log(1-np.exp(-1*h*c*vMOparax/kB/T)) \
)
T=200
SvRO200 = h/kB/T*(c*vMOx*1/(np.exp(h/kB/T*c*vMOx)-1)-np.log(1-np.exp(-1*h*c*vMOx/kB/T)) \
+ c*vMOparax*1/(np.exp(h/kB/T*c*vMOparax)-1)-np.log(1-np.exp(-1*h*c*vMOparax/kB/T)) \
+ c*vMOparax*1/(np.exp(h/kB/T*c*vMOparax)-1)-np.log(1-np.exp(-1*h*c*vMOparax/kB/T)) \
)
T=400
SvRO400 = h/kB/T*(c*vMOx*1/(np.exp(h/kB/T*c*vMOx)-1)-np.log(1-np.exp(-1*h*c*vMOx/kB/T)) \
+ c*vMOparax*1/(np.exp(h/kB/T*c*vMOparax)-1)-np.log(1-np.exp(-1*h*c*vMOparax/kB/T)) \
+ c*vMOparax*1/(np.exp(h/kB/T*c*vMOparax)-1)-np.log(1-np.exp(-1*h*c*vMOparax/kB/T)) \
)
plt.figure(6)
plt.figure(figsize=(16,10),dpi=500)

plt.plot(vMOx, SvRO298,'-g',lw=2)
plt.plot(vMOx, SvRO200,'-r',lw=2)
plt.plot(vMOx, SvRO400,'-b',lw=2)
plt.plot(vMO,SvRO,'og',markersize=12)
#plt.title('The Scaling of M-O Entropy',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('Sv/R of M-O',size=20, fontweight='bold')
plt.legend(['Sv/R 298K using linear fits','Sv/R 200K using linear fits','Sv/R 400K using linear fits'],loc=1,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)

Marker = []
for i in range(0,len(vMO)):
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO, SvRO, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show

#Plotting vibrations and M-O2 Entropy

plt.figure(6)
plt.figure(figsize=(16,10),dpi=500)
plt.plot(vMOx, SvRO2atoplin,'-g',lw=2)
plt.plot(vMOx, SvRO2bridgelin,'-b',lw=2)
plt.plot(vMO,SvRO2atop,'og',markersize=12)
plt.plot(vMO,SvRO2bridge,'ob',markersize=12)
#plt.title('The Scaling of M-O Entropy',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('Sv/R of M-O2',size=20, fontweight='bold')
plt.legend(['Atop (using linear v(M-O) fits)','Bridge (using linear v(M-O) fits)'],loc=1,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)

Marker = []
for i in range(0,len(vMO)):
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,34])).rstrip('.0')+'}}$'))
for i in range(0,len(vMO)):
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,37])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)

SvR2O2 = np.concatenate((SvRO2atop,SvRO2bridge))
vM2O = np.concatenate((vMO,vMO))

texts = []
for x, y, s in zip(vM2O,SvR2O2, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True)
plt.show

#A-factor ratio
AO2toOatop = np.exp(2*SvRO-SvRO2atop)
AO2toObridge = np.exp(2*SvRO-SvRO2bridge)
AO2toOatoplin = np.exp(2*SvROpredicted-SvRO2atoplin)
AO2toObridgelin = np.exp(2*SvROpredicted-SvRO2bridgelin)

plt.figure(7)
plt.figure(figsize=(16,10),dpi=500)
plt.plot(vMOx,AO2toOatoplin,'-g',lw=2)
plt.plot(vMOx,AO2toObridgelin,'-b',lw=2)
plt.plot(vMO,AO2toOatop,'og',markersize=12)
plt.plot(vMO,AO2toObridge,'ob',markersize=12)
#plt.title('A-factor ratio (M-O2)/(M-O)',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('A-factor ratio',size=20, fontweight='bold')

plt.legend(['Atop (using linear v(M-O) fits)','Bridge (using linear v(M-O) fits)'],loc=1,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)

Marker = []
for i in range(0,len(vMO)):
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,34])).rstrip('.0')+'}}$'))
for i in range(0,len(vMO)):
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,37])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)

AO2toO = np.concatenate((AO2toOatop,AO2toObridge))
vM2O = np.concatenate((vMO,vMO))

texts = []
for x, y, s in zip(vM2O,AO2toO, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True)
plt.show

#PLotting v(M-O) vs v(M-NHx)
idx = np.isfinite(vMO) & np.isfinite(vMN)
m1,b1 = np.polyfit(vMO[idx], vMN[idx], 1) 
idx = np.isfinite(vMO) & np.isfinite(vMNH)
m2,b2 = np.polyfit(vMO[idx], vMNH[idx], 1) 
idx = np.isfinite(vMO) & np.isfinite(vMNH2)
m3,b3 = np.polyfit(vMO[idx], vMNH2[idx], 1) 

plt.figure(8)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMO,vMN,'o',markersize=12)
plt.plot(vMO,vMNH,'o',markersize=12)
plt.plot(vMO,vMNH2,'o',markersize=12)
plt.plot(vMO, m1*vMO+b1,'-k')
plt.plot(vMO, m2*vMO+b2,'-k') 
plt.plot(vMO, m3*vMO+b3,'-k')
#plt.title('Experimental frequencies in adsorbed O2 vs O',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('v(M-NHx) $cm^{-1}$',size=20, fontweight='bold')
plt.legend(['v(M-N): %sx+%s' %(round(m1,2),round(b1,2)), \
'v(M-NH) atop: %sx+%s' %(round(m2,2),round(b2,2)) ,\
'v(M-NH2) bridge: %sx+%s' %(round(m3,2),round(b3,2))],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)

vMO_2 = np.concatenate((vMO,vMO))
vMO_3 = np.concatenate((vMO_2,vMO))
vMN_2 = np.concatenate((vMN,vMNH))
vMN_3 = np.concatenate((vMN_2,vMNH2))



mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,43])).rstrip('.0')+'}}$'))
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,46])).rstrip('.0')+'}}$'))
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,49])).rstrip('.0')+'}}$'))    

mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO_3, vMN_3, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#PLotting v(M-O) vs v(M-OHx)
idx = np.isfinite(vMO) & np.isfinite(vMOH)
m1,b1 = np.polyfit(vMO[idx], vMOH[idx], 1) 
idx = np.isfinite(vMO) & np.isfinite(vMOH2)
m2,b2 = np.polyfit(vMO[idx], vMOH2[idx], 1) 

plt.figure(9)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMO,vMOH,'o',markersize=12)
plt.plot(vMO,vMOH2,'o',markersize=12)
plt.plot(vMO, m1*vMO+b1,'-k')
plt.plot(vMO, m2*vMO+b2,'-k') 
#plt.title('Experimental frequencies in adsorbed O2 vs O',size=12, fontweight='bold')
plt.xlabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('v(M-OHx) $cm^{-1}$',size=20, fontweight='bold')
plt.legend(['v(M-OH): %sx+%s' %(round(m1,2),round(b1,2)), \
'v(M-OH2): %sx+%s' %(round(m2,2),round(b2,2))],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)

vMO_2 = np.concatenate((vMO,vMO))
vMOHx = np.concatenate((vMOH,vMOH2))



mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,16])).rstrip('.0')+'}}$'))
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,13])).rstrip('.0')+'}}$'))

mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMO_3, vMOHx, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,only_move={'points':'y','text':'y'},autoalign=True, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#PLotting v(M-O) vs v(deltaE)
plt.figure(10)
plt.figure(figsize=(16,10),dpi=500)
EJosh = Data[:,5]
EAbildPedersen = Data[:,7]
idx = np.isfinite(EJosh) & np.isfinite(vMO)
m1,b1 = np.polyfit(EJosh[idx], vMO[idx], 1)
plt.plot(EJosh[idx], m1*EJosh[idx]+b1,'-b')
idx = np.isfinite(EAbildPedersen) & np.isfinite(vMO)
m2,b2 = np.polyfit(EAbildPedersen[idx], vMO[idx], 1) 
plt.plot(EAbildPedersen[idx], m2*EAbildPedersen[idx]+b2,'-g')


plt.plot(EJosh,vMO,'o',markersize=12)
plt.plot(EAbildPedersen,vMO,'o',markersize=12)
#plt.title('Experimental frequencies in adsorbed O2 vs O',size=12, fontweight='bold')
plt.ylabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.xlabel('Eads',size=20, fontweight='bold')
plt.legend(['EJosh: %sx+%s' %(round(m1,2),round(b1,2)), \
'EAbildPedersen$^{1}$: %sx+%s' %(round(m2,2),round(b2,2))],loc=3,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)
vMO_2 = np.concatenate((vMO,vMO))
Eads = np.concatenate((EJosh,EAbildPedersen))



mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))
for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))

mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(Eads, vMO_2, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,only_move={'points':'x','text':'x'},autoalign=True, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#PLotting v(M-O) vs d-band center
idx = np.isfinite(Dband) & np.isfinite(vMO)
m1,b1 = np.polyfit(Dband[idx], vMO[idx], 1) 

plt.figure(12)
plt.figure(figsize=(16,10),dpi=500)
plt.plot(Dband,vMO,'og',markersize=12)
plt.plot(Dband[idx], m1*Dband[idx]+b1,'-g', lw=2)
#plt.title('Experimental frequencies in adsorbed O2 vs O',size=12, fontweight='bold')
plt.xlabel('d-band center (eV)',size=20, fontweight='bold')
plt.ylabel('v(M-O) $cm^{-1}$',size=20, fontweight='bold')
plt.legend(['v(M-O): %sx+%s' %(round(m1,2),round(b1,2))],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)

mat.rc('text', usetex = True)
Marker = []

for i in range(0,len(vMO)):
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+'}}$'))

mat.rc('text', usetex = False)

mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(Dband, vMO, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

#Plotting parallel vibrations for Hydrogen
vMHhor = []
NN = Data[:,0]
Ms = Data[:,2]
cosa = Data[:,63]
sina = Data[:,64]
HN2overN = Data[:,68]
HNN2cos = Data[:,69]
HNN2sin = Data[:,70]
for i in range(0,len(vMHhorlit)):
    if Data[i,0]==3:
        vMHhor.append(vMH[i]*((1+NN[i]*(Ms[i]/MassH)*cosa[i]**2)/((1+NN[i]/2*(Ms[i]/MassH)*sina[i]**2)))**0.5)
         #vMHhor.append(vMH[i]/2**0.5*(sina[i] + HNN2sin[i]/HN2overN[i])/(cosa[i] + HNN2cos[i]/HN2overN[i]))
    elif Data[i,0]==2:
        #vMOhor.append(vMO[i]*((1+NN[i]*(Ms[i]/MassO)*cosa[i]**2)/((1+NN[i]*(Ms[i]/MassO)*sina[i]**2)))**0.5)
        vMHhor.append(float('nan'))   
    else:
        vMHhor.append(float('nan'))  
vMHhor = np.array(vMHhor)
 
idx = np.isfinite(vMH) & np.isfinite(vMHhor)
 
pHOR = np.polyfit(vMH[idx], vMHhor[idx], 1)
plt.figure(3)
plt.figure(figsize=(16,8),dpi=500)
plt.plot(vMH,vMHhor,'o',markersize=12)
plt.plot(vMH, pHOR[0]*vMH +pHOR[1],'-k')
#plt.title('Parallel and Perendicular frequencies for Oxygen',size=12, fontweight='bold')
plt.xlabel('v(M-H) parallel Exp $cm^{-1}$',size=20, fontweight='bold')
plt.ylabel('v(M-H) parallel Central Force $cm^{-1}$',size=20, fontweight='bold')
plt.legend(['v(M-O) parallel: %sx' %(round(pHOR[0],2))],loc=2,prop={'size':20})
plt.xticks(size=16)
plt.yticks(size=16)
mat.rc('text', usetex = True)
Marker = []
for i in range(0,len(vMHhorlit)):
    #Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,10])).rstrip('.0')+','+str(round(Data[i,31])).rstrip('.0')+'}}$'))
    Marker.append(''.join(Metal_Info[i,0]+'${^{'+str(round(Data[i,60])).rstrip('.0')+'}}$'))
mat.rc('text', usetex = False)
texts = []
for x, y, s in zip(vMH, vMHhor, Marker):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=20, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,only_move={'points':'y','text':'y'},autoalign=True, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()





#for i, txt in enumerate(Metal_Info[:,0]):
#    plt.annotate(txt, (vMO[i],Mr[i]), xytext = (-10, 10),
#        textcoords = 'offset points',ha = 'right', va = 'bottom',
#        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'green', alpha = 0.5),
#        arrowprops = dict(arrowstyle = '->'))
        
#plt.get_current_fig_manager().window.raise_()
'''
vMO2 = []
Mr2 = []
for i in range(0,len(vMO)):
    if  ~np.isnan(vMO[i]) and ~np.isnan(Mr[i]):
        vMO2.append(vMO[i])
        Mr2.append(Mr[i])
'''