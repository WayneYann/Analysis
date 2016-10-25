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
from scipy import stats
from ase.thermochemistry import IdealGasThermo
from ase import Atoms
from scipy import stats
from statsmodels.stats import diagnostic
import statmech
c = 29979245800 #cm/s
h = 6.6260693*10**(-34) #planks constant
kB = 1.3806505*10**(-23) #boltzman's constant
JeV = 1.60217653*10**(-19) #eV to Joules
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['legend.numpoints'] = 1
mat.rcParams['lines.linewidth'] = 5
mat.rcParams['lines.markersize'] = 20
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Data.csv')
#with open(vibration_file,'rb') as csvfile:
# newfile = csv.reader(csvfile, delimiter=',')
Data = pd.read_csv(vibration_file,sep=',',header=0)
cols = Data.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
Data.columns = cols
vibration_file2 = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/O_H_N_Quals_Data.csv')
Data2 = pd.read_csv(vibration_file2,sep=',',header=0)
#Available Adsorbates: OHx, NHx, CHx, O2, OOH, N2, CO

Ad1 = 'O'; Ad2 = 'OH'; Ad3 = 'CO'; Ad4 = 'OO'; 
Surf1 = 111; NN1 = 3; BadMetals = ['Cr']
freqs = Data.Zfrequency; masses = Data.Adsorbate_mass 
#freqs = Data.frequency_corrected; masses = Data.Mr 

"""Dont change below this line"""
#predicted and actual slope values
v1 = [];E1 = [];L1 = [];v2 = [];E2 = [];L2 = []; M1 = []; M2 = []; G1 = []; G2 = []
v3 = [];E3 = [];L3 = [];v4 = [];E4 = [];L4 = []; M3 = []; M4 = []; G3 = []; G4 = []
MetalLabels=[]
DataPoints = len(Data.frequency_corrected)
for i in range(0,DataPoints):
    if (Data.Surface[i] == Surf1 and Data.Substrate[i] not in BadMetals and Data.NN[i] == NN1):
        if Data.Adsorbate[i] ==Ad1:
            v1.append(freqs[i])
            E1.append(Data.Eads[i])
            L1.append(Data.mindistance[i])
            M1.append(masses[i])
            G1.append(Data.vibGibbs[i])
        elif Data.Adsorbate[i] ==Ad2:
            v2.append(freqs[i])
            E2.append(Data.Eads[i])
            L2.append(Data.mindistance[i])
            M2.append(masses[i])
            G2.append(Data.vibGibbs[i])
        elif Data.Adsorbate[i] ==Ad3:
            MetalLabels.append(Data.Substrate[i])
            v3.append(freqs[i])
            E3.append(Data.Eads[i])
            L3.append(Data.mindistance[i])
            M3.append(masses[i])
            G3.append(Data.vibGibbs[i])
        elif Data.Adsorbate[i] ==Ad4:
            v4.append(freqs[i])
            E4.append(Data.Eads[i])
            L4.append(Data.mindistance[i])
            M4.append(masses[i])
            G4.append(Data.vibGibbs[i])

v1 = np.array(v1); v2=np.array(v2); E1 = np.array(E1); E2 = np.array(E2); L1 = np.array(L1); L2 = np.array(L2); M1 = np.array(M1); M2 = np.array(M2); G1 = np.array(G1); G2 = np.array(G2)
v3 = np.array(v3); v4=np.array(v4); E3 = np.array(E3); E4 = np.array(E4); L3 = np.array(L3); L4 = np.array(L4); M3 = np.array(M3); M4 = np.array(M4); G3 = np.array(G3); G4 = np.array(G4)
#redefining O2 eenergy change from molecular oxygen
E1 = E1+9.5884297/2-1.9555224

mat.rc('text', usetex = False)
idx1 = np.isfinite(E3) & np.isfinite(G3)
idx2 = np.isfinite(E4) & np.isfinite(G4)
p = np.polyfit(E3[idx1],G3[idx1], 1) 
m1 = p[0]
b1=p[1]
p = np.polyfit(E4[idx2],G4[idx2], 1) 
m2 = p[0]
b2=p[1]
idx3 = np.isfinite(E1) & np.isfinite(G1)
p = np.polyfit(E1[idx3],G1[idx3], 1) 
m3 = p[0]
b3=p[1]
plt.figure(1)
plt.figure(figsize=(14,10))
plt.plot(E1[idx3], G1[idx3],'rs')
plt.plot(E4, G4,'bo')
plt.plot(E3, G3,'g^')
plt.plot(np.array([np.min(E1[idx3]),np.max(E1[idx3])]), m3*np.array([np.min(E1[idx3]),np.max(E1[idx3])])+b3,'-r')
plt.plot(np.array([np.min(E4[idx2]),np.max(E4[idx2])]), m2*np.array([np.min(E4[idx2]),np.max(E4[idx2])])+b2,'-b')
plt.plot(np.array([np.min(E3[idx1]),np.max(E3[idx1])]), m1*np.array([np.min(E3[idx1]),np.max(E3[idx1])])+b1,'-g')
plt.xlabel('${\Delta}$E$_{ads}$ [eV]',size=32)
plt.ylabel('$G_{vib}$ at 298K [eV]',size=32)
plt.legend(['O: $G_{vib}$=%4.3f${\Delta}$E$_{ads}$ + %4.3f eV' %(m3,b3),'O$_{2}$: $G_{vib}$=%4.3f${\Delta}$E$_{ads}$ + %4.3f eV' %(m2,b2), '%s: $G_{vib}$=%4.3f${\Delta}$E$_{ads}$ + %1.0f eV' %(Ad3,m1,abs(b1))],loc=3,prop={'size':28},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.ylim([-.03,0.15])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
M2 = np.concatenate((MetalLabels,MetalLabels,np.array(MetalLabels)[idx3]))
E34 = np.concatenate((E3,E4,E1[idx3]))
G34 = np.concatenate((G3,G4,G1[idx3]))
for x, y, s in zip(E34, G34, M2):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
                       arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

p = np.polyfit(v1,v2, 1) 
m1 = p[0]
b1=p[1]
p = np.polyfit(v1,v4, 1) 
m2 = p[0]
b2=p[1]


plt.figure(2)
plt.figure(figsize=(14,10))
plt.plot(v1, v2,'g^')
plt.plot(v1, v4,'bo')
plt.plot(np.array([np.min(v1),np.max(v1)]), m1*np.array([np.min(v1),np.max(v1)])+b1,'-g')
plt.plot(np.array([np.min(v1),np.max(v1)]), m2*np.array([np.min(v1),np.max(v1)])+b2,'-b')
plt.xlabel(r'$\mathbf{\nu}_{\perp (M-O)}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp (M-Ox)}$ [$cm^{-1}$]',size=32)
plt.legend([r'$\mathbf{\nu}_{\perp,OH}$=%s$\mathbf{\nu}_{\perp,O}$ + %s $cm^{-1}$' %(round(m1,2),(round(b1,1))),r'$\mathbf{\nu}_{\perp,O_{2}}$=%s$\mathbf{\nu}_{\perp,O}$ - %s $cm^{-1}$' %(round(m2,2),abs(round(b2,1)))],loc=2,prop={'size':28},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.ylim([90,400])
plt.xlim([295,500])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
M2 = np.concatenate((MetalLabels,MetalLabels))
v11 = np.concatenate((v1,v1))
v24 = np.concatenate((v2,v4))
for x, y, s in zip(v11, v24, M2):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=32, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, )
plt.show()

idx1 = np.isfinite(Data2['v(M-O)']) & np.isfinite(Data2['v(M-O2)'])
p = np.polyfit(Data2['v(M-O)'][idx1],Data2['v(M-O2)'][idx1], 1) 
m3 = p[0]
b3=p[1]
idx2 = np.isfinite(Data2['v(M-O)']) & np.isfinite(Data2['v(M-N)'])
p = np.polyfit(Data2['v(M-O)'][idx2],Data2['v(M-N)'][idx2], 1) 
m4 = p[0]
b4=p[1]

plt.figure(3)
plt.figure(figsize=(14,10))
plt.plot(Data2['v(M-O)'], Data2['v(M-O2)'],'g^')
plt.plot(Data2['v(M-O)'], Data2['v(M-N)'],'bo')
plt.plot(np.array([np.min(Data2['v(M-O)'][idx1]),np.max(Data2['v(M-O)'][idx1])]), m3*np.array([np.min(Data2['v(M-O)'][idx1]),np.max(Data2['v(M-O)'][idx1])])+b3,'-g')
plt.plot(np.array([np.min(Data2['v(M-O)'][idx2]),np.max(Data2['v(M-O)'][idx2])]), m4*np.array([np.min(Data2['v(M-O)'][idx2]),np.max(Data2['v(M-O)'][idx2])])+b4,'-b')
plt.xlabel(r'$\mathbf{\nu}_{\perp (M-O)}$ [$cm^{-1}$]',size=32)
plt.ylabel(r'$\mathbf{\nu}_{\perp (M-X)}$ [$cm^{-1}$]',size=32)
plt.legend([r'$\mathbf{\nu}_{\perp,O_{2}}$=%s$\mathbf{\nu}_{\perp,O}$ + %s eV' %(round(m3,2),(round(b3,1))),r'$\mathbf{\nu}_{\perp,N}$=%s$\mathbf{\nu}_{\perp,O}$ + %s eV' %(round(m4,2),round(b4,1))],loc=2,prop={'size':28},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
#plt.xlim([245,700])
plt.ylim([190,600])
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
M3 = np.concatenate((Data2['Surface'][idx1],Data2.Surface[idx2]))
vOO = np.concatenate((Data2['v(M-O)'][idx1],Data2['v(M-O)'][idx2]))
vO2N = np.concatenate((Data2['v(M-O2)'][idx1],Data2['v(M-N)'][idx2]))
for x, y, s in zip(vOO, vO2N, M3):
    texts.append(plt.text(x, y, s, bbox={'pad':0, 'alpha':0}, size=22, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'}, 
        arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

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
delGCOPt = G3[PtI]+E3[PtI]-COgibbs

delGCOAg = G3[AgI]+E3[AgI]-COgibbs
delGCOPtAg = G3[PtI]+E3[AgI]-COgibbs
delGOPt = G1[PtI]+E1[PtI]-O2gibbs/2
delGOAg = G1[AgI]+E1[AgI]-O2gibbs/2
delGOPtAg = G1[PtI]+E1[AgI]-O2gibbs/2


KcoPt = np.exp(-delGCOPt/(kB*298/JeV))
KcoAg = np.exp(-delGCOAg/(kB*298/JeV))
KcoPtAg = np.exp(-delGCOPtAg/(kB*298/JeV))
KoPt = np.exp(-delGOPt/(kB*298/JeV))
KoAg = np.exp(-delGOAg/(kB*298/JeV))
KoPtAg = np.exp(-delGOPtAg/(kB*298/JeV))

import pylab
idx1 = np.isfinite(E1) & np.isfinite(G1)
idx2 = np.isfinite(E3) & np.isfinite(G3)
p = np.polyfit(E1[idx1],G1[idx1], 1) 
mO = p[0]
bO=p[1]
p = np.polyfit(E3[idx2],G3[idx2], 1) 
mCO = p[0]
bCO=p[1]

pCO = 0.5

xxmin = -0.7
xxmax = 0
yymin = -.7
yymax = -0.3
xx = pylab.linspace(xxmin, xxmax, 100)
yy = pylab.linspace(yymin, yymax, 100)
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
pylab.figure(4)
pylab.pcolor(xx, yy, zz)
#pylab.title('Coverage of CO with Scaling of Gvib')
pylab.ylabel('CO Eads [eV]')
pylab.xlabel('O Eads [eV]')
pylab.colorbar()
pylab.ylim([yymin,yymax])
pylab.xlim([xxmin,xxmax])
pylab.show()

def thetaCO(EadsO, EadsCO):
    delGCO = EadsCO-COgibbs +mCO*EadsCO+bCO+G3[PtI]
    Kco = np.exp(-1*delGCO/(kB*298/JeV))
    delGO = EadsO-O2gibbs/2 +G1[PtI]
    Ko = np.exp(-1*delGO/(kB*298/JeV))
    theta = Kco*pCO/(1+Kco*pCO+(Ko*(1-pCO))**0.5)
    return theta

for i in xrange(len(xx)):
    for j in xrange(len(yy)):
        zz[j, i] = thetaCO(xx[i], yy[j])
pylab.figure(5)
pylab.pcolor(xx, yy, zz)
pylab.ylim([yymin,yymax])
pylab.xlim([xxmin,xxmax])
#pylab.title('Coverage of CO without Vibrational Scaling')
pylab.ylabel('CO Eads [eV]')
pylab.xlabel('O Eads [eV]')
pylab.colorbar()
pylab.show()



"""Plot of Gvib/delG"""
mat.rc('text', usetex = False)
idx1 = np.isfinite(E1) & np.isfinite(G1)
idx3 = np.isfinite(E3) & np.isfinite(G3)
p = np.polyfit(E3[idx3],G3[idx3], 1) 
m3 = p[0]
b3=p[1]
p = np.polyfit(E1[idx1],G1[idx1], 1) 
m1 = p[0]
b1=p[1]

EO = np.linspace(-2,2,2000)
ECO = EO
delGO = EO*(1+m1)+b1-O2gibbs/2
delGCO = ECO*(1+m3)+b3-COgibbs
GvibO = EO*m1+b1
GvibCO = ECO*m3+b3
for i in range(0,len(EO)):
    if delGO[i] > -0.1 and delGO[i] < 0.1:
        delGO[i] = 'nan'
        GvibO[i] = 'nan'
    if delGCO[i] > -0.1 and delGCO[i] < 0.1:
        delGCO[i] = 'nan'
        GvibCO[i] = 'nan'


plt.figure(10)
plt.figure(figsize=(14,10))
plt.plot(delGO,(GvibO/delGO),'b')
plt.plot(delGCO,(GvibCO/delGCO),'g')
plt.xlim([-1,1])
plt.ylim([-0.5,0.5])
plt.xlabel('${\Delta}$G$_{ads}$ [eV]',size=32)
plt.ylabel('G$_{vib}$/${\Delta}$G$_{ads}$',size=32)
plt.legend(['O','CO'],loc=2,prop={'size':28},frameon=False)
plt.xticks(size=28)
plt.yticks(size=28)
plt.show()

#T test
mvp = np.array([-0.038736455,-0.136159028,-0.186696644,-0.069105119,-0.225828899,-0.377797055,-0.151966429,0.138888407,0.177621218,-0.149457563,-0.164333338,-0.069042195,])
mE= np.array([-0.2,-0.274,-0.12,-0.22,-0.33,-0.57,-0.34,0.13,0.05,-0.26,-0.03,-0.12,])


"""geting means and sample variances"""
svp = np.std(mvp,ddof=1)
sE = np.std(mE,ddof=1)
meanvp = np.mean(mvp)
meanE = np.mean(mE)

"""Test for normality"""
print('Anderson Darling Test for normality')
print(diagnostic.normal_ad(mvp))
print(diagnostic.normal_ad(mE))

print(stats.ttest_1samp(mvp, 0))
print(stats.ttest_1samp(mE, 0))
#print(stats.anderson(pvMO_Hres,dist='norm'))

"""Test for Different Means two-tailed t-test"""
print('Two-tailed T-test for different means')
two_sample_diff_var = stats.ttest_ind(mvp, mE, equal_var=False)
print(two_sample_diff_var)

Pairs = ['O/OH','O/OH$_{2}$','O/O$_{2}$','O/OOH','N/NH','N/NH$_{2}$','NH/NH$_{2}$','C/CH','C/CH$_{2}$','C/CH$_{3}$','C/CO','O/CH$_{3}$']
me = [0.43,0.096,0.53,0.74,0.65,0.29,0.47,0.83,0.47,0.24,0.57,0.66]
mv = [0.63,0.37,0.65,0.96,0.98,0.86,0.81,0.7,0.42,0.5,0.6,0.78]
mvp = [0.59,0.23,0.46,0.89,0.75,0.48,0.66,0.84,0.60,0.35,0.44,0.71]
plt.figure(7)
plt.figure(figsize=(12,10))
p = np.polyfit(E1, E2, 1)
m1 = p[0]
b1=p[1]
R2 = 1 - sum(((m1*E1+b1)-E2)**2)/sum((E2-np.mean(E2))**2)
#plt.figure(figsize=(12,10))
plt.plot(mv, me,'bo')
plt.plot(mv, mvp,'go')
plt.plot(mv,mv,'k-')
#plt.title('Comparing The Predicted Frequency and Energy Slopes',size=20)
plt.xlabel('Actual Frequency Slope (m$_{V}$)',size=24)
plt.ylabel('m$_{E}$ and m$_{P}$',size=24)
plt.legend(['m$_{E}$','m$_{P}$','Parity'],loc=2,prop={'size':20})
plt.xticks(size=20)
plt.yticks(size=20)

texts = []
for a, b, s in zip(mv,mv,Pairs):
    texts.append(plt.text(a, b, s, bbox={'pad':0, 'alpha':0}, size=24, fontweight='bold',style='normal',name ='Calibri'))
adjustText.adjust_text(texts,autoalign=True,only_move={'points':'y','text':'y'},
    arrowprops=dict(arrowstyle="-", color='k', lw=2))
plt.show()

vCu = np.array([1802.14,254.26,220.5,218.67,107.87,104.36])
vPt = np.array([1731.9,349.97,322.15,319.11,161.45,160.01])
vAg = np.array([1887.8,184.08,136.43,135.39,28.845,17.72])
vRh = np.array([1730.3,340.98,267.67,266.45,146.44,145.42])
vAl = np.array([1591.86,193.18,169.06,168.5,71.08,70.13])
vIr = np.array([1704.79,358.09,255.558,248.724,150.6,142.86])
vNi = np.array([1751.16,348.58,276.47,273.86,138.5,134.67])
vPd = np.array([1750.07,351.6,348.65,323.906,164.009,162.306])


Pairs = ['Cu','Pt','Ag','Rh','Al','Ir','Ni','Pd']
vibgibbs = [statmech.vibgibbs(vCu,298),statmech.vibgibbs(vPt,298),statmech.vibgibbs(vAg,298),statmech.vibgibbs(vRh,298),statmech.vibgibbs(vAl,298),statmech.vibgibbs(vIr,298),statmech.vibgibbs(vNi,298),statmech.vibgibbs(vPd,298)]
vibgibbsTaylor = [statmech.vibgibbsTaylor(vCu,298,vCu),statmech.vibgibbsTaylor(vPt,298,vCu),statmech.vibgibbsTaylor(vAg,298,vCu),statmech.vibgibbsTaylor(vRh,298,vCu),statmech.vibgibbsTaylor(vAl,298,vCu),statmech.vibgibbsTaylor(vIr,298,vCu),statmech.vibgibbsTaylor(vNi,298,vCu),statmech.vibgibbsTaylor(vPd,298,vCu)]
plt.figure(8)
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