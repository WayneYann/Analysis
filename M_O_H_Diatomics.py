# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 19:46:06 2016

@author: lansford
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import adjustText
import os
from statsmodels.stats import diagnostic
from scipy import stats
from scipy import integrate
mat.rcParams['mathtext.default'] = 'regular'
mat.rcParams['lines.linewidth'] = 7
mat.rcParams['lines.markersize'] = 14
mat.rcParams['legend.numpoints'] = 1
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Diatomic hydrogen and oxygen metals.csv')
Metal_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,0]
Data_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,1:]
Data = np.genfromtxt(vibration_file, delimiter=',')[1:,1:]
Metal_Info = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[1:,0]
DataLabels = []
for i in range(0,len(Data_Labels)):
    DataLabels.append(Data_Labels[i])
    DataLabels.append(i)

percent = 100

"""Getting Frequencies and residuals"""
vMOres =Data[:,2]  - Data[:,3]
vMHres = Data[:,7] - Data[:,8]
idxO = np.isfinite(vMOres)
vMOres = vMOres[idxO]
vMO = Data[:,3][idxO]
idxH = np.isfinite(vMHres)
vMHres = vMHres[idxH]
vMH = Data[:,8][idxH]

"""Getting Lengths and Residuals"""
rMOres = Data[:,4] - Data[:,5]
rMHres = Data[:,9] - Data[:,10]  
idxO = np.isfinite(rMOres)
rMOres = rMOres[idxO]
rMO = Data[:,5][idxO]
idxH = np.isfinite(rMHres)
rMHres = rMHres[idxH]
rMH = Data[:,10][idxH]

"""Getting percentage errors"""
pvMOres = vMOres/vMO*percent
pvMHres = vMHres/vMH*percent
prMOres = rMOres/rMO*percent
prMHres = rMHres/rMH*percent
prMO_Hres = np.concatenate((prMOres,prMHres))
pvMO_Hres = np.concatenate((pvMOres,pvMHres))

"""geting means and sample variances"""
sO = np.std(pvMOres,ddof=1)
sH = np.std(pvMHres,ddof=1)
sO_H = np.std(pvMO_Hres,ddof=1)
muO = np.mean(pvMOres)
muH = np.mean(pvMHres)
muO_H = np.mean(pvMO_Hres)

"""T-test for non-zero mean"""
print('Two-tailed t-test for non-zero means')
print(stats.ttest_1samp(pvMOres, 0))
print(stats.ttest_1samp(pvMHres, 0))
print(stats.ttest_1samp(prMOres, 0))
print(stats.ttest_1samp(prMHres, 0))

"""Test for normality"""
print('Anderson Darling Test for normality')
print(diagnostic.normal_ad(pvMOres))
print(diagnostic.normal_ad(pvMHres))
print(diagnostic.normal_ad(pvMO_Hres))
print(diagnostic.normal_ad(prMOres))
print(diagnostic.normal_ad(prMHres))
print(diagnostic.normal_ad(prMO_Hres))
#print(stats.anderson(pvMO_Hres,dist='norm'))

"""Test for Different Means two-tailed t-test"""
print('Two-tailed T-test for different means')
two_sample_eq_var = stats.ttest_ind(pvMOres, pvMHres, equal_var=True)
print(two_sample_eq_var)
two_sample_diff_var = stats.ttest_ind(pvMOres, pvMHres, equal_var=False)
print(two_sample_diff_var)

"""Testing for different variances p-value"""
print('p-value for F-test of different Standard Deviations')
vO = len(pvMOres)-1
vH = len(pvMHres)-1
F = (sO/sH)**2
print(2*stats.f.cdf(F,vO,vH))

rO = len(prMOres)-1
rH = len(prMHres)-1
F = (np.std(prMOres)/np.std(prMHres))**2
print(2*stats.f.cdf(F,rO,rH))

#Making histograms of data
binnum=5
pvMOhist, pvMOedges = np.histogram(pvMOres,bins=binnum)
pvMHhist, pvMHedges = np.histogram(pvMHres,bins=binnum)
  
#plotting pvMOres data
"""Generating Theortical Distributions"""
binlength = (pvMOedges[1]-pvMOedges[0])
xO = np.arange(pvMOedges[0],pvMOedges[-1],0.0001*percent)
lenx = len(xO)
pvMOTheory = np.zeros(lenx)
for i in range(0,lenx):
    y = np.arange(xO[i]-binlength/2,xO[i]+binlength/2,.0001*percent)
    pvMOf = stats.norm.pdf(y,muO,sO)
    pvMOTheory[i] = integrate.trapz(pvMOf,y)
plt.figure(1)
plt.figure(figsize=(16,10),dpi=500)
hist = np.histogram(pvMOres,bins=binnum)
xticks = pvMOedges[0:binnum]+binlength/2
plt.plot(xO,pvMOTheory*sum(hist[0][:]),'g')
plt.bar(xticks,hist[0][:],width=0.02*percent, align='center')
plt.xticks(xticks)
plt.xlabel('Percent Error in Frequency',size=20, fontweight='bold')
plt.ylabel('Occurances',size=20, fontweight='bold')
plt.title('5 bin Histogram of Diatomic Metal-Oxygen Frequency Data',size=20, fontweight='bold')
plt.xticks(size=16)
plt.yticks(size=16)
plt.legend(['Theoretical Distribution','Occurances'],prop={'size':20})
plt.show()

#plotting pvMHres data
"""Generating Theortical Distributions"""
binlength = (pvMHedges[1]-pvMHedges[0])
xH = np.arange(pvMHedges[0],pvMHedges[-1],0.0001*percent)
lenx = len(xH)
pvMHTheory = np.zeros(lenx)
for i in range(0,lenx):
    y = np.arange(xH[i]-binlength/2,xH[i]+binlength/2,.0001*percent)
    pvMHf = stats.norm.pdf(y,muH,sH)
    pvMHTheory[i] = integrate.trapz(pvMHf,y)
plt.figure(2)
plt.figure(figsize=(16,10),dpi=500)
hist = np.histogram(pvMHres,bins=binnum)
xticks = pvMHedges[0:binnum]+binlength/2
plt.plot(xH,pvMHTheory*sum(hist[0][:]),'g')
plt.bar(xticks,hist[0][:],width=0.02*percent, align='center')
plt.xticks(xticks)
plt.xlabel('Error in Frequency (%)',size=20, fontweight='bold')
plt.ylabel('Occurances',size=20, fontweight='bold')
plt.xticks(size=16)
plt.yticks(size=16)
plt.title('5 bin Histogram of Diatomic Metal-Hydrogen Frequency Data',size=20, fontweight='bold')
plt.legend(['Theoretical Distribution','Occurances'],prop={'size':20})
plt.show()

binnum=9
pvMO_Hhist, pvMO_Hedges = np.histogram(pvMO_Hres,bins=binnum)
#plotting pvMO_Hres data
"""Generating Theortical Distributions"""
binlength = (pvMO_Hedges[1]-pvMO_Hedges[0])
xO_H = np.arange(pvMO_Hedges[0],pvMO_Hedges[-1],0.0001*percent)
lenx = len(xO_H)
pvMO_HTheory = np.zeros(lenx)
for i in range(0,lenx):
    y = np.arange(xO_H[i]-binlength/2,xO_H[i]+binlength/2,.0001*percent)
    pvMO_Hf = stats.norm.pdf(y,muO_H,sO_H)
    pvMO_HTheory[i] = integrate.trapz(pvMO_Hf,y)
plt.figure(3)
plt.figure(figsize=(16,10),dpi=500)
hist = np.histogram(pvMO_Hres,bins=binnum)
xticks = pvMO_Hedges[0:binnum]+binlength/2
plt.plot(xO_H,pvMO_HTheory*sum(hist[0][:]),'g')
plt.bar(xticks,hist[0][:],width=0.02*percent, align='center')
plt.xticks(xticks)
plt.xlabel('Error in Frequency (%)',size=20, fontweight='bold')
plt.ylabel('Occurances',size=20, fontweight='bold')
plt.xticks(size=16)
plt.yticks(size=16)
plt.title('Histogram of Diatomic Metal-Oxygen and Metal-Hydrogen Frequency Data',size=20, fontweight='bold')
plt.legend(['Theoretical Distribution','Occurances'],prop={'size':20})
plt.show()