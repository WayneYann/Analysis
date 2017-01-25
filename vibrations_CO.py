# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:12:15 2016

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
mat.rcParams['lines.linewidth'] = 5
mat.rcParams['lines.markersize'] = 20

vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Energy_Frequency_Slopes_CO.csv')
Data = pd.read_csv(vibration_file,sep=',',header=0)
cols = Data.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
Data.columns = cols


nS=6;nP=8
rA = Data.rA; rB = Data.rB
s_frac_C = Data.s_frac_C;

#s_frac_C = 0;
p_frac_C = 1-s_frac_C; s_frac_CO = Data.s_frac_CO; p_frac_CO = 1-s_frac_CO
#s_frac_CO = s_frac_C; p_frac_CO = p_frac_C
RnA = (s_frac_C*nS*(nS+1)*rA**nP+p_frac_C*nP*(nP+1)*rA**nS)/(rA**2*(s_frac_C*rA**nP+p_frac_C*rA**nS))
RnB = (s_frac_CO*nS*(nS+1)*rB**nP+p_frac_CO*nP*(nP+1)*rB**nS)/(rB**2*(s_frac_CO*rB**nP+p_frac_CO*rB**nS))
mVpred = Data.mE_no_Disp*Data.MA/Data.MB*(RnB/RnA)*(Data.vA_Pt/Data.vB_Pt)
mV = Data.mv
Marker = np.array(Data.Surface)
mVCO = mV[0:3]
mVpredCO= mVpred[0:3]
print(np.mean(abs(mVpredCO-mVCO)))

rB = float(Data.rB[(Data.Surface=='111') & (Data.Adsorbates=='C/CO')])
rA = float(Data.rA[(Data.Surface=='111') & (Data.Adsorbates=='C/CO')])
s_frac_CO = np.linspace(0,1,50); p_frac_CO = 1-s_frac_CO
RnB = (s_frac_CO*nS*(nS+1)*rB**nP+p_frac_CO*nP*(nP+1)*rB**nS)/(rB**2*(s_frac_CO*rB**nP+p_frac_CO*rB**nS))
RnA = (s_frac_CO*nS*(nS+1)*rA**nP+p_frac_CO*nP*(nP+1)*rA**nS)/(rA**2*(s_frac_CO*rA**nP+p_frac_CO*rA**nS))


plt.figure(1)
plt.figure(figsize=(12,10))
plt.plot(s_frac_CO,RnA,'g-',s_frac_CO,RnB,'b--')
plt.xlabel(r'$f_{s}$',size=28)
plt.ylabel(r'$R_{n}$',size=28)
plt.xticks(size=28)
plt.yticks(size=28)
#plt.xlim([0.56,0.7])
plt.ylim([11.5,25])
plt.legend(['C','CO']
,loc=1,prop={'size':24},frameon=False)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
texts = []
