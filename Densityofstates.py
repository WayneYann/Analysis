#!python2
# -*- coding: utf-8 -*-

"""
Created on Wed Jan 11 14:33:54 2017

@author: lansford
"""
from __future__ import division
import numpy as np
from pykit.densityofstates import dosread
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import os

Atomelement = 'Pt'
Atomindex = 53
Atomelement = 'O'
Atomindex = 0

filepath = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/VASP_Files/Projected_Densities/O')

filepath.replace('\\','/')

if filepath[filepath.__len__()-1] != '/':
    filepath += '/'


# read and return densityofstates object
dos = dosread(filepath + 'DOSCAR', posfile=filepath + 'CONTCAR')
dos.spin
# Get energy grid
energies = dos.energies
fermi = dos.efermi
idx = (np.abs(energies-0)).argmin()
# get_atom_dos will return the DOS for a specific atom label as a numpy array
# Each column corresponds to the s, p, and d orbitals
Atom = Atomelement + str(Atomindex)
d_atom = dos.get_atom_dos(Atom)
d = (dos.get_total_dos())[1]
if dos.spin == True:
    d = np.sum(d,axis=1)
    s = d_atom[:,0] + d_atom[:,1]
    p = d_atom[:,2] + d_atom[:,3]
    dband = d_atom[:,4] + d_atom[:,5]
    d_atom = dband+s+p
if dos.spin == False:
    s = d_atom[:,0]
    p = d_atom[:,1]
    dband = d_atom[:,2]
    d_atom = s + p + dband

s[np.isnan(s)==True]=0
p[np.isnan(p)==True]=0
d[np.isnan(d)==True]=0
dband[np.isnan(dband)==True]=0
energies[np.isnan(energies)==True]=0
print(dos.get_band_energy(band='s',atomtype='C',include='C0'))
print(dos.get_band_filling(band='p',atomtype='',include='C0'))
filling = trapz(d_atom[0:idx],energies[0:idx])
print('filling' + str(filling))
Energy = trapz(dband[0:idx]*energies[0:idx],energies[0:idx])
d_filling = trapz(dband[0:idx],energies[0:idx])
Ed = Energy/d_filling
print('d energy and filling:' + str(Ed) + '/' + str(d_filling))
print('s energy and percent filled')
Energy = trapz(s[0:idx]*energies[0:idx],energies[0:idx])
s_filling = trapz(s[0:idx],energies[0:idx])
Es = Energy/s_filling
Energy = trapz(p[0:idx]*energies[0:idx],energies[0:idx])
p_filling = trapz(p[0:idx],energies[0:idx])
Ep = Energy/p_filling
ps = s_filling/(s_filling+p_filling)
print(str(Es) + '/' + str(ps))
print('p energy and percent filled')
pp = p_filling/(s_filling+p_filling)
print(str(Ep) + '/' + str(pp))
maxs = max(s)
maxp = max(p)
maxd = max(dband)
s[s==0]=np.nan
p[p==0]=np.nan
d[d==0]=np.nan
dband[dband==0]=np.nan
if Atomelement == 'C':
    COorPt = 1000
else:
    COorPt = 1
plt.figure(1)
plt.title('Pt atom that binds to Carbon')
plt.plot(energies,s,'o',zorder=2)
plt.plot(energies,p,'o',zorder=3)
plt.plot(energies,dband,'o',zorder=COorPt)
plt.plot((energies[idx],energies[idx]),[0,max(maxs,maxp,maxd)],'k-',lw=3,zorder=20)
plt.legend(['s','p','d','fermi energy'],loc=2)
plt.ylabel('Projected Density [states/eV]')
plt.xlabel('Energy Level [eV]')
plt.show()
#plt.figure(2)
#plt.plot(energies,d,'o')