# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:20:17 2016
　
@author: lansford
　
"""
from __future__ import division
from ase.thermochemistry import IdealGasThermo
from ase import Atoms
OOH = Atoms('OOH')
OOHfreqs = [.484418503,.146088952,.003153291,.002971872]
OOHThermo = IdealGasThermo(OOHfreqs,'linear',atoms=OOH,symmetrynumber=1,spin=0.5,natoms=3)
OOHgibbs = OOHThermo.get_gibbs_energy(298,101325);
OOHZPE = OOHThermo.get_enthalpy(0);
print(OOHgibbs)
print(OOHZPE)