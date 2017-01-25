# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:20:15 2016

@author: lansford
"""

from __future__ import division
import os

import Text_Parse
import shutil

directory=os.path.expanduser('C:\Users\lansf\Documents\VASP_Files\\551K500eV')

directory_v02 = ''.join([directory,'_compressed'])

if not os.path.exists(directory_v02):
    os.makedirs(directory_v02)
full_file_paths = Text_Parse.get_filepaths(directory)
full_file_paths = [i for i in full_file_paths if not ('Ru0001' in i)]
#full_file_paths = [i for i in full_file_paths if not ('Gases' in i)]
CONTCAR_FILES = Text_Parse.list_contains(full_file_paths,"CONTCAR")
OSZICAR_FILES = Text_Parse.list_contains(full_file_paths,"OSZICAR")
OUTCAR_FILES = Text_Parse.list_contains(full_file_paths,"OUTCAR")
CONTCAR_FILES = [i for i in CONTCAR_FILES if not ('freq' in i)]
CONTCAR_FILES = [i for i in CONTCAR_FILES if not ('CONTCAR_old' in i)]
OSZICAR_FILES = [i for i in OSZICAR_FILES if not ('freq' in i)]
OSZICAR_FILES = [i for i in OSZICAR_FILES if not ('OSZICAR_old' in i)]
OUTCAR_FILES = [i for i in OUTCAR_FILES if not ('freq' in i)]
freq_FILES = Text_Parse.list_contains(full_file_paths,"freq\\OUTCAR")
freq_FILES = [i for i in freq_FILES if not ('OUTCAR_015' in i)]
freq_FILES = [i for i in freq_FILES if not ('OUTCAR_025' in i)]
OUTCAR_FILES = [i for i in OUTCAR_FILES if not ('OUTCAR_015' in i)]
OUTCAR_FILES = [i for i in OUTCAR_FILES if not ('OUTCAR_025' in i)]

Freq_Names = Text_Parse.between_values(freq_FILES,"\\551K500eV\\","\\freq")
OUTCAR_Names = Text_Parse.between_values(OUTCAR_FILES,"\\551K500eV\\","\\OUTCAR")
Surface_Names = Text_Parse.between_values(OSZICAR_FILES,"\\551K500eV\\","\\OSZICAR")
for i in range(0,len(OSZICAR_FILES)):
    Surface_Names[i] = Surface_Names[i].replace("\\","_")
    shutil.copy(OSZICAR_FILES[i],directory_v02+"\\"+Surface_Names[i]+"_OSZICAR")
Surface_Names = Text_Parse.between_values(CONTCAR_FILES,"\\551K500eV\\","\\CONTCAR")   
for i in range(0,len(CONTCAR_FILES)):
    Surface_Names[i] = Surface_Names[i].replace("\\","_")
    shutil.copy(CONTCAR_FILES[i],directory_v02+"\\"+Surface_Names[i]+"_CONTCAR")
                    
for i in range(0,len(Freq_Names)):
    Freq_Names[i] = Freq_Names[i].replace("\\","_")
    freq_file2 = directory_v02+"\\"+Freq_Names[i]+"_freq"
    g = open(freq_file2, 'w+')
    with open(freq_FILES[i], 'r') as f:
        for line in f:
            if 'Eigenvectors and eigenvalues of the dynamical matrix' in line:                
                for line in f: # now you are at the lines you want
                    g.write(line)
    g.close()
    
for i in range(0,len(OUTCAR_Names)):
    OUTCAR_Names[i] = OUTCAR_Names[i].replace("\\","_")
    OUTCAR_file2 = directory_v02+"\\"+OUTCAR_Names[i]+"_OUTCAR"
    g = open(OUTCAR_file2, 'w+')
    for line in reversed(open(OUTCAR_FILES[i]).readlines()):
        g.write(line)
        if 'DFTD3 V3.0 Rev 1' in line:                
            break
    g.close()