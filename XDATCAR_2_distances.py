# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:58:40 2016

@author: Josh
"""
from __future__ import division
import os
import numpy as np
import Text_Parse
#os.listdir('C:/Users/Josh/Desktop/XSD files')
# Run the above function and store its results in a variable.
filename=os.path.expanduser('~\Box Sync\Synced_Files\Coding\Research\Analysis\O2_PW91_XDATCAR')
#labels1212 = np.genfromtxt('./Lansford_Stats_HW3_data/12.12.csv', delimiter=',', dtype=str)[0,:]
#data1212 = np.genfromtxt('C:\Users\Josh\Desktop\XSD files', delimiter=',')[1:,:]
Surface_Distances = []

Multiplier = Text_Parse.file2listfixed(filename,2,4,0,None)
Multiplier = Text_Parse.list_split(Multiplier," ")
Multiplier = Text_Parse.string2float(Multiplier)
Multiplier = np.array(Multiplier)

Coordinates = Text_Parse.file2listflex(filename,3,['0'],0,None)
Coordinates = Text_Parse.list_split(Coordinates," ")
Coordinates = Text_Parse.string2float(Coordinates)
Coordinates = np.array(Coordinates)
Coordinates = np.dot(Coordinates,Multiplier)