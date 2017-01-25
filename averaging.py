# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 14:15:58 2016

@author: lansford
"""
from __future__ import division
a = 0.3; b = 1-a;
c = 6*a; d = 8*b;

x = 2; y = 3; w = 4; z = 3;

val = (a*x+b*y)/(c*x+d*y)*(c*w+d*z)/(a*w+b*z)
print(val)