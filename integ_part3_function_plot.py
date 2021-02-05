# -*- coding: utf-8 -*-
""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Original code: TrapMethods.py
    
    2019: Extended by Lev Kaplan to include Simpson integration and loop over number of points"""

# integ.py: trapezoid and Simpson integration, a<x<b, N pts, N-1 intervals 
import numpy as np
import matplotlib.pyplot as plt

from pylab import * # import numpy and matplotlib

def func(x):          # function to be integrated
    return sin(exp(5.0)*x+7.0)

A = 0.0
B = 1.0

exact = -cos(exp(5.0)*1.0+7.0)/exp(5.0) - (-cos(exp(5.0)*0.0+7.0)/exp(5.0))

xpoints = np.arange(0.0,1.0,0.001)
fxpoints = [func(x) for x in np.arange(0.0,1.0,0.001)]

plot(xpoints, fxpoints)
title('Function '+r'$ sin(x e^5 + 7)$')
savefig('part3_funcx.png')
show()
