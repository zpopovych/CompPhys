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
    return cos(5*x + x*x)/x**(3/2)

def funcy(y):  # function to be integrated
    return cos(y**(-4)+5*y**(-2))


def trapezoid(A,B,N):   #integrate from A to B using N points
    h = (B - A)/(N - 1)                     # step size 
    sum = (func(A)+func(B))/2               # (1st + last)/2
    for i in range(1, N-1):        # i goes from 1 to (N-1)-1
       sum += func(A+i*h)
       #sum=float32(sum)            # to simulate single-precision (32 bit) calculation
    return h*sum  

def simpson(A,B,N):
    if ((N-1)%2==1):     #  if number of intervals odd
        print("Simpson's rule requires even number of intervals")
        return 0
    h = (B - A)/(N - 1)                     # step size 
    sum = (func(A)+func(B))/3               # (1st + last)/3
    for i in range(1, N-1,2):        # i loops over odd integers  from 1 to (N-1)-1
       sum += 4/3*func(A+i*h)
       #sum=float32(sum)
    for i in range(2, N-1,2):        # i loops over even integers starting with 2
       sum += 2/3*func(A+i*h)
       #sum=float32(sum)
    return h*sum

def trapezoidy(A,B,N):   #integrate from A to B using N points
    h = (B - A)/(N - 1)                     # step size
    sum = (funcy(A)+funcy(B))/2               # (1st + last)/2
    for i in range(1, N-1):        # i goes from 1 to (N-1)-1
       sum += funcy(A+i*h)
       #sum=float32(sum)            # to simulate single-precision (32 bit) calculation
    return h*sum

def simpsony(A,B,N):
    if ((N-1)%2==1):     #  if number of intervals odd
        print("Simpson's rule requires even number of intervals")
        return 0
    h = (B - A)/(N - 1)                     # step size
    sum = (funcy(A)+funcy(B))/3               # (1st + last)/3
    for i in range(1, N-1,2):        # i loops over odd integers  from 1 to (N-1)-1
       sum += 4/3*funcy(A+i*h)
       #sum=float32(sum)
    for i in range(2, N-1,2):        # i loops over even integers starting with 2
       sum += 2/3*funcy(A+i*h)
       #sum=float32(sum)
    return h*sum

def romberg(A,B,N):
    return (16*simpson(A,B,2*N-1) - simpson(A,B,N))/15

def rombergy(A,B,N):
    return (16*simpsony(A,B,2*N-1) - simpsony(A,B,N))/15

#num * 2 for num in lst

xpoints = np.arange(0.001,1,0.001)
fxpoints = [func(x) for x in np.arange(0.001,1,0.001)]

plot(xpoints, fxpoints)
title('Function '+r'$ cos(5x + x^2)/x^{3/2}$')
savefig('part5_funcx.png')
show()

ypoints = np.arange(0.001,1,0.001)
fypoints = [funcy(x) for x in np.arange(0.001,1,0.001)]
plot(ypoints, fypoints)
title('Function '+r'$ cos(y^{-4}+5y^{-2})$')
savefig('part5_funcy.png')
show()

i_simpsonx = simpson(1.0e-15, 1.0, 1000001)
#i_trapezoidx = trapezoid(1.0e-15, 1.0, 10000001)

i_rombergx = romberg(1.0e-15, 1.0, 1000001)

i_simpsony = simpsony(1.0e-15, 1.0, 1000001)
#i_trapezoidy = trapezoidy(1.0e-15, 1.0, 10000001)

i_rombergy = rombergy(1.0e-15, 1.0, 1000001)

print('Simpson method for x:', i_simpsonx)
print('Simpson method for y:', i_simpsony)

print('Romberg aproximation for Simpson method for x:', i_rombergx)
print('Romberg aproximation for Simpson method for y:', i_rombergy)

i_mathematica = 0.0329952

print(' Mathematica gives:', i_mathematica)
print(' Estimated error of Simpson:', abs(i_simpsony-i_mathematica))
print(' Estimated error of Romberg aproximation for Simpson method:', abs(i_rombergy-i_mathematica))

#print('Trapezoid method for x:', i_trapezoidx )
#print('Trapezoid method for y:', i_trapezoidy )
