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
    return exp(-x)

def trapezoid(A,B,N):   #integrate from A to B using N points
    h = (B - A)/(N - 1)                     # step size 
    sum = (func(A)+func(B))/2               # (1st + last)/2
    for i in range(1, N-1):        # i goes from 1 to (N-1)-1
       sum += func(A+i*h)
       sum=float32(sum)            # to simulate single-precision (32 bit) calculation
    return h*sum

def trapezoid_d(A,B,N):   #integrate from A to B using N points
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
       sum=float32(sum)
    for i in range(2, N-1,2):        # i loops over even integers starting with 2
       sum += 2/3*func(A+i*h)
       sum=float32(sum)
    return h*sum

def simpson_d(A,B,N):
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
              
A = 0.0
B = 1.0

maxpoints = 1e7

Nvalues = []    #empty lists
traperror = []
simpson_error = []

traperror_d = []
simpson_error_d = []

exact=1-exp(-1)

N=3
while N<maxpoints:    # loop over number of points
    print(N)    #just so you can see while code is running
    Nvalues.append(N)
    traperror.append(abs(trapezoid(A,B,N)-exact)/exact) # Error in trapezoid method
    traperror_d.append(abs(trapezoid_d(A, B, N) - exact) / exact)  # Error in trapezoid method

    simpson_error.append(abs(simpson(A, B, N) - exact) / exact)  # Error in simpson method
    simpson_error_d.append(abs(simpson_d(A, B, N) - exact) / exact)  # Error in simpson method

    N=int(N*1.1)+1    # N grows roughly by 1.1 factor each time
    if N%2 == 0:
        N=N+1    # make sure N is odd
        
loglog(Nvalues,traperror,label="Trapezoid (single )")   # log log plot of error in trapezoid method
loglog(Nvalues,simpson_error,label="Simpson (single)")   # log log plot of error in simpson method

loglog(Nvalues,traperror_d,label="Trapezoid (double)")   # log log plot of error in trapezoid method
loglog(Nvalues,simpson_error_d,label="Simpson (double)")   # log log plot of error in simpson method

loglog(Nvalues,0.1*array(Nvalues)**(-2.0),label="0.1/N^2", alpha=.5)   # plot 0.1/N^2 for comparison
loglog(Nvalues,0.07*array(Nvalues)**(-5.0),label="0.07/N^5", alpha=.5)   # plot 0.1/N^5 for comparison
loglog(Nvalues,0.005*array(Nvalues)**(-4.0),label="0.005/N^4", alpha=.5)   # plot 0.1/N^5 for comparison

ylim([1e-21,1])   # set range of y values
title('Numerical integration error for different methods \n (single vs double precession)')
xlabel('Number of points')
ylabel('Total error')
legend(loc="upper right")
savefig('Integation_'+str(maxpoints)+"_points_sing_vs_double.png")
show()    #show the graph
