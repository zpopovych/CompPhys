""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Simplified and adapted by Lev Kaplan 2019"""

# rk4.py 4th order Runge Kutta


from pylab import *
import numpy as np

ydumb = zeros(2);       y = zeros(2)
fReturn = zeros(2);     k1 = zeros(2, float)
k2 = zeros(2, float);   k3 = zeros((2), float) 
k4 = zeros((2), float)

def f( t, y):      # Force function: component 0 is position, component 1 is velocity
    fReturn[0] = y[1]       # d theta/dt  = omega
    fReturn[1] = - sin(y[0])     # d omega/dt = - (g/L) * sin (theta)
    return fReturn


def rk4(t,h,n):          # take one 4th order RK step          
                          # evolve y from time t to time t+h
                          # n is number of variables to evolve
                          # k1,k2,k3,k4 are four estimates for Delta y
    k1 = h*f(t, y)       # note that we use vector (array) notation instead of using loops
                         # see text for how to do this the long way with loops
    ydumb = y + k1/2     #estimate for midpoint using Euler
    k2 = h*f(t+h/2, ydumb) 
    ydumb = y + k2/2     #another estimate for midpoint
    k3 = h*f(t+h/2, ydumb)
    ydumb = y + k3       #estmate for endingn poitn of interval
    k4 = h*f(t+h, ydumb) 
    ynew = y + (k1 + 2*(k2 + k3) + k4)/6
    return ynew   

thetas = []
periods = []
theor_periods = []

for theta in np.arange(0.1, pi-0.1, 0.01):

    a = 0.  # evolve from time a to time b in n steps
    b = 100.
    n = 150
    t = a
    h = (b - a) / n

    y[0] = theta
    y[1] = 0  # initialize position and velocity

    time = []
    position = []

    roots = []
    pre_x = 0

    while (t < b):  # Time loop
        if ((t + h) > b):
            h = b - t  # Last step
        y = rk4(t, h, 2)
        t = t + h

        # Problem 1: Prepare lists of positions, velocities, energy and time
        x = y[0]
        v = y[1]
        position.append(x)
        time.append(t)

        # Problem 2: find roots of position function to evaluate period
        if len(position) > 1: pre_x = position[-2]
        if x * pre_x < 0: roots.append((x * t - pre_x * (t - h)) / (x - pre_x))
    thetas.append(theta)
    periods.append(2*average(np.diff(roots)))
    theor_periods.append( 2*pi*( 1 + (1/4)*(sin(theta/2))**2 + (9/64)*(sin(theta/2))**4) )

plot(thetas, periods, label='rk4', alpha=.5)
plot(thetas, theor_periods, label='analytical 4-th order', alpha=.5)
axhline(y=2*pi, linewidth=1, label='harmonic oscillator', alpha=.5, c ='green')
title('Period of pendulum oscillations ')
xlabel(r'Initial angle, $\theta$')
ylabel(r'Period, $T$')
legend()
savefig('06_03.png')
show()