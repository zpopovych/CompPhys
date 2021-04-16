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
    gama = 0
    Fdr = .01
    fReturn[0] = y[1]       # d theta/dt  = omega
    fReturn[1] = - sin(y[0]) - gama*y[1] + Fdr*sin(omega_dr*t)    # d omega/dt = - (g/L) * sin (theta)
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

omegas = []
phase_delta = []
phase_delta_theor = []

for omega_dr in np.arange(0.0, 2., 0.01):

    a = 0.  # evolve from time a to time b in n steps
    b = 10.
    n = 100
    t = a
    h = (b - a) / n

    init_position = 0.1

    y[0] = init_position
    y[1] = 0  # initialize position and velocity

    time = []
    position = []

    while (t < b):  # Time loop
        if ((t + h) > b):
            h = b - t  # Last step
        y = rk4(t, h, 2)
        t = t + h
        x = y[0]

        position.append(x)
        time.append(t)

    omegas.append(omega_dr)
    phase_delta.append((init_position*position[n-2])%(2*pi) - (omega_dr*time[n-2])%(2*pi))

plot(omegas, phase_delta, label='rk4', alpha=.5)
title('Relative pahse of pendulum from driving ')
xlabel(r'Driving frequency, $\omega_{dr}$')
ylabel(r'Phase difference, $\Delta \varphi$')
savefig('06_04_phases.png')
show()