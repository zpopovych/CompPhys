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

gama = 0.05

def f( t, y):      # Force function: component 0 is position, component 1 is velocity
    #gama = 0.01
    Fdr = .1
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

phase_delta = []
theor_phase_delta = []

omegas = []
periods = []
theor_periods = []

for omega_dr in np.arange(0.5, 1.5, 0.0001):

    a = 0.  # evolve from time a to time b in n steps
    b = 200.
    n = 1000
    t = a
    h = (b - a) / n

    y[0] = 0
    y[1] = 0  # initialize position and velocity

    time = []
    position = []

    roots = []

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
        if len(position) > int(n/2):
            pre_x = position[-2]
            if x * pre_x < 0:
                if x > pre_x: up_flag = True
                else: up_flag = False
                roots.append((x * t - pre_x * (t - h)) / (x - pre_x))

    omegas.append(omega_dr)
    periods.append(2*average(np.diff(roots)))
    theor_periods.append(6/omega_dr+.2)

    if up_flag: back_step=1
    else: back_step = 2
    phase_pend = 2*pi*((time[-1] - roots[-back_step])%periods[-1])/periods[-1]
    #print(phase_pend)
    delta = phase_pend - (omega_dr*time[-1])%(2*pi)
    if delta > pi: delta = delta - 2*pi
    if delta < -pi: delta = delta + 2*pi
    phase_delta.append(abs(delta))
    if omega_dr < 1 : theor_phase_delta.append(abs(arctan(gama*omega_dr/(omega_dr**2-1))))
    if omega_dr > 1: theor_phase_delta.append(pi - abs(arctan(gama * omega_dr / (omega_dr ** 2 - 1))))


plot(omegas, periods, marker=',', linestyle=None, label='rk4', alpha=.5)
plot(omegas, theor_periods, label='trend', alpha=.3)
title('Period of pendulum oscillations under external driving ')
xlabel(r'Driving frequency, $\omega_{dr}$')
ylabel(r'Period, $T$')
savefig('06_04.png')
show()

plot(omegas, phase_delta, marker=',', linestyle='None')
plot(omegas, theor_phase_delta, label='trend', alpha=.3)
#plot(omegas, phase_drv, label='driving', alpha=.5)
title(r'Relative phase of pendulum from driving, with damping $\gamma$ = '+str(gama))
xlabel(r'Driving frequency, $\omega_{dr}$')
ylabel(r'Phase difference, $|\Delta \varphi|$')
savefig('06_04_phases_g'+str(gama*100)+'.png')
show()