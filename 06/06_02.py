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
    fReturn[0] = y[1]       # dx/dt  = v                                     
    fReturn[1] = -y[0]       # dv/dt = -x
    return fReturn

def pendulum_f( t, y):      # Force function: component 0 is position, component 1 is velocity
    g = 1; L = 1
    fReturn[0] = y[1]       # d theta/dt  = omega
    fReturn[1] = - (g/L) * sin(y[0])     # d omega/dt = - (g/L) * sin (theta)
    return fReturn


def energ(v,x):
    k=1
    m=1
    return (m*v**2 + k * x**2)/2


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

a = 0.     # evolve from time a to time b in n steps
b = 65.
n = 200
t = a;       h = (b-a)/n;

y[0] = 1;   y[1] = 0    #initialize position and velocity

time = []
position = []
theor_pos = []
velocity = []
theor_vel = []
energy = []
theor_energy = []
roots = []


pre_x = 0

while (t < b):                                              # Time loop
    if ((t + h) > b):
        h = b - t                                           # Last step
    y = rk4(t,h,2)
    t = t + h

    # Problem 1: Prepare lists of positions, velocities, energy and time
    x = y[0]
    v = y[1]
    position.append(x)
    theor_pos.append(sin(pi/2+t))
    velocity.append(v)
    theor_vel.append(sin(pi+t))
    energy.append(energ(v,x))
    theor_energy.append(.5)
    time.append(t)

    # Problem 2: find roots of position function to evaluate period
    if len(position)>1: pre_x = position[-2]
    if x*pre_x <0: roots.append( (x*t - pre_x *(t-h))/(x-pre_x))

# Problem 1: Plot dynamics of position and energy
figure(figsize=(8,3))
subplot(211)
tight_layout(pad=2.0)
title(r'Energy $E = \frac{m v^2}{2} + \frac{k x^2}{2}$')
plot(time, energy, '--', label='rk4' )
plot(time, theor_energy, label='analitical', alpha=.5)
legend(loc='lower right')
ylim(.49, .51)
subplot(212)
title('Position dynamics')
plot(time, position, '--')
plot(time, theor_pos, alpha=.5)
axhline(linewidth=1, c='black')
axvline(x=pi/2, c='green', alpha=.5)
axvline(x=3*pi/2, c='green', alpha=.5)
axvline(x=2*pi+pi/2, c ='green', alpha=.5, label='roots')
axvline(x=3*pi+pi/2, c ='green', alpha=.5)
axvline(x=4*pi+pi/2, c ='green', alpha=.5)
axvline(x=5*pi+pi/2, c ='green', alpha=.5)
axvline(x=6*pi+pi/2, c ='green', alpha=.5)
axvline(x=7*pi+pi/2, c ='green', alpha=.5)
axvline(x=8*pi+pi/2, c ='green', alpha=.5)
legend(loc='lower right')
locs, labels = xticks()
locs = np.arange(0, 60, step=pi)
labels = [ str(n)+r'$\pi$' for n in range(0, int(60/pi)+1) ]
xticks(locs, labels)
savefig('06_02.png')
show()


# Problem 2: caclulate average period
periods = np.diff(roots)
T = average(periods)
print('Average period:', T)
    
print('Final positions:', y)   #print out final position and velocity
                                