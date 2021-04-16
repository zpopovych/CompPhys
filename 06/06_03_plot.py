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

a = 0.     # evolve from time a to time b in n steps
b = 30.
n = 200
t = a;       h = (b-a)/n;

y[0] = 0.1;   y[1] = 0    #initialize position and velocity

time = []
position = []
theor_pos = []

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
    time.append(t)

    # Problem 2: find roots of position function to evaluate period
    if len(position)>1: pre_x = position[-2]
    if x*pre_x <0: roots.append( (x*t - pre_x *(t-h))/(x-pre_x))

# Problem 1: Plot dynamics of position and energy

title('Comparison of pendulum and harmonic oscillator')
plot(time, position, '--', label='pendulum')
plot(time, theor_pos, alpha=.5, label='harmonic')
axhline(linewidth=1, c='black')
legend()
locs, labels = xticks()
locs = np.arange(0, 30, step=pi)
labels = [ str(n)+r'$\pi$' for n in range(0, int(30/pi)+1) ]
xticks(locs, labels)
savefig('06_03_div.png')
show()

print('Average period:', average(np.diff(roots)))
    
print('Final positions:', y)   #print out final position and velocity
                                