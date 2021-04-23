"""
Finding quantum eigenvalue in quantum well using bisection
Based on "A SURVEY OF COMPUTATIONAL PHYSICS", Python eBook Version
   by RH Landau, MJ Paez, and CC Bordeianu
   Copyright Princeton University Press, Princeton, 2011; Book  Copyright R Landau, 
   Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2011.
   Support by National Science Foundation 
   
Lev Kaplan 2019
"""

from pylab import *
from math import *

# find eigenvalues of a quantum particle in a one-dimensional box 
# the potential is V(x)=V1 for |x|<a, V(x)=V2 for |x|>a 
# work in units hbar=1, 2m=1, a=1 

def f(E):     # f(E)=0 when E is eigenvalue 
    k = sqrt(E-V1)   # wave number inside well 
    r = sqrt(V2-E)   # decay constant outside well 
    return k*tan(k)-r


def bisect(E1,E2):
    tolerance = 1e-12  # bisect until bracket becomes this small
    while E2-E1 > tolerance:
        Enew=(E1+E2)/2    # midpoint of bracket
        if f(Enew)*f(E1) > 0:    # f(Enew) has same sign as f(E1)
            E1 = Enew
        else:
            E2 = Enew
    return E1

def newton(x):
    counter = 0
    eps = 1e-14
    dx = 1e-14
    while abs(f(x)) > eps:
        df = (f(x + dx/2) - f(x - dx/2))/ dx
        dx = - f(x)/df
        x += dx
        counter += 1
    #print('Number of iterrations:', counter)
    return x, counter

print('Machine precission:', np.finfo(float).eps)

V1 = 10.0  # bottom of well 
V2 = 70.0  # top of well
epsilon = 1e-8

print('Enengy:', newton(25.)[0])

E0 = []
E = []
I = []

for E_i in np.arange(24.01, 27.0, 0.1):
    E0.append(E_i)
    Ec, Ic = newton(E_i)
    E.append(Ec)
    I.append(Ic)

plot(E0, I)
title('Convergence of Newton-Raphson algorithm')
xlabel(r'Starting point $E_0$')
ylabel('Number of iterations')
savefig('07_03.png')
show()




