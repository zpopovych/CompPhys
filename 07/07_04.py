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
    k = sqrt(abs(E-V1))   # wave number inside well
    r = sqrt(abs(V2-E))   # decay constant outside well
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
        if df > 1e-8: dx = - f(x)/df
        x += dx
        counter+=1
        if counter > 100: break
    #print('Number of iterrations:', counter)
    return x

V1 = 10.0  # bottom of well 
epsilon = 1e-8

V=[]
E = []
V_2_list = np.arange(70., 10000., .1)

for V2 in V_2_list:
    print(V2)
    V.append(V2)
    E_candidte = newton(25.)
    if (E_candidte < V2) and (E_candidte > V1):
        E.append(E_candidte)
    else:
        E.append(E[-1])
    #E.append(newton(V1+(V2-V1)/5))

loglog(V, E)
title('Energy dependence on the depth of the well')
xlabel(r'Potential $V_2$')
ylabel(r'Energy $E_2$')
savefig('07_04.png')
show()







