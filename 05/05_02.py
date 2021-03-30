""" From "A SURVEY OF COMPUTATIONAL PHYSICS", Python eBook Version
   by RH Landau, MJ Paez, and CC Bordeianu
   Copyright Princeton University Press, Princeton, 2011; Book  Copyright R Landau, 
   Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2011.
   Support by National Science Foundation , Oregon State Univ, Microsoft Corp
   
   Adapted by Lev Kaplan 2019"""  

from pylab import *
from random import *

N     = 100                                         # number of spins
J     = 1                        # Exchange energy (J>0 ferromagnetic)
k     = 1.                                            # Boltmann constant

seed()                                     # Seed random generator

def energy (S):                                 # Method to calc energy
    sum = 0.                       # Sum  energy (interaction term only)
    for  i in range(0,N-1):    # i goes from 0 to N-2
        sum += S[i]*S[i+1]
    return -J*sum

def energy_dynanics (T):

    state = array([+1] * N)  # spin state: +1 up or -1 (down)
    # start with all spins up

    upspinsx = []  # record up and down spins for plot
    downspinsx = []
    upspinst = []
    downspinst = []

    ES = energy(state)  # energy of initial state
    mag = sum(state)  # initial magnetization
    E = []
    E.append(ES)
    # Here is the Metropolis algorithm

    for j in range(1, 5000):  # was 200 steps
        for i in range(0, N):
            if state[i] == 1:
                upspinsx.append(i)  # record position and "time" for each up spin
                upspinst.append(j)
            else:
                downspinsx.append(i)  # and for each down spin
                downspinst.append(j)
        r = int(N * random());  # randomly choose which spin to flip
        state[r] *= -1  # temporarily flip that spin
        ET = energy(state)  # finds energy of the test config.
        p = exp((ES - ET) / (k * T))  # test with Boltzmann factor
        if p < random():  # reject change
            state[r] *= -1  # go back
        else:
            E.append(ES)
            ES = ET  # update energy and magnetization
            mag = sum(state)
    return (E)


E100 = 	energy_dynanics(100.)
E1 = energy_dynanics(1.)
title('Energy dynamic')
plot(E100[:1000], label=r'$T=100$')
plot(E1[:1000], label=r'$T=1$')
xlabel('steps')
ylabel(r'Energy $E(t)$')
legend(loc='center right')
savefig('05_02.png')
show()