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
T     = 100.                                              # Temperature
state = array([+1]*N)                 # spin state: +1 up or -1 (down)  
                                    # start with all spins up            
seed()                                     # Seed random generator

def energy (S):                                 # Method to calc energy
    sum = 0.                       # Sum  energy (interaction term only)
    for  i in range(0,N-1):    # i goes from 0 to N-2
        sum += S[i]*S[i+1]
    return -J*sum 
		
ES = energy(state)                  # energy of initial state
mag = sum(state)                 # initial magnetization
                                    # Here is the Metropolis algorithm
for j in range (1,200):
    r = int(N*random());          # randomly choose which spin to flip
    state[r] *= -1                 # temporarily flip that spin
#    ET = energy(state)             # finds energy of the test config.
    ET = ES                           # find change in energy
    if (r>=1): ET+=  -2*J*state[r-1]*state[r]
    if (r<=N-2): ET+=  -2*J*state[r]*state[r+1]
    deltamag = 2*state[r]    #find change in magnetization
    
    p = exp((ES-ET)/(k*T))        # test with Boltzmann factor
    if p < random():             # reject change
        state[r] *= -1    # go back and keep previous energy
    else:
        ES = ET       #update energy and magnetization
        mag+=deltamag
