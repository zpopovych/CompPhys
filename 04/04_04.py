""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Adapted by Lev Kaplan 2019"""  

# Walk.py  Random walk in 2 dimensions
    
from pylab import *
from random import *

seed(None)                  # Seed generator, None => use system clock
nwalks = 100   # of walks (trials)
nsteps = 10000   # of stepsfor each walk

ravg=zeros(nsteps)
ravgt=zeros(nsteps)

for w in range(0,nwalks):   # iterate over trials
    x=0; y=0   # Start at origin

    for i in range(0, nsteps):
        x += (random() - 0.5)*2.                        # -1 =< dx =< 1
        y += (random() - 0.5)*2.                        # -1 =< dy =< 1
        r = sqrt(x*x+y*y)   # distance after i steps
        ravg[i] += r     # compute average distance after i steps

        ravgt[i] += sqrt(i)/sqrt(2)

ravg = ravg / nwalks
ravgt = ravgt / nwalks

plot(ravg, label='experimental')
plot(ravgt, label='theoretical', alpha=0.7)
legend()
title('Average distance from the origin as a function of N')
xlabel('N')
ylabel(r'$r_{avg}$')
#show()
savefig('04_04.png')