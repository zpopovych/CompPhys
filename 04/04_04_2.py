""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Adapted by Lev Kaplan 2019"""  

# Walk.py  Random walk in 2 dimensions
    
from pylab import *
from random import *

from numpy.random import choice

seed(None)                  # Seed generator, None => use system clock
nwalks = 100   # of walks (trials)
nsteps = 10000   # of stepsfor each walk

ravg=zeros(nsteps)
ravgt=zeros(nsteps)

for w in range(0,nwalks):   # iterate over trials
    x=0; y=0   # Start at origin

    for i in range(0, nsteps):

        directions = [[0., 1.], [0., -1.], [1., 0.], [-1., 0.]]
        dir_choice = choice(range(4))
        dx = directions[dir_choice][0]
        dy = directions[dir_choice][1]

        x += dx                        # -1 =< dx =< 1
        y += dy                        # -1 =< dy =< 1
        r = sqrt(x*x+y*y)   # distance after i steps
        ravg[i] += r     # compute average distance after i steps

        ravgt[i] += sqrt(i)*.9

ravg = ravg / nwalks
ravgt = ravgt / nwalks

plot(ravg, label='experimental')
plot(ravgt, label=r'theoretical: $r_{avg} = \sqrt{\frac{N}{2}}$', alpha=0.7)
legend()
title('Average distance from the origin as a function of N')
xlabel('N')
ylabel(r'$r_{avg}$')
show()
#savefig('04_04_2.png')