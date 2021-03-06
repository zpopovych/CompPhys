"""
    Random number generation using linear congruent method
    
    Lev Kaplan 2019
"""

# rand1.py: experimenting with random numbers 

from pylab import *

def drand48():
    global rnd
    rnd = (0o273673163155 * rnd + 11) % 2**48    # 0o means octal notation
    return rnd/2**48   # return number between 0 and 1

def RANDU():
    global rnd
    rnd = 65539 * rnd % 2**31
    return rnd/2**31
    
rnd = 1   # set seed to 1

#N = 10000000
N = int(1e8)

xlist = []
ylist = []

x1 = RANDU()
for i in range(0,N):   # collect N pairs of adjacent random numbers
    x2 = x1
    x1 = RANDU()
    if x1<=0.01 and x2<0.01:
        xlist.append(x1)
        ylist.append(x2)
        
scatter(xlist,ylist,s=1)   # scatter plot with points of size 1
title('Random pairs generatred with RANDU()')
show()
    