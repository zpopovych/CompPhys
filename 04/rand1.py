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
    
rnd = 1   # set seed to 1

N = 10000000

xlist = []
ylist = []

x1 = drand48()
for i in range(0,N):   # collect N pairs of adjacent random numbers
    x2 = x1
    x1 = drand48()
    if x1<=0.01 and x2<0.01:
        xlist.append(x1)
        ylist.append(x2)
        
scatter(xlist,ylist,s=1)   # scatter plot with points of size 1
show()
    