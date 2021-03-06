"""
    Ten dimensional Monte Carlo inteaartion
    
    Lev Kaplan 2019
"""

# rand1.py: 10 dimensional monte calro integration

from pylab import *
from random import *

max = 2**16   # max number of points

nlist = []
values = []

y = 0.0
n = 1

for i in range(0,max+1):
    x = 0  # reset x
    for j in range(0,10):  # x1 + x2 + ... + x10
        x = x + random()
    y += x*x  # square and sum up
    
    if i == 2**n:  # sve average after 2, 4, 8, 16, ... points
        nlist.append(n)
        values.append(y/i)
        n +=1

        
    