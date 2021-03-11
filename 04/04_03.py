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
    global rndr
    rndr = 65539 * rndr % 2**31
    return rndr/2**31
    
rnd = 1   # set seed to 1
rndr = 1   # set seed to 1

#N = 10000000
N = int(1e8)

xlist_r = []
ylist_r = []
xlist_d = []
ylist_d = []

x1_r = RANDU()
x1_d = drand48()
Sd = 0
Sr = 0

for i in range(0,N):   # collect N pairs of adjacent random numbers
    x2_r = x1_r
    x2_d = x1_d
    x1_r = RANDU()
    x1_d = drand48()
    Sr += x1_r**2 * x2_r**2
    Sd += x1_d ** 2 * x2_d ** 2
    if x1_r<=0.01 and x2_r<0.01:
        xlist_r.append(x1_r)
        ylist_r.append(x2_r)
    if x1_d <= 0.01 and x2_d < 0.01:
        xlist_d.append(x1_d)
        ylist_d.append(x2_d)
Qr=Sr/N
Qd=Sd/N

print("Everage Q(RANDU) =", Qr)
print("Everage Q(drand48) =", Qd)

figure(figsize=(10,20))

subplot(2,1,1)

fnts = 18

scatter(xlist_r,ylist_r,s=1)   # scatter plot with points of size 1
title(r'Random pairs of RANDU(),  $Q_{avr}=$'+str(Qr), fontsize=fnts)

subplot(2,1,2)
scatter(xlist_d,ylist_d,s=1)   # scatter plot with points of size 1
title(r'Random pairs of drand48(), $Q_{avr}=$'+str(Qd), fontsize=fnts)

tight_layout()
# tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#show()
savefig('04_03.png')