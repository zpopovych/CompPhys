""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    Adapted by Lev Kaplan 2019"""

# Lagrange.py: Langrange interpolation tabulated data
    
from pylab import *

def legendrepol (x,beg,finish):          # poly interpolation at x 
    y = 0.                               # using input points from beg to finish
    for  i in range(beg,finish+1): 
       lambd = 1.0;
       for j in range(beg,finish+1):
           if i != j:                       #Lagrange polynom formed here
              lambd = lambd * ((x - xin[j])/(xin[i] - xin[j]))
       y += yin[i] * lambd
    return y

xin = np.arange(-1, 1, .1)
yin = 1.0/(1.0 + 25*xin**2)

m = xin.size-1

xvalues=np.arange(-1.1, 1.1, .01)

yvalues = []
yexact = []

numpoints = 2 # use linear interpolation
firstpoint = 0 # first point to use for interpolation

point_num = 0

for x in xvalues:         # now interpolate
    yvalues.append(legendrepol(x, point_num, point_num + 1))
    yexact.append(1.0 / (1.0 + 25 * x ** 2))
    if (point_num < m-1)&(x > xin[point_num+1]): point_num+=1

plot(xvalues,yexact ,"-",label="exact function")
plot(xin[0:m],yin[0:m],"o",label="points used for interpolation")
plot(xvalues,yvalues,"-",label="linear interpolation")
title(r'Extrapolation of $\frac{1}{1+25x^2}$')
xlabel(r'$x$')
ylabel(r'$y$')
ylim(-2.5,2.5)
xlim(-1.1,1.1)
legend(loc="lower center")
savefig("a02p4_2_interpolation.png")
show()

plot(xvalues,yexact ,"-",label="exact function")
plot(xin[0:m],yin[0:m],"o",label="points used for interpolation")
plot(xvalues,yvalues,"-",label="linear interpolation")
title(r'Extrapolation of $\frac{1}{1+25x^2}$')
xlabel(r'$x$')
ylabel(r'$y$')
ylim(-0.5,1.25)
xlim(-1.1,1.1)
legend(loc="lower center")
savefig("a02p4_2_magnified_interpolation.png")
show()
