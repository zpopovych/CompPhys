""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation

    Adapted by Lev Kaplan 2019"""

# Lagrange.py: Langrange interpolation tabulated data

from pylab import *
import numpy as np


def legendrepol(x, beg, finish):  # poly interpolation at x
    y = 0.  # using input points from beg to finish
    for i in range(beg, finish + 1):
        lambd = 1.0;
        for j in range(beg, finish + 1):
            if i != j:  # Lagrange polynom formed here
                lambd = lambd * ((x - xin[j]) / (xin[i] - xin[j]))
        y += yin[i] * lambd
    return y


xin = np.arange(-1, 1, .1)
yin = 1.0 / (1.0 + 25 * xin ** 2)

m = xin.size

xvalues = np.arange(-1.1, 1.1, .01)

yvalues = []
yexact = []

numpoints = m  # use all points for interpolation
firstpoint = 0  # first point to use for interpolation

for x in xvalues:  # now interpolate
    yvalues.append(legendrepol(x, firstpoint, firstpoint + numpoints - 1))
    yexact.append(1.0 / (1.0 + 25 * x ** 2))

# plot(xvalues,yexact ,"-",label="exact function")
# plot(xin[0:m],yin[0:m],"o",label="points used for interpolation")
plot(xvalues, abs(array(yvalues) - array(yexact)), "-", label="Error in interpolation")
title('Error in Global Lagrange Interpolation \n'+r'of function $\frac{1}{1+25x^2}$')
xlabel(r'$x$')
ylabel(r'$Error, \epsilon$')
ylim(-0.1, 2.5)
xlim(-1, 1)
legend(loc="upper right")
savefig("a02p4_1_error_epsilon_vs_x.png")
show()