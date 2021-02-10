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

NMAX = 100  # max number of input points

xin = zeros(NMAX)  # each is array of length NMAX, all elements set to zero
yin = zeros(NMAX)

inputfile = open("lagrange.dat","r")  # read in the input x,y values
r = inputfile.readlines()  # read the whole file into list (one item per line)
inputfile.close()
        # input has the form: x0 y0
        #                     x1 y1
        #                     ...

m = 0
for line in r:
    #print(line)
    s = line.split() # split line and split into list of items(assume items separated by spaces)
    xin[m] = s[0] # first number in each line is the x value
    yin[m] = s[1]
    print(xin[m],yin[m])
    m+=1         # m is total number of input data points
                 # will be stored in xin[0]..xin[m-1],yin[0].yin[m-1]

#xvalues=range(0,201,5)
# Extend the x range of the output to extrapolate up to E=220 MeV
# Make sure you are using a sufficiently fine grid for the output
# so that your output curve appears smooth, and plot the result along with the original 9 points.
xvalues=range(0,220,1)

yvalues = []

numpoints = m # use all points for interpolation
firstpoint = 0 # first point to use for interpolation

for x in xvalues:         # now interpolate      
    yvalues.append(legendrepol(x,firstpoint,firstpoint+numpoints-1))
    
plot(xin[0:m],yin[0:m],"o",label="input data")
plot(xvalues,yvalues,"-",label="Interpolated data")
title('Neutron scattering cross section \n (extrapolation of experimental data with all available points)')
xlabel('Energy of incident neutrons, '+ r'$E$ (MeV)')
ylabel('Cross section ' + r'$\sigma_S$'+'(particles/second)')
legend(loc="upper right")
savefig("a02p1_interpolation.png")
show()
