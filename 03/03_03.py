""" Example of curve fitting using scipy.optimize package
    Performs chi squared fit of data insample.txt file to qudratic function
    
    Lev Kaplan 2019"""

# Nonlinfit.py: Langrange interpolation tabulated data
    
from pylab import *
from scipy.optimize import curve_fit   # chi squared fitting

def func(x,a,b,c):   # define functional form
    return a + b*x + c*x*x

NMAX = 1000  # max number of input points

xin = zeros(NMAX)  # each is array of length NMAX, all elements set to zero
yin = zeros(NMAX)
sig = zeros(NMAX)

inputfile = open("lagrange.dat","r")  # read in the input x,y values
r = inputfile.readlines()  # read the whole file into list (one item per line)
inputfile.close()
        # input has the form: x0 y0 sig0
        #                     x1 y1 sig1
        #                     ...

m = 0
for line in r:
    print(line)
    s = line.split() # split line and split into list of items(assume items separated by spaces)
    xin[m] = s[0] # first number in each line is the x value
    yin[m] = s[1]
    sig[m] = s[2]
    m+=1         # m is total number of input data points
                 # will be stored in xin[0]..xin[n-1],yin[0]..yin[n-1]

popt,pcov = curve_fit(func, xin[0:m], yin[0:m], p0=[1,2,0], sigma=sig[0:m])

print("best fit parameters a,b,c =",popt)

print("covariance matrix for the parameters a,b,c =\n",pcov)

print("uncertainties in parameters =",sqrt(diag(pcov)))

xvalues = linspace(0,7,100)
yvalues = func(xvalues,popt[0],popt[1],popt[2])

errorbar(xin[0:m],yin[0:m],sig[0:m],fmt="o",label="input data")
plot(xvalues,yvalues,"b-",label="quadratic fit")
legend(loc="upper right")
show()
