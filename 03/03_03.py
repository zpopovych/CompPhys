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
# print(r)

m = 0
for line in r:
    print(line)
    s = line.split() # split line and split into list of items(assume items separated by spaces)
    xin[m] = s[0] # first number in each line is the x value
    yin[m] = s[1]
    sig[m] = sqrt(yin[m]) # was sig[m] = s[2]
    m+=1         # m is total number of input data points
                 # will be stored in xin[0]..xin[n-1],yin[0]..yin[n-1]



popt,pcov = curve_fit(func, xin[0:m], yin[0:m], p0=[1,2,0], sigma=sig[0:m])

print("best fit parameters a,b,c =", popt)

print("covariance matrix for the parameters a,b,c =\n",pcov)

print("uncertainties in parameters =",sqrt(diag(pcov)))

xvalues = linspace(0,200,1000)
yvalues = func(xvalues, popt[0], popt[1], popt[2])

chi2 = 0
yfit = func(xin, popt[0], popt[1], popt[2])
for i in range(0, m):
    sig2 = sig[i] * sig[i]
    chi2 += ((yin[i] - yfit[i]) ** 2) / sig2

print('chi^2 = ', chi2)

errorbar(xin[0:m],yin[0:m],sig[0:m],fmt="o",label="experiment data")
plot(xvalues,yvalues,"b-",label="quadratic fit")
title('Nonlinear Fit of Neutron Scattering Experimental Data \n'+r' with function  $y = A + B x + C x^2$')
ylabel( 'Cross section, '+ r'$g(E_i)$ [mb]')
xlabel( 'Energy of neutron, $E$ [MeV] '+ r'   $\chi^2 \approx $' + str(round(chi2,2)))
legend(loc="upper right")
savefig('03_03_quadratic_abc.png')
show()
