""" Example of curve fitting using scipy.optimize package
    Performs chi squared fit of data insample.txt file to qudratic function
    
    Lev Kaplan 2019"""

# Nonlinfit.py: Langrange interpolation tabulated data
    
from pylab import *
from scipy.optimize import curve_fit   # chi squared fitting

def func_poly_quad(x,a,b,c):   # define functional form
    return a + b*x + c*x*x

def func_abc(x,a,b,c):   # define functional form
    return a/(1+b*(x-c)**2)

def func_abcd(x,a,b,c,d):   # define functional form
    return a/(d+b*(x-c)**2)
    # return a + b*x + c*x*x

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

chi2k = []
kk = []

for k in np.arange(start=0.3, stop=1.0, step=0.01):
    #print(k, end="")

    m = 0
    for line in r:
        #print(line)
        s = line.split() # split line and split into list of items(assume items separated by spaces)
        xin[m] = s[0] # first number in each line is the x value
        yin[m] = s[1]
        sig[m] = k*sqrt(yin[m]) # was sig[m] = s[2]
        m+=1         # m is total number of input data points
                     # will be stored in xin[0]..xin[n-1],yin[0]..yin[n-1]


    init_abcd = [80,1,75,1]

    popt,pcov = curve_fit(func_abcd, xin[0:m], yin[0:m], p0=init_abcd, sigma=sig[0:m])

    xvalues = linspace(0,200,1000)
    yvalues = func_abcd(xvalues, popt[0], popt[1], popt[2], popt[3])

    chi2=0
    yfit = func_abcd(xin, popt[0], popt[1], popt[2], popt[3])
    for i in range(0, m):
        sig2 = sig[i] * sig[i]
        chi2 += ((yin[i] - yfit[i])**2)/sig2

    chi2k.append(chi2)
    kk.append(k)

kk=np.array(kk)
chi2k=np.array(chi2k)

print(kk.shape)
print(chi2k.shape)

plot(kk, chi2k, label= r'$\chi^2(k)$')
plot(kk, np.linspace(6.0, 6.0, num=kk.shape[0]), alpha=0.5, label=r'optimum $\chi^2 \approx 6.0 $')
plot(np.linspace(.463, .463, num=chi2k.shape[0]), chi2k, alpha=0.5, label=r'optimum $k \approx 0.463$')
title(r'$\chi^2(k)$')
ylabel( r'$\chi^2$')
xlabel( r'$k$')
plt.grid(b=True, which='both')
legend(loc="upper right")
savefig('03_chi_of_k.png')
show()
