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

m = 0
k = 0.4624  # manually adjusted
for line in r:
    print(line)
    s = line.split() # split line and split into list of items(assume items separated by spaces)
    xin[m] = s[0] # first number in each line is the x value
    yin[m] = s[1]
    sig[m] = k*sqrt(yin[m]) # was sig[m] = s[2]
    m+=1         # m is total number of input data points
                 # will be stored in xin[0]..xin[n-1],yin[0]..yin[n-1]

init_abc = [80,1,75]
init_abcd = [80,1,75,1]

popt_poly,pcov_poly = curve_fit(func_poly_quad, xin[0:m], yin[0:m], p0=init_abc, sigma=sig[0:m])

popt_abc,pcov_abc = curve_fit(func_abc, xin[0:m], yin[0:m], p0=init_abc, sigma=sig[0:m])

popt,pcov = curve_fit(func_abcd, xin[0:m], yin[0:m], p0=init_abcd, sigma=sig[0:m])

print("best fit parameters a,b,c =", popt_abc)
print("best fit parameters a,b,c,d =", popt)

print("covariance matrix for the parameters a,b,c =\n",pcov_abc)
print("uncertainties in parameters =",sqrt(abs(diag(pcov_abc))))

print("covariance matrix for the parameters a,b,c,d =\n",pcov)
print("uncertainties in parameters =",sqrt(abs(diag(pcov))))

xvalues = linspace(0,200,1000)
yvalues_poly = func_poly_quad(xvalues, popt_poly[0], popt_poly[1], popt_poly[2])
yvalues_abc = func_abc(xvalues, popt_abc[0], popt_abc[1], popt_abc[2])
yvalues = func_abcd(xvalues, popt[0], popt[1], popt[2], popt[3])

chi2_poly=0
yfit_poly = func_poly_quad(xin, popt_poly[0], popt_poly[1], popt_poly[2])
for i in range(0, m):
    sig2 = sig[i] * sig[i]
    chi2_poly += ((yin[i] - yfit_poly[i])**2)/sig2

chi2_abc=0
yfit_abc = func_abc(xin, popt_abc[0], popt_abc[1], popt_abc[2])
for i in range(0, m):
    sig2 = sig[i] * sig[i]
    chi2_abc += ((yin[i] - yfit_abc[i])**2)/sig2

chi2=0
yfit = func_abcd(xin, popt[0], popt[1], popt[2], popt[3])
for i in range(0, m):
    sig2 = sig[i] * sig[i]
    chi2 += ((yin[i] - yfit[i])**2)/sig2

print('quadratic polynomial fit chi^2 = ', chi2_poly)
print('a,b,c fit chi^2 = ', chi2_abc)
print('a,b,c,d fit chi^2 = ', chi2)

errorbar(xin[0:m],yin[0:m],sig[0:m],fmt="o",label="experiment data")
plot(xvalues,yvalues,"b-",label=r'$ g = A / (D + B {(E-C)}^2 )$', alpha=0.5)
plot(xvalues,yvalues_abc,"g--",label=r'$ g = A / (1 + B {(E-C)}^2 )$', alpha=0.5)
plot(xvalues,yvalues_poly,"y-",label=r'$ g = A + B x + C x^2 $')

title('Nonlinear Fit of Neutron Scattering Experimental Data \n (different functions)')
ylabel( 'Cross section, '+ r'$g(E_i)$ [mb]')
xlabel( 'Energy of neutron, $E$ [MeV].   ' + r' Best $\chi^2 \approx $' + str(round(chi2,2)) + r', best $k \approx 0.4624 $')
legend(loc="upper right")
savefig('03_03-04-05.png')
show()
