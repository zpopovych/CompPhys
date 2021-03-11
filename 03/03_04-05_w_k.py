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
sig1 = zeros(NMAX)

inputfile = open("lagrange.dat","r")  # read in the input x,y values
r = inputfile.readlines()  # read the whole file into list (one item per line)
inputfile.close()
        # input has the form: x0 y0 sig0
        #                     x1 y1 sig1
        #                     ...
# print(r)

m = 0
#k = 0.4624  # manually adjusted

k = 10

for line in r:
    print(line)
    s = line.split() # split line and split into list of items(assume items separated by spaces)
    xin[m] = s[0] # first number in each line is the x value
    yin[m] = s[1]
    sig[m] = k*sqrt(yin[m]) # was sig[m] = s[2]
    sig1[m] = 1 * sqrt(yin[m])  # was sig[m] = s[2]
    m+=1         # m is total number of input data points
                 # will be stored in xin[0]..xin[n-1],yin[0]..yin[n-1]

init_abc = [80,1,75]

popt_abc1, pcov_abc1 = curve_fit(func_abc, xin[0:m], yin[0:m], p0=init_abc, sigma=sig1[0:m])

popt_abc, pcov_abc = curve_fit(func_abc, xin[0:m], yin[0:m], p0=init_abc, sigma=sig[0:m])



print("best fit for k = 1 parameters a,b,c =", popt_abc1)
print("best fit for k = 0.4624 parameters a,b,c =", popt_abc)

print("for k = 1 covariance matrix for the parameters a,b,c =\n",pcov_abc1)
print("for k = 0.4624 covariance matrix for the parameters a,b,c =\n",pcov_abc)

print("for k = 1 uncertainties in parameters =",sqrt(abs(diag(pcov_abc1))))
print("for k = 0.4624 uncertainties in parameters =",sqrt(abs(diag(pcov_abc))))


xvalues = linspace(0,200,1000)

yvalues_abc = func_abc(xvalues, popt_abc[0], popt_abc[1], popt_abc[2])

yvalues_abc1 = func_abc(xvalues, popt_abc1[0], popt_abc1[1], popt_abc1[2])

chi2_abc=0
yfit_abc = func_abc(xin, popt_abc[0], popt_abc[1], popt_abc[2])
for i in range(0, m):
    sig2 = sig[i] * sig[i]
    chi2_abc += ((yin[i] - yfit_abc[i])**2)/sig2

chi2_abc1=0
yfit_abc1 = func_abc(xin, popt_abc1[0], popt_abc1[1], popt_abc1[2])
for i in range(0, m):
    sig2 = sig[i] * sig[i]
    chi2_abc1 += ((yin[i] - yfit_abc[i])**2)/sig2

print('k = 0.4624  fit chi^2 = ', chi2_abc)
print('k = 1 fit chi^2 = ', chi2_abc1)

errorbar(xin[0:m],yin[0:m],sig[0:m],fmt="o",label="experiment data")

plot(xvalues,yvalues_abc,"b--",label=r'$ k = $'+str(k), alpha=0.5)
plot(xvalues,yvalues_abc1,"g--",label=r'$ k = 1.0 $', alpha=0.5)

title('Fit with different k parameter of Neutron Scattering Experimental Data \n'+ r'$ g = A / (1 + B {(E-C)}^2 )$')
ylabel( 'Cross section, '+ r'$g(E_i)$ [mb]')
xlabel( 'Energy of neutron, $E$ [MeV].   ')
legend(loc="upper right")
savefig('03_04-05_w_k.png')
show()
