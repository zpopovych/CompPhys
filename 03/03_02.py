""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    adapted by Lev Kaplan in 2019 for linear decay fit y(x) = a + b x """
 
from pylab  import*


x = range(5,120,10)    # time from 5 to 115 in steps of 10 (12 points)
Nd = len(x)   # number of data points
y = log([32,17,21,7,8,6,5,2,2,0.1,4,1])   # log of number of counts
sig = 1 / sqrt(np.array([32, 17, 21, 7, 8, 6, 5, 2, 2, 0.1, 4,1]))  # error bars all set to 1/sqrt(n)
print(sig)

print("Number of data points:", Nd)

plot(x, y, 'bo', label="experimantal data" )         # Plot data in blue

errorbar(x,y,sig)                                     # Plot error bars
title(r'Linear Least Squares Fit with adjusted uncertainty $\sigma$')                        # Plot figure
xlabel( 't [ns]' )                                            # Label axes
ylabel( r'$\ln{(\Delta N)}$' )
grid(True)                                               # plot grid
xlim(0,120)                                              # x range for plot


ss = sx = sxx = sy = sxy = 0   # initialize various sums

for i in range(0, Nd):         # compute various sums over data points                              
    sig2 = sig[i] * sig[i]
    ss += 1. / sig2;    sx += x[i]/sig2;        sy += y[i]/sig2
    sxx += x[i] * x[i]/sig2;    sxy += x[i]*y[i]/sig2;
         
delta = ss*sxx-sx*sx;
slope = (ss*sxy-sx*sy) / delta    #calculate best fit slope
inter = (sxx*sy-sx*sxy) / delta   # calculate best fit intercept
      
print('Linear Fit Final Results\n') 
print('y(x) = a + b x')                          # Desired fit
print('a = ', inter, '+/-', sqrt(sxx/delta))                  
print('b = ', slope, '+/-', sqrt(ss/delta))
print('correlation =',-sx/sqrt(sxx*ss))

print('tau = ' + str(-1/slope)+'[ns]')
print('[ns] = 1E-9 [s]')

tau = -1/slope/10
delta_slope = sqrt(ss/delta)
delta_tau = (delta_slope/slope) * tau
print('tau = ' + str(round(tau,2))+'+/-'+ str(round(abs(delta_tau),2)) + ' E-8 [s]')

chi2=0

for i in range(0, Nd):
    sig2 = sig[i] * sig[i]
    chi2 += ((y[i] - (inter + slope * x[i]))**2)/sig2

print('chi^2 = ', chi2)

# red line is the fit, red dots the fits at y[i]
t = range(0,120,1)
curve  = inter + slope*t
points = inter + slope*x

plot(t, curve,'r', label="OLS fit" )
plot(x, points, 'ro')
legend(loc="lower left")
savefig('03_02_OLS_fit_w_sigma_adjusted.png')
show()
