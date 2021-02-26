""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation
    
    adapted by Lev Kaplan in 2019 for linear decay fit y(x) = a + b x """
 
from pylab  import*

x = range(5,120,10)    # time from 5 to 115 in steps of 10 (12 points) - consider \delta t = 10
Nd = len(x)   # number of data points
y = log([32,17,21,7,8,6,5,2,2,0.1,4,1])   # log of number of counts - consider it \delta N

print("Number of data points:", Nd)

sig = [1] * Nd   # error bars all set to 1

plot(x, y, 'bo', label="experimantal data" )                                   # Plot data in blue

errorbar(x,y,sig)                                     # Plot error bars
title('Linear Least Squares Fit for Exponential Decay')                        # Plot figure
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
print('tau = ' + str(round(-1/slope/10,2))+' E-8 [s]')


# red line is the fit, red dots the fits at y[i]
t = range(0,120,1)
curve  = inter + slope*t
points = inter + slope*x
plot(t, curve,'r', label="OLS fit" )
plot(x, points, 'ro')
xlabel( 't [ns], '+r'$\tau \approx$ '+ str(round(-1/slope/10,2))+ r'$ \times 10^{-8}$ s.')                                           # Label axes
ylabel( r'$\ln{(\Delta N)}$' )
legend(loc="lower left")
savefig('03_01_OLS_fit.png')
show()

#Plot in linear scale

deltaN = [32,17,21,7,8,6,5,2,2,0.1,4,1]
bar(x, deltaN, label="experiment", width=10, alpha=0.5)
plot(t, exp(curve),'r', linestyle='dashed', label="fit" )
title('Exponential Decay \n (linear scale)')
ylabel( 'Number of decays detected, '+ r'$\Delta N$')
xlabel( 'Time, t [ns]')
legend(loc="upper right")
savefig('03_01_linear_scale.png')
show()
