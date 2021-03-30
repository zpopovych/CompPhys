"""
    Ten dimensional Monte Carlo inteaartion
    
    Lev Kaplan 2019
"""

# rand1.py: 10 dimensional monte calro integration

from pylab import *
from random import *

max = 2**16   # max number of points

nlist = []
values = []


y = 0.0
n = 1

for i in range(0,max+1):
    x = 0  # reset x
    for j in range(0,10):  # x1 + x2 + ... + x10
        x = x + random()
    y += x*x  # square and sum up
    
    if i == 2**n:  # sve average after 2, 4, 8, 16, ... points
        nlist.append(n)
        values.append(y/i)
        n +=1

theor = np.empty(len(nlist))
theor.fill(155/6)

error_values = np.abs(values - theor)
theor_error = 10./np.sqrt(2**np.array(nlist))

figure(figsize=(5,10))
subplot(2,1,1)

plot(nlist, values, label = 'Monte-Carlo numerical estimation')
plot(nlist, theor, label='theoretical value 155/6', alpha=0.5)

title('Monte-Carlo Estimation of Multidimentional Integral \n' +
      r'$\int_0^1 dx_1 \int_0^1 dx_2 ... \int_0^1 dx_{10} \left(x_1 + x_2 + ... + x_{10}\right)^2$')
ylabel('Value of integral ')
#xlabel(r'$log_{10} ( N )$')
legend()

subplot(2,1,2)
plot(nlist, error_values, label='actual error')
plot(nlist, theor_error, '--', label=r'theoretical error trend $1/\sqrt{N}$', alpha=0.5)
ylabel('Error ')
xlabel(r'$log_{2} ( N )$')
legend()

savefig('04_05_multi_int_error.png')
show()