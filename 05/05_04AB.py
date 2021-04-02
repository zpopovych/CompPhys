""" From "A SURVEY OF COMPUTATIONAL PHYSICS", Python eBook Version
   by RH Landau, MJ Paez, and CC Bordeianu
   Copyright Princeton University Press, Princeton, 2011; Book  Copyright R Landau, 
   Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2011.
   Support by National Science Foundation , Oregon State Univ, Microsoft Corp
   
   Adapted by Lev Kaplan 2019"""  

from pylab import *
from random import *

N     = 100                                         # number of spins
J     = 1                        # Exchange energy (J>0 ferromagnetic)
k     = 1.                                            # Boltmann constant
T     = 1.                                              # Temperature
B     = .001



state = array([+1]*N)                 # spin state: +1 up or -1 (down)
                                    # start with all spins up            
seed()                                     # Seed random generator

def energy (S):                                 # Method to calc energy
    sum = 0.  # Sum  energy (interaction term only)
    for  i in range(0,N-1):    # i goes from 0 to N-2
        sum += S[i]*S[i+1]

    return -J*sum + np.sum(S)*B

def averages(T, A, Bmax):

    state = array([+1] * N)  # spin state: +1 up or -1 (down)
    # start with all spins up

    ES = energy(state)  # energy of initial state
    mag = sum(state)  # initial magnetization
    # Here is the Metropolis algorithm
    E_avg = ES
    M_avg = mag
    E_sq_avg = ES ** 2

    a = 0
    b = 0

    for j in range(1, (A+Bmax)*100):
        r = int(N * random());  # randomly choose which spin to flip
        state[r] *= -1  # temporarily flip that spin
        #    ET = energy(state)             # finds energy of the test config.
        ET = ES  # find change in energy
        if (r >= 1): ET += -2 * J * state[r - 1] * state[r]
        if (r <= N - 2): ET += -2 * J * state[r] * state[r + 1]
        deltamag = 2 * state[r]  # find change in magnetization

        p = exp((ES - ET) / (k * T))  # test with Boltzmann factor
        if p < random():  # reject change
            state[r] *= -1  # go back and keep previous energy
        else:
            a += 1
            ES = ET  # update energy and magnetization
            mag += deltamag
            if a > A:
                b +=1
                if b > Bmax: break
                E_avg = (E_avg * (b - 1) + ES) / b
                M_avg = (M_avg * (b - 1) + mag) / b
                E_sq_avg = (E_sq_avg * (b - 1) + ES * ES) / b
            else:
                E_avg = ES
                M_avg = mag
                E_sq_avg = ES ** 2


    C = (E_sq_avg - E_avg ** 2) / (N * T ** 2)

    return E_avg, M_avg, E_sq_avg, C, b


E_avg = []
M_avg = []
E_sq_avg = []
C = []

E_avg_r = []
M_avg_r = []
E_sq_avg_r = []
C_r = []

E_theor = []
M_theor = []
C_theor = []


T = []
A = 5000

for t in np.arange(0.01, 5., 0.01):

    E, M, E_sq, c, b0 = averages(t, A=0, Bmax=1000 )
    Er, Mr, E_sqr, cr, b = averages(t, A=10000, Bmax=1000 )

    E_avg.append(E)
    M_avg.append(M)
    C.append(c)

    E_avg_r.append(Er)
    M_avg_r.append(Mr)
    C_r.append(cr)

    E_theor.append(-N*J*tanh(J/(k*t)))
    M_theor.append(N*exp(J/k/t)*sinh(B/k/t)/(sqrt(exp(2*J/k/t)*sinh(B/k/t)**2+exp(-2*J/k/t))))
    C_theor.append(((J/(k*t))/cosh(J/(k*t)))**2)

    T.append(t)

print('Relaxation of A= %1d steps' % A)
print('Averaging over B= %1d steps' % b)


figure(figsize=(5,8))
subplot(3,1,1)
plot(T, E_avg, alpha=.7)
plot(T, E_avg_r, alpha=.7)
plot(T, E_theor)
title(r'Everage energy $\langle E \rangle$')
subplot(3,1,2)
plot(T, M_avg, alpha=.7)
plot(T, M_avg_r, alpha=.7)
plot(T, M_theor)
title(r'Everage magnetization $\mathcal{M}$')
subplot(3,1,3)
plot(T, C, label='sim w/o relax', alpha=.7)
plot(T, C_r, label='sim w relax A='+str(A), alpha=.7)
plot(T, C_theor, label='theoretical')
title(r'Specific heat $C = \frac{\langle E^2 \rangle - \langle E \rangle}{N T^2}$')
xlabel(r'Temperature, $T$')
ylim(0,1)
legend(loc='center right')
tight_layout()
savefig('05_04AB.png')
show()