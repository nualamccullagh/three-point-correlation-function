import numpy as N
import matplotlib.pyplot as plt

import threeptfn as T

def plot_realspace(b1, b2, bs):
    h=0.7
    rr=N.arange(105.0*h, 220.0*h, 1.0)
    zeta_eq_real=T.zeta_eq_real(rr, b1, b2, bs)
    print zeta_eq_real.size
    print rr.size
    plt.plot((rr/h), (rr/h)**3*(zeta_eq_real), 'k', label='Real-space')
    
    plt.xlabel('$r$ [Mpc/$h$]')
    plt.ylabel(r'$r^3 \zeta_{eq}(r)$')
    plt.legend()
    
    plt.show()

def make_plots(b1, b2, bs, f):
    
    h=0.7
    
    #range of r's to compute the 3pt function over, assuming r1=r2=r3
    rr=N.arange(105.0*h, 250.0*h, 1.0)
    
    #for now we will just take theta1 = theta2
    theta1=N.arange(30.0, 90.0, 10.0)
    
    #convert to radians
    theta1r=N.pi*theta1/180.0
    
    #convert to my coordinate system
    t1, phi1=T.mu1mu2_to_thetaphi(rr, rr, rr, N.cos(theta1r), N.cos(theta1r))
    
    zeta_eq=T.zeta_eq_rs(rr, t1, phi1, b1, b2, bs, f)
    zeta_eq_mon=T.zeta_eq_mon(rr, b1, b2, bs, f)
    plt.plot((rr/h), (rr/h)**3*zeta_eq_mon, 'k', label='RS Monopole')
    for i in N.arange(zeta_eq.shape[0]):
        plt.plot(rr/h, (rr/h)**3*zeta_eq[i, :], label=r'$\theta_1=\theta_2=%.1f$'%theta1[i])
    plt.xlabel('$r$ [Mpc/$h$]')
    plt.ylabel(r'$r^3 \zeta_{eq}(r)$')
    plt.ylim([-40000, 80000])
    plt.legend(ncol=1)
    plt.show()