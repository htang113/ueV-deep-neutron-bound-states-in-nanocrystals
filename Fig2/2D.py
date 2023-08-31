# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:04:26 2023

@author: 17000
"""

import numpy as np;
from itertools import product;
import matplotlib.pyplot as plt;
from scipy.special import j0;

a = 4.00017; #lattice constant
thickness = 40; #nm
Nz = int(thickness*10//a);

pos0 = [(0,0,0),(1/2,1/2,1/2),(1/2,1/2,0),(0,0,1/2)];
polarization = 1;
bH = -3.739011+2/np.sqrt(3/4)*(-1/4)*25.2749*polarization;  #fm
bLi = -2.22#-2/np.sqrt(15/4)*(-3/4)*2.49;
b0 = [bH]*2+[bLi]*2;  #scattering length, unit fm
h2m = 2.07212*10**-3; #hbar^2/2m, unit eV*A^2

zl, bl = [],[];
for nz in range(Nz):
    for at in range(4):
        zl.append((pos0[at][2]+nz)*a);
        bl.append(b0[at]);
N_atom = len(zl);

Rc = 10;
klist = [(0.002*(10-i),0,0) for i in range(10)]+[(0,0.002*i/np.sqrt(2),0.002*i/np.sqrt(2)) for i in range(11)];
klist = [(0,0,0)];
El = [];

def core(z,kappa,k, accuracy = 10**-5):
    rc = -np.log(accuracy)/kappa;
    interp = int(1/np.sqrt(accuracy));
    rl = np.linspace(0,rc,interp)+rc/interp*0.5;
    res = sum([np.exp(-kappa*np.sqrt(z**2+r**2))/np.sqrt(z**2+r**2)*r*j0(k*r) for r in rl])
    return res*rc/interp;

for k0 in klist:
    k0 = np.array(k0);
    k0n = np.linalg.norm(k0);
    
    def M(kp):
        c_list = [core(a/2*i, kp, k0n) for i in range(2*Nz)];
        mat = np.zeros([N_atom,N_atom]);
        for alpha in range(N_atom):
            for beta in range(N_atom):
                zab = abs(zl[alpha] - zl[beta]);
                ind = int(round(zab/(a/2)))
                mat[alpha,beta] = 2*np.pi*bl[beta]*10**-5/(a**2/2)*c_list[ind];
                
        return mat;
    
    def crit(kp):
        mat = np.eye(N_atom)+M(kp);
        return np.real(np.linalg.eigvals(mat));
    
    for i in range(5):
        kl, kh = 0,0.05;
        if(crit(0.0001)[i]<0):
            while(kh-kl>10**-5):
                km = (kh+kl)/2;
                if(crit(km)[i]>0):
                    kh = km;
                else:
                    kl = km;
        else:
            km = 0;
    
        Eb = -h2m*km**2*10**6;
        El.append(Eb);
    print(El);
    
