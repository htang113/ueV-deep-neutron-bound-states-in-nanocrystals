#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:29:58 2023

@author: ubuntu
"""

import numpy as np;
from itertools import product;
import matplotlib.pyplot as plt;
from scipy.special import j0;
import scipy;

a = 4.00017; #lattice constant
output = [];

R = 30;

Rc = R/2;

Radius = 10*R;

Nz = int(Radius//Rc+1);

polarization = 1;
bH = -3.739011+2/np.sqrt(3/4)*(-1/4)*25.2749*polarization;  #fm
bLi = -2.22#-2/np.sqrt(15/4)*(-3/4)*2.49;
b0 = (bH*4+bLi*4)/a**3*Rc**3;  #scattering length, unit fm
h2m = 2.07212*10**-3; #hbar^2/2m, unit eV*A^2

rl = [];
for nx in range(-Nz,Nz+1):
    for ny in range(-Nz,Nz+1):
        r = [Rc*nx,Rc*ny];
        if(np.linalg.norm(r)<Radius):
            rl.append([nx,ny]);

N_atom = len(rl);
rl = np.array(rl);

klist = [(0.002*i,0,0) for i in range(21)];

def core(z,kappa,k, accuracy = 10**-3):
    rc = -np.log(accuracy)/kappa;
    f = lambda r: np.exp(-kappa*np.sqrt(z**2+r**2))/np.sqrt(z**2+r**2)*2*np.cos(k*r);
    res = scipy.integrate.quad(f,a/2,rc)[0];
    if(z!=0):
        res = res + f(0)/2;
    return res;

for k0 in klist:
    El = [];
    k0 = np.array(k0);
    k0n = np.linalg.norm(k0);
    
    def M(kp):
        c_list = [[core(Rc*np.sqrt(i**2+j**2), kp, k0n) for j in range(i+1)] for i in range(2*Nz)];
        mat = np.zeros([N_atom,N_atom]);
        for alpha in range(N_atom):
            for beta in range(N_atom):
                zab = np.abs(rl[alpha]-rl[beta]);
                ind = [np.max(zab),np.min(zab)];

                mat[alpha,beta] = b0*10**-5/Rc*c_list[ind[0]][ind[1]];
                
        return mat;
    
    def crit(kp):
        mat = np.eye(N_atom)+M(kp);
        return np.real(np.linalg.eigvals(mat));
    
    for i in range(10):
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
    output.append(El);


    with open('1D_energy.txt','w') as file:
        file.write('r\t E1\t E2\t E3\t E4\t E5\n');
        for i in range(len(klist)):

            file.write(str(klist[i])+'\t');

            for j in range(10):
                file.write(str(output[i][j])+'\t');
            file.write('\n');

