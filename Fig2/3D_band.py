# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:04:26 2023

@author: 17000
"""

import numpy as np;
from itertools import product;
import matplotlib.pyplot as plt;

a = 4.00017; #lattice constant
caa = [(0,0,0),(0,1/2,1/2),(1/2,0,1/2),(1/2,1/2,0)];
cab = [(1/2,0,0),(0,1/2,0),(0,0,1/2),(1/2,1/2,1/2)];

polarization = 1;
bH = -3.739011+2/np.sqrt(3/4)*(-1/4)*25.2749*polarization;  #fm
bLi = -2.22#-2/np.sqrt(15/4)*(-3/4)*2.49;

b = [bH]*4+[bLi]*4;  #scattering length, unit fm
sigma = [0.1]*4+[0.1]*4;

h2m = 2.07212*10**-3; #hbar^2/2m, unit eV*A^2

Rc = 10;
res = [];
klist = [(0.002*(10-i),0,0) for i in range(10)]+[(0,0.002*i/np.sqrt(2),0.002*i/np.sqrt(2)) for i in range(11)]
klist = [(0,0,0)]
El = [];

for k0 in klist:
    k0 = np.array(k0);
    k0n = np.linalg.norm(k0);
    def M(kp,key):
        sm = 0;
        if(key == 'aa'):
            cl = caa;
        else:
            cl = cab;
        for i,j,k in product(range(-Rc-1,Rc+2),range(-Rc-1,Rc+2),range(-Rc-1,Rc+2)):
            for c in cl:
                r_vec = a*(np.array([i,j,k])+np.array(c));
                r = np.linalg.norm(r_vec);
                if(r<a*Rc and np.abs(r)>10**-3):
                    sm += np.exp(1j*np.dot(k0,r_vec))*np.exp(-kp*r)/r;
        if(k0n != 0):
            sm += 4*np.pi*4/a**3*np.exp(-kp*Rc*a)/(k0n**2+kp**2)*(np.cos(k0n*Rc*a)+kp/k0n*np.sin(k0n*Rc*a));
        else:
            sm += 4*np.pi*4/a**3*np.exp(-kp*Rc*a)/(k0n**2+kp**2)*(np.cos(k0n*Rc*a)+kp*Rc*a);

        return np.real(sm);
    
    def crit(kp):
        mat = np.matrix([[1+bH*M(kp,'aa')*10**-5,bLi*M(kp,'ab')*10**-5],
                         [bH*M(kp,'ab')*10**-5,1+bLi*M(kp,'aa')*10**-5]]);
        return np.linalg.det(mat);
    
    kl, kh = 0,0.05;
    if(crit(0.0001)<0):
        while(kh-kl>10**-5):
            km = (kh+kl)/2;
            if(crit(km)>0):
                kh = km;
            else:
                kl = km;
    else:
        km = 0;
    
    Eb = -h2m*km**2;
    El.append(Eb);
    print(k0);
    
# =============================================================================
# plt.plot(res);
# 
# with open('Eq3.txt','w') as file:
#     for i in range(21):
#         file.write(str(0.02*(-1+0.1*i))+'\t'+str(El[i])+'\n');
# =============================================================================
       
