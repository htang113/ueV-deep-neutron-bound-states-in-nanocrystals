# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 11:38:17 2022

@author: 17000
"""

import numpy as np;
from pymatgen import core;
from itertools import product;
import matplotlib.pyplot as plt;

a = 4.00017; #lattice constant
coord  = [[1/2,0,0],[0,1/2,0],[0,0,1/2],[1/2,1/2,1/2],
          [0,0,0],[1/2,1/2,0],[1/2,0,1/2],[0,1/2,1/2]];

polarization = 1;
bH = -3.739011+2/np.sqrt(3/4)*(-1/4)*25.2749*polarization;  #fm
bLi = -2.22#-2/np.sqrt(15/4)*(-3/4)*2.49;

b = [bH]*4+[bLi]*4;  #scattering length, unit fm

fLi = 2.650221980500000;
fH = 2.475592140000000;
sLi = (1/np.sqrt(6.941*fLi))**0.5*0.179797;
sH = (1/np.sqrt(1.0079*fH))**0.5*0.179797;

sigma = [sH]*4+[sLi]*4;
natom = len(coord);

h2m = 2.07212*10**-3; #hbar^2/2m, unit eV*A^2

N = 4;
res = [];
klist = [(0.002*(10-i),0,0) for i in range(10)]+[(0,0.002*i/np.sqrt(2),0.002*i/np.sqrt(2)) for i in range(11)]
for k0 in klist:
    k0 = np.array(k0);
    hk = np.zeros([(2*N+1)**3]*2)*1j;
    for i,j,k in product(range(-N,N+1),range(-N,N+1),range(-N,N+1)):
        for l,m,n in product(range(-N,N+1),range(-N,N+1),range(-N,N+1)):
            x = i*(2*N+1)**2+j*(2*N+1)+k;
            y = l*(2*N+1)**2+m*(2*N+1)+n;
            if(hk[y,x] == 0):            
                G = np.array([i-l, j-m, k-n])*(2*np.pi)/a;
                if(x==y):
                    T = np.linalg.norm(np.array([i,j,k])*(2*np.pi/a)+k0)**2;
                else:
                    T = 0;
                
                V = 4*np.pi/a**3*10**-5;
                V *= sum([b[q]*np.exp(-(sigma[q]*np.linalg.norm(G))**2/2-1j*np.dot(coord[q],G)) for q in range(natom)]);
                hk[x,y] = T+V;
            else:
                hk[x,y] = np.conjugate(hk[y,x]);
            
    res.append( np.min(np.linalg.eigvalsh(hk))*h2m*10**6);
    print(k0);
plt.plot(res);

with open('direct.txt','w') as file:
    for i in range(21):
        file.write(str(0.02*(-1+0.1*i))+'\t'+str(res[i])+'\n');
       
