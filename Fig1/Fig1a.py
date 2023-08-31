# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 11:38:17 2022

@author: 17000
"""

import numpy as np;
from pymatgen import core;
from itertools import product;

with open('energy.txt','w') as file:
    file.write('R\t E1s \t E1p \t E2\t E2s\t E2p\n');
    
    
for R in [10*i for i in range(45,46)]:
    with open('energy.txt','a') as file:
        print('calculating R = '+str(R));
        file.write(str(R)+'\t');
        
        d = R/10;    
        Mat = d*np.eye(3);
        na,nb,nc = int(R//Mat[0,0]+1), int(R//Mat[1,1]+1), int(R//Mat[2,2]+1);
        L = R*1.2*2;
        pos_list = [];
        atom_list= [];
        for a,b,c in product(range(-na,na),range(-nb,nb),range(-nc,nc)):
            r1 = np.dot([a,b,c],Mat);
            if(np.linalg.norm(r1)<R):
                pos_list.append(r1/L+1/2);
                atom_list.append(1);
        lat = L*np.eye(3);
        pos = core.IStructure( lat, atom_list,pos_list)
        
        
        
        n = len(pos);
        a = 4;
        bH = -3.739011+2/np.sqrt(3/4)*(-1/4)*25.2749;  #fm
        bLi = -2.22#-2/np.sqrt(15/4)*(-3/4)*2.49;
        be = 4*(bH+bLi)/a**3*d**3*10**-5;
        typ = pos.atomic_numbers;
        bi  = be*np.ones(n)
        
        Rij = pos.distance_matrix;
        
        def f(kappa):
            M = np.eye(n)+(np.exp(-kappa*Rij)-np.eye(n))/(Rij+np.eye(n))*bi;
            res = np.linalg.eigvalsh(M);
            return res;
        
        for orbit in [0,1,4,6,9]:
            k0,k1 = 0,0.02;
            print('orbit = '+str(orbit));
            if(f(k0)[orbit]<0):
                while(k1-k0>10**-5):
                    km = (k0+k1)/2;
                    if(f(km)[orbit]>0):
                        k1 = km;
                    else:
                        k0 = km;
                Eb = (k0+k1)**2/4*24.046*10**3
                file.write(str(Eb)+'\t');
            else:
                file.write('0 \t');
        file.write('\n');
   
# =============================================================================
# x = np.ones(n)/np.sqrt(n);
# de,eig2 = 10,10;
# while(de>0.00001):
#     y = np.linalg.solve(M,x);
#     eig,eig2 = eig2, 1/np.dot(x,y);
#     x = y/np.linalg.norm(y);
#     print(1/eig2);
#     de = abs(eig2-eig)
# =============================================================================
