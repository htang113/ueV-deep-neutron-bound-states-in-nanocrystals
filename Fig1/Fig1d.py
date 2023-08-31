# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 11:38:17 2022

@author: 17000
"""

import numpy as np;
from pymatgen import core;
from itertools import product;

R = 300;

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

for orbit in [1]:
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

kappa = (k0+k1)/2;
M = np.eye(n)+(np.exp(-kappa*Rij)-np.eye(n))/(Rij+np.eye(n))*bi;
e,w = np.linalg.eigh(M);
wave = (w.T)[1];

with open('wave_p.txt','w') as file:
    for a,b in product(range(-2*na,2*na),range(-2*nb,2*nb)):
        r1 = np.dot([a,b,0],Mat)+L/2;
        s = 0;
        for i in range(n):
            ri = np.array([pos[i].x,pos[i].y,pos[i].z]);
            rij= np.linalg.norm(r1-ri);
            if(rij>0.01):
                s += np.exp(-kappa*rij)/rij*be*wave[i];
        file.write(str(r1[0]/10)+'\t'+str(r1[1]/10)+'\t'+str(s)+'\n');
        if(b==0):
            print(a);
