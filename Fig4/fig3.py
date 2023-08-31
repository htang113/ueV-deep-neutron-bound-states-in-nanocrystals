# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 21:37:22 2022

@author: 17000
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 11:38:17 2022

@author: 17000
"""

import numpy as np;
from pymatgen import core;
from itertools import product;

R = 400;

d = R*0.16;    
Mat = d*np.eye(3);
na,nb,nc = int(R//Mat[0,0]+1), int(R//Mat[1,1]+1), int(R//Mat[2,2]+1);
L = R*1.2*2;
pos_list = [];
atom_list= [];
for a,b,c in product(range(-na,na+1),range(-nb,nb+1),range(-nc,nc+1)):
    r1 = np.dot([a,b,c],Mat);
    if(np.linalg.norm(r1)<=R):
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

for orbit in [0]:
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

kappa_s = (k0+k1)/2;
M = np.eye(n)+(np.exp(-kappa_s*Rij)-np.eye(n))/(Rij+np.eye(n))*bi;
e,w = np.linalg.eigh(M);
w_s = (w.T)[0];

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
        Eb = (k0+k1)**2/4*24.046*10**3;
kappa_p = (k0+k1)/2;
M = np.eye(n)+(np.exp(-kappa_p*Rij)-np.eye(n))/(Rij+np.eye(n))*bi;
e,w = np.linalg.eigh(M);
w_p = (w.T)[1];

Na, Nb, Nc = int(1.2*na),int(1.2*nb),int(1.2*nc);

s = np.zeros([2*Na+1,2*Nb+1,2*Nc+1]);
p = np.zeros([2*Na+1,2*Nb+1,2*Nc+1]);
M = np.zeros(3);
rec = 0;

for a,b,c in product(range(-Na,Na+1),range(-Nb,Nb+1),range(-Nc,Nc+1)):
    r1 = np.dot([a,b,c],Mat)+L/2;
    if(np.linalg.norm(r1)<R):
        s[a,b,c] = w_s[rec];
        p[a,b,c] = w_p[rec];
        M += s[a,b,c]*p[a,b,c]*r1;
        rec += 1;
    else:
        for i in range(n):
            ri = np.array([pos[i].x,pos[i].y,pos[i].z]);
            rij= np.linalg.norm(r1-ri);
            if(rij>0.01):
                cs = np.exp(-kappa_s*rij)/rij*be;
                cp = np.exp(-kappa_p*rij)/rij*be;
                s[a,b,c] += cs*w_s[i];
                p[a,b,c] += cp*w_p[i];
        M += s[a,b,c]*p[a,b,c]*r1;
    if(c==0 and b==0):
        print(str(a))

res = M/np.linalg.norm(s)/np.linalg.norm(p);
print(res)

mH = 1.0078;
mLi = 6.9410;
M_sphere = n*d**3/a**3*(mH+mLi)*4*1.66054*10**-27;
dE = 

import json;




