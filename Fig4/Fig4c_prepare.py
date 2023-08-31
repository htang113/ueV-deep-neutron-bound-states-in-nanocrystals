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

d = R*0.1;    
Mat = d*np.eye(3);
na,nb,nc = int(R//Mat[0,0]+1), int(R//Mat[1,1]+1), int(R//Mat[2,2]+1);
L = R*1.2*2;
pos_list = [];
atom_list= [];
for a,b,c in product(range(-na,na+1),range(-nb,nb+1),range(-nc,nc+1)):
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

kappa = [];
wl = [];
for orbit in range(10):
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
        kappa.append((k0+k1)/2);
        M = np.eye(n)+(np.exp(-kappa[-1]*Rij)-np.eye(n))/(Rij+np.eye(n))*bi;
        e,w = np.linalg.eigh(M);
        wl.append((w.T)[orbit]);

Na, Nb, Nc = int(1.2*na),int(1.2*nb),int(1.2*nc);

wave = np.zeros([10,2*Na+1,2*Nb+1,2*Nc+1]);
M = np.zeros([10,10,3]);
rec = 0;

for a,b,c in product(range(-Na,Na+1),range(-Nb,Nb+1),range(-Nc,Nc+1)):
    r1 = np.dot([a,b,c],Mat)+L/2;

    if(np.linalg.norm(r1-L/2)<R):
        for u in range(10):
            wave[u,a,b,c] = wl[u][rec];
        rec += 1;
    else:
        for i in range(n):
            ri = np.array([pos[i].x,pos[i].y,pos[i].z]);
            rij= np.linalg.norm(r1-ri);
            if(rij>0.01):
                for u in range(10):
                    c1 = np.exp(-kappa[u]*rij)/rij*be;
                    wave[u,a,b,c] += c1*wl[u][i];
    for u in range(10):
        for v in range(10):
            M[u,v] += wave[u,a,b,c]*wave[v,a,b,c]*(r1-L/2);
    if(c==0 and b==0):
        print(str(a))

rmn = [[np.linalg.norm(M[u,v]/np.linalg.norm(wave[u])/np.linalg.norm(wave[v])) for v in range(10)] for u in range(10)];
print(rmn)

mH = 1.0078;
mLi = 6.9410;
M_sphere = n*d**3/4**3*(mH+mLi)*4*1.66054*10**-27;
#dE = (kappa_s**2-kappa_p**2)*24.046*10**3
#alpha_s = 1/np.linalg.norm(s);
#alpha_p = 1/np.linalg.norm(p);

#out = {'rmn':np.linalg.norm(res),'dE':dE,
#               'M_sphere':M_sphere,'as':alpha_s,'ap':alpha_p};
#import json;

#with open('rabi.json','w') as file:
#    json.dump(out,file);


#rmn, M = out['rmn'],out['M_sphere'];
#cs, cp, P = [1],[0],[1];

#R = 400;
V = 1; # voltage: V
q = 4*np.pi*8.85419*10**-12*V*4*R*10**-10; #C
E = 1*10**5; #V/m
mn = 1.674927*10**-27;
k = 1.38065*10**-23;
hbar = 1.0545716*10**-34;
Omega = [[q*E/M_sphere*mn/hbar*rmn[u][v]*10**-10*10**-6 for v in range(10)] for u in range(10)];

import json;

with open('omega.json','w') as file:
    json.dump(Omega,file);