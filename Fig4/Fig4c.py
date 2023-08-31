# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 15:40:30 2022

@author: 17000
"""
import json;
import numpy as np;

with open('omega.json','r') as file:
    data = json.load(file);

El = np.linspace(0,5,10);

sp = np.mean(data[0][1:4]);

with open('rmn.json','r') as file:
    rmn = json.load(file);

R = 400;
mH = 1.0078;
mLi = 6.9410;
M_sphere = 2.185276*10**-19;

V = 1; # voltage: V
q = 4*np.pi*8.85419*10**-12*V*4*R*10**-10; #C
E = 1*10**5; #V/m
mn = 1.674927*10**-27;
k = 1.38065*10**-23;
hbar = 1.0545716*10**-34;
Omega = [[q*E/M_sphere*mn/hbar*rmn[u][v][0]*10**-10*10**-6 for v in range(10)] for u in range(10)];

sp = np.linalg.svd(np.array(Omega)[1:4,0:1])[1][0];
pe = np.linalg.svd(np.array(Omega)[1:4,4:6])[1][0];
pt = np.linalg.svd(np.array(Omega)[1:4,6:9])[1][0];
ps = np.linalg.svd(np.array(Omega)[1:4,9:10])[1][0];

with open('E_d.txt','w') as file:
    for i in range(len(El)):
        file.write(str(El[i])+'\t');
        file.write(str(sp*El[i]/2/np.pi)+'\t');
        file.write(str(pe*El[i]/2/np.pi)+'\t');
        file.write(str(pt*El[i]/2/np.pi)+'\t');
        file.write(str(ps*El[i]/2/np.pi)+'\n');

