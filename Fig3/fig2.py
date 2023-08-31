# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 14:06:48 2022

@author: 17000
"""

import json;
import numpy as np;

with open('H_data.json','r') as file:
    data = json.load(file);
    
b_dic = {'H':-3.739011+2/np.sqrt(3/4)*(-1/4)*25.2749,
         'He':3.263,
         'Li':-2.22,
         'Be':7.79,
         'B':6.654,
         'C':6.646,
         'N':9.36,
         'O':0.583,
         'F':5.654,
         'Ne':4.566,
         'Na':3.632,
         'Mg':5.375,
         'Al':3.449,
         'Si':4.1491,
         'P':5.13,
         'S':2.847,
         'Cl':3.08,
         'Ar':1.909,
         'K':3.67,
         'Ca':4.7,
         'Sc':12.29,
         'Ti':-3.438,
         'V':-0.3824,
         'Cr':3.6357,
         'Mn':-3.732,
         'Fe':9.45,
         'Co':2.49,
         'Ni':10.31,
         'Cu':7.718,
         'Zn':5.68,
         'Ga':7.288,
         'Ge':8.185,
         'As':6.58,
         'Se':7.483,
         'Br':6.795,
         'Kr':8.12,
         'Rb':7.09,
         'Sr':7.02,
         'Y':7.75,
         'Zr':7.16,
         'Nb':7.054,
         'Mo':6.715,
         'Tc':6.83,
         'Ru':7.03,
         'Rh':5.88,
         'Pd':5.91,
         'Ag':5.923,
         'Cd':4.875,
         'In':4.065,
         'Sn':5.225,
         'Sb':5.573,
         'Te':5.803,
         'I':5.282,
         'Xe':4.923,
         'Ba':5.07,
         'Cs':5.42
         }

a_dic = {'H' : 0.33267,
         'He': 0,
         'Li': 0.04543,  
         'Be':0.007,
         'B':0.0055,
        'C':0.0035,
        'N':1.9,
        'O':0.000192,
        'F':0.00965,
        'Ne':0.0394,
        'Na':0.5305,
        'Mg':0.0633,
        'Al':0.2313,
        'Si':0.1713,
        'P':0.1726,
        'S':0.531,
        'Cl':0.4336,
        'Ar':0.6759,
        'K':2.11,
        'Ca':0.432,
        'Sc':27.52,
        'Ti':6.0913,
        'V':5.084,
        'Cr':3.058,
        'Mn':13.32,
        'Fe':2.563,
        'Co':37.186,
        'Ni':4.4916,
        'Cu':3.782,
        'Zn':1.112,
        'Ga':2.75,
        'Ge':2.2,
        'As':4.51,
        'Se':0.615,
        'Br':6.92,
        'Kr':0.0032,
        'Rb':0.381,
        'Sr':1.286,
        'Y':1.282,
        'Zr':0.1853,
        'Nb':1.155,
        'Mo':2.484,
        'Tc':20.1,
        'Ru':2.5613,
        'Rh':144.87,
        'Pd':6.94,
        'Ag':63.34,
        'Cd':2520,
        'In':193.8,
        'Sn':0.6269,
        'Sb':4.915,
        'Te':4.71,
        'I':6.156,
        'Xe':23.9,
        'Ba':1.1,
        'Cs':29};


h = 6.62607015*10**-34; #international
mn = 1.67492749804*10**-27;  #kg
k = 1.3806504*10**-23; #C
e = 1.602176*10**-19;
coef = h**2/(2*np.pi*mn)*10**-15*10**30/e*10**6;  #eV*A^3
c_a = 220004.9571;

T = [];
E = [];
Name = [];

T1 = [];
E1 = [];
Name1 = [];

for dp in data:
    n_ele = len(dp['Ele']);
    if(not sum([(dp['Ele'][i] not in b_dic.keys()) for i in range(n_ele)])):   
        if(dp['stable']):
            b_sum = sum([b_dic[dp['Ele'][i]]*dp['N'][i] for i in range(n_ele)]);
            b_ave = -coef*b_sum/dp['V'];
            
            a_sum = sum([a_dic[dp['Ele'][i]]*dp['N'][i] for i in range(n_ele)]);
            tau   = (a_sum/dp['V'])**-1/c_a*10**3;
            if(b_ave*tau>0):
                T.append(tau);
                E.append(b_ave);
                Name.append(dp['name']);

for i in range(len(T)):
    add = True;
    rem = [];
    for j in range(len(T1)):
        if(T1[j]>T[i] and E1[j]>E[i]):
            add = False;
        if(T1[j]<T[i] and E1[j]<E[i]):
            rem.append(j);
    for r in rem[::-1]:
        del T1[r];
        del E1[r];
        del Name1[r]
    if(add):
        T1.append(T[i]);
        E1.append(E[i]);
        Name1.append(Name[i]);
            
import matplotlib.pyplot as plt;
SMALL_SIZE = 10
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.figure(figsize=[12,8]) 
plt.scatter(T,E)
plt.scatter(T1,E1)
plt.axis([0,0.7,0,0.5])
for i, label in enumerate(Name1):
    plt.text(T1[i], E1[i], label, ha='center', va='bottom');

# =============================================================================
# with open('zoomin.txt','w') as file:
#     for i in range(len(T)):
#         file.write(str(T[i])+'\t'+str(E[i]));
#         if(i<len(T1)):
#             file.write('\t'+str(T1[i])+'\t'+str(E1[i])+'\t'+Name1[i]);
#         file.write('\n')
# =============================================================================


with open('zoomin.txt','w') as file:
    for i in range(len(T)):
        if(T[i]*E[i]>0.055):
            file.write(str(T[i])+'\t'+str(E[i])+'\t'+Name[i]+'\n');