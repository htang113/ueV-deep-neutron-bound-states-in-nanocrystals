# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 10:31:02 2021

@author: 17000
"""
import tkinter as tk;
from tkinter import ttk;
from tkinter import *;
import pymatgen as pmg;
import os;
import shutil;

top = Tk();
top.geometry('1500x800+100+100')

tk.Label(top, text='TD spin-transition', font = ('Times New Roman',20)).grid(row = 0);

def set_sel(name, pos, values,lab=True):
    functional = ttk.Combobox(top,font = ('Times New Roman',20));
    functional['values'] = values;
    functional.current(0);
    functional.grid(row=pos[0]+1,column=pos[1])
    if(lab):
        tk.Label(top,text = name, font= ('Times New Roman',20)).grid(row = pos[0],column=pos[1]);
    return functional;

def set_fil(name,pos, default,lab=True):
    if(lab):
        tk.Label(top, text=name, font = ('Times New Roman',20)).grid(row = pos[0],column=pos[1]);
    A = Entry(top, font = ('Times New Roman',20));
    A.insert(0,default)
    A.grid(row=pos[0]+1,column=pos[1]);
    return A;

def set_chk(name,pos):
    var1 = tk.BooleanVar()
    c1 = tk.Checkbutton(top, text=name,font = ('Times New Roman',20), 
                        variable=var1, onvalue=True, offvalue=False)
    c1.grid(row=pos[0],column=pos[1]);
    return var1;

pos = pmg.core.IStructure.from_file('POSCAR');
natm = len(pos);
El = [];
for sp in pos.species:
    if(sp not in El):
        El.append(sp);
    
c = {
    'functional':set_sel('functional', (1,0),
                   ('GGA-PBE','mGGA-SCAN','HSE06')),
    'PREC': set_sel('precision', (3,0),
                   ('Accurate','Normal','high','medium','low')),
    'charge' :set_sel('charge (e)', (5,0),
                      ('-1','0','1','-2','2','3','-3')),
    'kpoint' :set_fil('k-point (Nx Ny Nz)',(1,1),'3 3 3'),
    'POT_dir': set_fil('POTCAR dir', (3,1),'potpaw_PBE.54'),
    'prefix': set_fil('prefix', (5,1),'task'),
    'NPAR' :set_fil('NPAR', (1,2), '8'),
    'LDIPOL' : set_chk('LDIPOL',(12,0)),
    'step'   :set_fil('step dq_i (A)',(3,2),'0.1'),
    'quantities' :set_fil('quantities to compute',(5,2),'D A Q'),
    'atom': set_fil('atom (nuclear spin)',(12,2),'127'),
    'zaxis':set_fil('z axis (for Azz)', (12,1),'0.816 0.471 0.333'),
    'T_start':set_fil('start temperature (K)',(14,0),'4'),
    'T_end':set_fil('end temperature (K)',(14,1),'1000'),
    'N_T':set_fil('No. of Temperature',(14,2),'997')
    };

tk.Label(top, text='I', font = ('Times New Roman',20)).grid(row = 8);
tk.Label(top, text='mass (u)', font = ('Times New Roman',20)).grid(row = 9);
tk.Label(top, text='gyromagnetic (MHz/T)', font = ('Times New Roman',20)).grid(row = 10);
tk.Label(top, text='quadrupole (mb)', font = ('Times New Roman',20)).grid(row = 11);

for i in range(len(El)):
    string = str(list(El)[i]);
    c['I_'+string] = set_fil(string,(7,i+1),'0.5');
    c['m_'+string] = set_fil('',(8,i+1),'12');
    c['g_'+string] = set_fil('',(9,i+1),'3.077');
    c['q_'+string] = set_fil('',(10,i+1),'16.6');
    

def generate():
########## modification ########################
    pre = c['prefix'].get();
    if(not os.path.exists(pre)):
        os.mkdir(pre);
    with open(pre+'/KPOINTS','w') as file:
        file.write('Automatic mesh\n 0\n')
        file.write('Gamma\n'+c['kpoint'].get()+'\n')
        file.write('0. 0. 0.  ');
    shutil.copyfile('file/INCAR',pre+'/INCAR_rel');
    shutil.copyfile('file/INCAR',pre+'/INCAR_dfpt');
    shutil.copyfile('file/INCAR',pre+'/INCAR_daq');
    shutil.copyfile('file/TD_transition.py',pre+'/TD_transition.py');
    potdir = c['POT_dir'].get();
    ne = [];
    lines = [];
    for spc in list(El):
        if(os.path.exists(potdir+'/'+str(spc)+'/POTCAR')):
            nm = potdir+'/'+str(spc)+'/POTCAR';
        else:
            nm = potdir+'/'+str(spc)+'/'+str(spc);
        with open(nm,'r') as file:
            qt = file.readlines()[:];
            ne.append(float(qt[1][:-1]));
            lines += qt;
    with open(pre+'/POTCAR','w') as file:
        file.writelines(lines);
    
    i_idx = 0;
    nelect = 0;
    for i in range(natm):
        nelect += ne[i_idx];
        if(i!=natm-1 and pos.atomic_numbers[i]!=pos.atomic_numbers[i+1]):
            i_idx += 1;
    nelect -= int(c['charge'].get());
    with open(pre+'/INCAR_rel','a') as file:
        file.write('PREC      = '+c['PREC'].get()+'\n');
        file.write('EDIFF     = 1E-6\n');
        file.write('NELECT    = '+str(nelect)+'\n');
        file.write('ALGO      = Fast\n');
        file.write('LREAL     = Auto\n');
        file.write('ISYM      = 2\n');
        file.write('EDIFFG    = -0.01\n');
        file.write('IBRION    = 2\n');
        file.write('NSW       = 100\n');
        file.write('POTIM     = 0.05\n');
        file.write('NFREE     = 10\n');
        file.write('NPAR      ='+c['NPAR'].get()+'\n');
    with open(pre+'/INCAR_dfpt','a') as file:
        file.write('PREC      = '+c['PREC'].get()+'\n');
        file.write('EDIFF     = 1E-8\n');
        file.write('NELECT    = '+str(nelect)+'\n');
        file.write('ALGO      = normal\n');
        file.write('LREAL     = False\n');
        file.write('ISYM      = 2\n');
        file.write('IBRION    = 8\n');
    with open(pre+'/INCAR_daq','a') as file:
        file.write('PREC      = '+c['PREC'].get()+'\n');
        file.write('EDIFF     = 1E-8\n');
        if('D' in c['quantities'].get()):
            file.write('LDMATRIX  = T\n');
        if('A' in c['quantities'].get()):
            file.write('LHYPERFINE = T\n');
        if('Q' in c['quantities'].get()):
            file.write('LEFG = T\n');
        if(c['functional'].get()=='mGGA-SCAN'):
            file.write('METAGGA = SCAN\n');
        elif(c['functional'].get()=='HSE06'):
            file.write('LHFCALC = T\n HFSCREEN = 0.2\n TIME   = 0.4\n');
        file.write('NELECT    = '+str(nelect)+'\n');
        file.write('ALGO      = normal\n');
        file.write('LREAL     = False\n');
        file.write('NSW      = 0\n');
        file.write('IBRION    = -1\n');
        if(c['LDIPOL'].get()):
            file.write('LDIPOL = T\n');
        if(c['functional'].get()=='HSE06' or ('D' in c['quantities'].get())):
            file.write('ISYM      = 3\n');
        else:
            file.write('ISYM      = 2\n');
        file.write('NPAR      ='+c['NPAR'].get()+'\n');
        
    shutil.copyfile('file/submit.sh',pre+'/submit.sh');
    shutil.copyfile('file/band.conf',pre+'/band.conf');
    shutil.copyfile('POSCAR',pre+'/POSCAR');
    shutil.copyfile('file/analysis.py',pre+'/analysis.py');
    shutil.copyfile('file/run.py',pre+'/run.py');
    with open(pre+'/run.py','a') as file:
        file.write('scheme.read_phonon(step='+c['step'].get()+');\n')
        file.write('scheme.sampling();\n')
        file.write("scheme.output("+str(c['quantities'].get().split())+",step="+c['step'].get()+");");
     
    with open(pre+'/analysis.py','a') as file:
        Ilist = str([float(c['I_'+str(q)].get()) for q in list(El)]);
        mlist = str([float(c['m_'+str(q)].get()) for q in list(El)]);
        glist = str([float(c['g_'+str(q)].get()) for q in list(El)]);
        qlist = str([float(c['q_'+str(q)].get()) for q in list(El)]);
        file.write('calc.set_atm('+c['atom'].get()+','+str([float(i1) for i1 in c['zaxis'].get().split()])+','+mlist);
        file.write(",gyromagnetic="+str(glist));
        file.write(",quadrupole="+str(qlist)+",I="+str(Ilist)+");\n");
        nT = c['N_T'].get();
        start = c['T_start'].get();
        coef = str((float(c['T_end'].get())-float(start))/(float(nT)-1));
        file.write("calc.calculate('"+c['quantities'].get()+"',Tlist=[i*"+coef+'+'+start+' for i in range('+nT+')]);\n');
    with open(pre+'/band.conf','a') as file:
        file.write('ATOM_NAME = '+''.join([str(i2)+' ' for i2 in list(El)]))
    top.destroy()
Button(top,text = 'Generate',command = generate, font = ('Times New Roman',20)).grid(row=16,column=0);
Button(top,text = 'Cancel',command = top.destroy, font = ('Times New Roman',20)).grid(row=16,column=2);
top.mainloop();
