# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 18:01:18 2022

@author: 17000
"""

import os;
import numpy as np;
import shutil;
import time;
import json;
from pymatgen import core;

class vasp_task(object):
    def __init__(self, name, route, jobs=100,username = 'th1543', sleep=5):
        self.name = name;
        self.jobs = jobs;
        self.username = username;
        self.sleep = sleep;
        self.route = route;
        os.mkdir(self.route+name);
        
    def set_pos(self,source):
        shutil.copyfile(source,self.route+self.name+'/POSCAR');
        
    def set_calculator(self, fl = ['INCAR','KPOINTS','POTCAR','submit.sh']):
        nl = ['INCAR','KPOINTS','POTCAR','submit.sh']
        for i in range(4):
            shutil.copyfile(self.route+fl[i],self.route+self.name+'/'+nl[i]);
    def submit(self):
        while(len(os.popen('squeue -u '+self.username).readlines())>self.jobs-1):
            time.sleep(self.sleep*60);
        os.chdir(self.route+self.name);
        os.system('sbatch submit.sh');
        os.chdir(self.route);
    def wait_complete(self):
        while(len(os.popen('squeue -u '+self.username).readlines())>1):
            time.sleep(self.sleep*60);

class ab_initio_scheme(object):
    def __init__(self,name_mat,sub,jobs = 100):
        self.mat = name_mat;
        self.sub = sub;
        self.route = os.getcwd()+'/';
        self.jobs = jobs;

    def relax(self):
        rel = vasp_task('rel',self.route,jobs=self.jobs);
        rel.set_pos(self.route+'POSCAR');
        rel.set_calculator(['INCAR_rel','KPOINTS','POTCAR',self.sub]);
        rel.submit();
        rel.wait_complete();
    
    def dfpt(self):
        rel = vasp_task('phonon',self.route,jobs = self.jobs);
        rel.set_pos(self.route+'rel/CONTCAR');
        rel.set_calculator(['INCAR_dfpt','KPOINTS','POTCAR',self.sub]);
        rel.submit();
        rel.wait_complete();
    
    def phonopy(self):
        shutil.copyfile(self.route+'band.conf',self.route+'phonon/band.conf');
        os.chdir(self.route+'phonon');
        os.system("phonopy -d --dim='1 1 1'");
        time.sleep(20);
        os.system("phonopy --fc vasprun.xml");
        time.sleep(20);
        os.system("phonopy --dim='1 1 1' -c POSCAR band.conf");
        time.sleep(20);
        os.chdir(self.route);
        
    def read_phonon(self,step=0.1):
        file =  open(self.route+'phonon/band.yaml','r');
        data = file.readlines();
        file.close();
        self.pos = core.IStructure.from_file(self.route+'phonon/POSCAR');
        self.natm = len(self.pos);
        self.step = step;
        i_ind = 0;
        self.freq = [];
        self.Nd  = [];
        while(data[i_ind]!='- q-position: [    0.3330000,    0.6670000,    0.0000000 ]\n'):
            content = data[i_ind];
            if(content[:6] == '  - # '):
                self.freq.append(float(data[i_ind+1].split()[1][:-1])*4.13567);
                self.Nd.append([]);
                for q in range(self.natm):
                    self.Nd[-1].append([]);
                    for p in range(3):
                        self.Nd[-1][-1].append(float(data[i_ind+q*4+p+4].split()[2][:-1]));
            i_ind += 1;
        indlist = [i for i in range(self.natm*3)];
        os.mkdir(self.route+'2nd');
        os.mkdir(self.route+'2nd/config');
        for ind in indlist:
            mode = self.Nd[ind];
            latt = self.pos.lattice.matrix;
            atm = self.pos.atomic_numbers;
            dc = {-1:'m',1:'p'}
            for ud in [-1,1]:
                crd = [np.array([self.pos[i].a,self.pos[i].b,self.pos[i].c])+self.step*ud*np.array(np.matrix(mode[i])*np.matrix(latt)**-1)[0] for i in range(self.natm)];
                pos1 = core.IStructure(latt,atm,crd);
                pos1.to('poscar','2nd/config/POSCAR_'+str(ind)+'_'+dc[ud]);
    
    def sampling(self):
        self.pos = core.IStructure.from_file(self.route+'phonon/POSCAR');
        self.natm = len(self.pos);
        name = '0';
        rel=vasp_task('2nd/'+name,self.route,jobs = self.jobs);
        rel.set_pos(self.route+'phonon/POSCAR');
        rel.set_calculator(['INCAR_daq','KPOINTS','POTCAR',self.sub]);
        rel.submit();
        for i in range(self.natm*3):	
            for ud in range(2):
                dc = {0:'m',1:'p'};
                name = str(i)+'_'+dc[ud];
                rel=vasp_task('2nd/'+name,self.route);
                rel.set_pos(self.route+'2nd/config/POSCAR_'+name);
                rel.set_calculator(['INCAR_daq','KPOINTS','POTCAR',self.sub]);
                rel.submit();

    def read_daq(self,folder,qtt):
        os.chdir(folder);
        out = {};
        if('D' in qtt):
            D = os.popen('grep D_diag -A 4 OUTCAR').readlines();
            D = [float(d1.split()[0]) for d1 in D[2:]];
            out['D'] = D;
        if('A' in qtt):
            Af = os.popen("grep 'Fermi contact' -A "+str(self.natm+3)+" OUTCAR").readlines();
            Ad = os.popen("grep Dipolar -A "+str(self.natm+3)+" OUTCAR").readlines();
            Af =[float(u.split()[-1]) for u in Af[4:]];
            A = [[(float(Ad[4+u].split()[1+v])+Af[u]*(v<3)) for v in range(6)] for u in range(self.natm)];
            out['A'] = A;
        if('Q' in qtt):
            EFG = os.popen("grep 'Electric field gradients after diagonalization' -A "+str(self.natm+4)+" OUTCAR").readlines();
            Q = [[float(q)*2.418*10**-3 for q in EFG[5:][kl].split()[1:4]] for kl in range(self.natm)];
            out['Q'] = [];
            for q in range(self.natm):
                out['Q'].append([(Q[q][0]-Q[q][1]-Q[q][2]),(-Q[q][0]+Q[q][1]-Q[q][2]),2*Q[q][2]]);

        os.chdir(self.route);
        return out;
    
    def output(self,quantities,step=0.1):
        self.step = step;
        self.pos = core.IStructure.from_file(self.route+'phonon/POSCAR');
        self.natm = len(self.pos);
        ls = os.listdir(self.route+'2nd/');
        nlist = [str(i) for i in range(10)];
        out = {};
        out['step'] = self.step;
        out['natm'] = self.natm;
        out['Nd'] = self.Nd;
        out['freq'] = self.freq;
        for fold in ls:
            if(fold[0] in nlist):
                out[fold] = self.read_daq(self.route+'2nd/'+fold,quantities);
        with open(self.route+'out.json','w') as file:
            json.dump(out,file);

class analysis(object):
    def __init__(self,file):
        with open(file,'r') as f1:
            self.data = json.load(f1);
            self.natm = int(self.data['natm']);
            self.pos = core.IStructure.from_file('POSCAR');
    def set_atm(self,at,zaxis,mass,gyromagnetic,quadrupole,I):
        zaxis = np.array(zaxis)/np.linalg.norm(zaxis);
        ls = {};
        for key in self.data:
            if(key not in ['step', 'natm', 'Nd', 'freq']):
                ls[key] = {};
                dp = self.data[key];
                if('A' in dp):
                    if('gyromagnetic' not in dir(self)):
                        self.gyromagnetic = [];
                        i_idx = 0;
                        for i in range(self.natm):
                            self.gyromagnetic.append(gyromagnetic[i_idx]);
                            if(i!=self.natm-1 and self.pos.atomic_numbers[i]!=self.pos.atomic_numbers[i+1]):
                                i_idx += 1;
                    atm = dp['A'][at-1];
                    mat = [[atm[0],atm[3],atm[4]],[atm[3],atm[1],atm[5]],[atm[4],atm[5],atm[2]]]
                    res = np.dot(zaxis*np.matrix(mat),zaxis)[0,0];
                    ls[key]['A'] = res*self.gyromagnetic[at-1];
                if('Q' in dp):
                    if('quadrupole' not in dir(self)):
                        self.quadrupole = [];
                        self.I = [];
                        i_idx = 0;
                        for i in range(self.natm):
                            self.quadrupole.append(quadrupole[i_idx]);
                            self.I.append(I[i_idx]);
                            if(i!=self.natm-1 and self.pos.atomic_numbers[i]!=self.pos.atomic_numbers[i+1]):
                                i_idx += 1;
                    Q = dp['Q'][at-1][2] #*16.6*2.418*10**-3*(3/4); #[P. S. PREGOSIN, E. W. RANDALL, and A. I. WHITE, Chem. Comm., 1602 (1971).]
                    ls[key]['Q'] = Q*self.quadrupole[at-1]/4/self.I[at-1]/(2*self.I[at-1]-1)*3/2;
                if('D' in dp):
                    D = dp['D'][2]*3/2/1000;
                    ls[key]['D'] = D;
        
        opt = {};
        for i in range(3*self.data['natm']):
            opt[i] = {};
            for qtt in list(self.data['0'].keys()):
                res = ls[str(i)+'_m'][qtt]+ls[str(i)+'_p'][qtt] - 2*ls['0'][qtt];
                opt[i][qtt] = res/self.data['step']**2;
        self.opt = opt;
        self.mass = [];
        i_idx = 0;
        for i in range(self.natm):
            self.mass.append(mass[i_idx]);
            if(i!=self.natm-1 and self.pos.atomic_numbers[i]!=self.pos.atomic_numbers[i+1]):
                i_idx += 1;
        self.M = [];
        for ind in range(3*self.natm):
            mode = self.data['Nd'][ind];
            self.M.append(sum([np.linalg.norm(mode[i])**2*self.mass[i] for i in range(len(mode))]));

    def f(self, qtt, T):
        coef = 4.18;
        res = 0;
        for i in range(3,3*self.natm):
            res += coef*self.opt[i][qtt]/(2*self.data['freq'][i]*self.M[i])/(np.exp(self.data['freq'][i]/T/0.08617)-1);
        return res;
    
    def df(self, qtt, T):
        coef = 4.18;
        res = 0;
        for i in range(3,3*self.natm):
            res += coef*self.opt[i][qtt]*np.exp(self.data['freq'][i]/T/0.08617)*self.data['freq'][i]/T**2/0.08617/(2*self.data['freq'][i]*self.M[i])/(np.exp(self.data['freq'][i]/T/0.08617)-1)**2;
        return res;
    
    def calculate(self, quantities, Tlist= [i*1 for i in range(4,1000)]):
        if('A' in quantities):
            self.Al = [self.f('A',T)*10**3 for T in Tlist];
            self.dA = [self.df('A',T)*10**6 for T in Tlist];
        if('D' in quantities):       
            self.Dl = [self.f('D',T)*10**3 for T in Tlist];
            self.dD = [self.df('D',T)*10**6 for T in Tlist];
        if('Q' in quantities):
            self.Ql = [self.f('Q',T)*10**3 for T in Tlist];
            self.dQ = [self.df('Q',T)*10**6 for T in Tlist];      
        
        with open('td_daq.txt','w') as file:
            string = str('T (K)\t');
            if('A' in quantities):
                string += 'Delta A (kHz)\t dA/dT (Hz/K)\t';
            if('D' in quantities):
                string += 'Delta D (MHz)\t dD/dT (kHz/K)\t';
            if('Q' in quantities):
                string += 'Delta Q (kHz)\t dQ/dT (Hz/K)\t';
            string = string[:-1]+'\n';
            file.write(string);
                
            for i in range(len(Tlist)):
                string = str(Tlist[i])+'\t';
                if('A' in quantities):
                    string += str(self.Al[i])+'\t'+str(self.dA[i])+'\t';
                if('D' in quantities):
                    string += str(self.Dl[i])+'\t'+str(self.dD[i])+'\t';
                if('Q' in quantities):
                    string += str(self.Ql[i])+'\t'+str(self.dQ[i])+'\t';
                string = string[:-1]+'\n';
                file.write(string);
            
