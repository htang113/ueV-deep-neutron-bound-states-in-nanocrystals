# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:55:40 2022

@author: 17000
"""

import os;
import numpy as np;
import shutil;
import time;
import json;
from pymatgen import core;
from TD_transition import vasp_task, ab_initio_scheme, analysis;

scheme = ab_initio_scheme('NV','submit.sh');
scheme.relax();
scheme.dfpt();
scheme.phonopy();

