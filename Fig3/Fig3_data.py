# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 09:00:35 2022

@author: 17000
"""
from mp_api.client import MPRester
import json;

with MPRester("2aZd3vZ4Ej86NEK5L14gH14WCCOLDFmV") as mpr:
    docs = mpr.summary.search(elements=["H"], 
                              fields=["theoretical","is_stable","formula_pretty","volume","composition","total_magnetization"])
    
    out = [];
    for doc in docs:
        composition = doc.composition.as_dict();
        elements = list(composition.keys());
        numbers = list(composition.values());
        volume = doc.volume;
        theoretical = doc.theoretical;
        is_stable = doc.is_stable;
        out.append({'V':volume, 'Ele':elements,'N':numbers,
                    'M':doc.total_magnetization,'name':doc.formula_pretty,
                    'theo':theoretical,'stable':is_stable});
    with open('H_data.json','w') as file:
        json.dump(out,file);
        
#    with open('hydride1.json','w') as file:
#        json.dump([doc.dict() for doc in docs],file);
    
    


