import numpy as np
import pandas as pd
from pycalphad import Database, equilibrium, Model, variables as v
from materialsmap.ref_data import eleweight
from collections import defaultdict
import requests
import json
import os

def ML_run(path,properties='melting_temperature'):
    """_summary_

    Args:
        path (str): path to open the setting and store the results
        properties (str): properteis of interested. Defaults to melting_temperature.
    """
    setting = np.load(f'{path}/setting.npy',allow_pickle=True)
    comps = []
    for i in setting[3]:
        comps.append(i.upper())
    comps.append('VA')
    compositions_list = []
    indepcomp_list = []
    all_comp_tot = []
    df = pd.read_excel(f'{path}/composition_for_feasibilityMap.xlsx')
    data_str = "[" # change the compositions to MAPP format
    for n in range(len(df.index)):
        comp_list={}
        all_comp = []
        nocomp_list = []
        indep = 0
        for o,i in enumerate(setting[3]):
            if indep == 0 and float(df.loc[n,[i]].values) > 0:
                indep  = 1
                indepcomp_list.append(i)
                all_comp.append(i)
                continue;
            elif float(df.loc[n,[i]].values) < 1E-5:
                nocomp_list.append(i)
                continue
            i = i.upper()
            all_comp.append(i)
            comp_list[v.W(i)]=float(df.loc[n,[i]].values)
        comp_list_mole = v.get_mole_fractions(comp_list,v.Species(indepcomp_list[n]),eleweight)
        if len(data_str) > 2: 
            data_str += ","
        if len(comp_list_mole) == 0:
            if len(all_comp[0]) > 1:
                comp_write = all_comp[0][0]+all_comp[0][1].lower()
            else:
                comp_write = all_comp[0]
            data_str += '{"9":"' + comp_write + '"}'
        else:
            val_dep = 1
            data_str += '{"9":"'
            for co,val in comp_list_mole.items():
                val_dep = val_dep - val
                cop = str(co)
                comp_write = cop.replace('X_','')
                if len(comp_write) > 1:
                    comp_write = comp_write[0]+comp_write[1].lower()
                else:
                    comp_write = comp_write        
                data_str += comp_write + str(val)
            if len(all_comp[0]) > 1:
                comp_write_dep = all_comp[0][0]+all_comp[0][1].lower()
            else:
                comp_write_dep = all_comp[0]        
            data_str += comp_write_dep + str(val_dep)
            data_str += '"}'   
    data_str += "]"
    # run with MAPP
    url = 'http://206.207.50.58:5007/MT_ML_Qijun_Hong_Predict_noNN'
    payload = data_str
    headers = {'content-type': 'application/json'}
    r = requests.post(url, data=payload, headers=headers)
    # results reformat and store
    df_mp = pd.DataFrame(str(r.content).split('{"melting temperature\": ')[1:], columns=['mp'])
    df_mp['melting_temperature_in_kelvin'] = df_mp.apply(lambda x: float( x['mp'].split(",")[0] ), axis=1)
    df_mp['standard_error_in_kelvin'] = df_mp.apply(lambda x: float( x['mp'].split(":")[1].split("}")[0] ), axis=1)
    data_ML = defaultdict(dict)
    for i in df_mp.index:
        data_ML['Point'+str(i)]['TK'] = df_mp['melting_temperature_in_kelvin'][i]
        data_ML['Point'+str(i)]['STD'] = df_mp['standard_error_in_kelvin'][i]
    isExist = os.path.exists(path+'/ML/'+properties+'/Result/')
    if not isExist:
        os.makedirs(path+'/ML/'+properties+'/Result/')
        print("The new directory is created!")
    output_File = json.dumps(data_ML)
    f = open(path+'/ML/'+properties+'/Result/data_ML.json','w')
    f.write(output_File)
    f.close()