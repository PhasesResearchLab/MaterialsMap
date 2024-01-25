import os
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pycalphad import Database, equilibrium, Model, variables as v
from pycalphad.core.calculate import _sample_phase_constitution
from pycalphad.core.errors import DofError
from pycalphad.core.utils import point_sample
from scheil import simulate_scheil_solidification
from multiprocessing import Pool
from collections import defaultdict
from materialsmap.ref_data import eleweight
import materialsmap.core.istarmap
from sys import argv
import tqdm    
import json

def pycalphad_eq(path):
    """running and store equlibrium calculations based on the settings on path

    Args:
        path (str): path to open the setting and store the results

    Returns:
        json: data_more.json that contains all eq results
    """

    setting = np.load(f'{path}/setting.npy',allow_pickle=True)
    dbf = Database(setting[6])
    comps = []
    for i in setting[3]:
        comps.append(i.upper())
    comps.append('VA')
    phases = list(dbf.phases.keys())
    potentials = {v.N: 1, v.T: setting[0], v.P: setting[7]}  # for equilibrium calculations
    compositions_list = []
    indepcomp_list = []
    all_comp_tot = []
    df = pd.read_excel(f'{path}/composition_for_feasibilityMap.xlsx')
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
        if len(comp_list_mole) == 0:
            comp_list_mole[v.X(nocomp_list[0])]=1E-5
            all_comp.append(nocomp_list[0])
        all_comp_tot.append(all_comp)
        compositions_list.append(comp_list_mole)

    iter_args_equilibrium = []
    for num, composition in enumerate(compositions_list):
        print(f"{composition} ({num+1}/{len(compositions_list)})")
        # Equilibrium calculation for feasibility
        
        for key,val in composition.items():
            composition[key] = float("{:.6f}".format(val))    
        conds = {**composition}
        conds.update(potentials)
        comps_new = all_comp_tot[num]
        comps_new.append('VA')
        iter_args_equilibrium.append((dbf, comps_new, phases, conds))
    # Multiprocessing step:
    cores = os.cpu_count() - 1
    eq_results = []
    with Pool(cores) as p:
        for eq_result in tqdm.tqdm(p.istarmap(equilibrium, iter_args_equilibrium),
                        total=len(iter_args_equilibrium)):
            eq_results.append(eq_result)
            pass
    # Chang eq result to dict
    equilibrium_result = defaultdict(dict)
    for num,eq in enumerate(eq_results):
        equilibrium_result['Point'+str(num)]['TK'] = list(eq.T.values)
        eq_phases_name = set(eq.Phase.values.flatten().tolist()) - {''}
        eq_phases = eq.Phase.values.squeeze().tolist()
        for eq_phase in eq_phases_name:
            equilibrium_result['Point'+str(num)][eq_phase] = []
            for n,va in enumerate(eq["NP"].values.squeeze()):
                count = eq_phases[n].count(eq_phase)
                if count == 0:
                    equilibrium_result['Point'+str(num)][eq_phase].append(0)
                elif count == 1:
                    index = eq_phases[n].index(eq_phase)
                    equilibrium_result['Point'+str(num)][eq_phase].append("{:.8f}".format(float(va[index])))
                elif count == 2:
                    index = [i for i, x in enumerate(eq_phases[n]) if x == eq_phase]
                    result = 0
                    for i in index:
                        result += float(va[i])
                    equilibrium_result['Point'+str(num)][eq_phase].append("{:.8f}".format(result))
                else:
                    print('error')
    isExist = os.path.exists(path+'/Pycalphad/Equilibrium Simulation/Result/')
    if not isExist:
        os.makedirs(path+'/Pycalphad/Equilibrium Simulation/Result/')
        print("The new directory is created!")
    output_File = json.dumps(equilibrium_result)
    f = open(path+'/Pycalphad/Equilibrium Simulation'+'/Result/data_mole.json','w')
    f.write(output_File)
    f.close()

def pycalphad_scheil(path,intial_temperature,liquid_name='LIQUID',step_temperature=1.0,eq_kwargs=None,stop=0.0001, verbose=False, adaptive=True):
    """running and store scheil simulations based on the settings on path

    Args:
        path (str): path to open the setting and store the results
        intial_temperature (_type_): the temperature of scheil
        liquid_name (str, optional): Defaults to 'LIQUID'.
        step_temperature (float, optional): Defaults to 1.0.
        eq_kwargs (_type_, optional): keywords pass to eq calculations. Defaults to None.
        stop (float, optional): the liqud fraction to stop scheil simulations. Defaults to 0.0001.
        verbose (bool, optional): Defaults to False.
        adaptive (bool, optional): Dynamic zoom in when close to stop. Defaults to True.
    """                   

    setting = np.load(f'{path}/setting.npy',allow_pickle=True)
    dbf = Database(setting[6])
    comps = []
    for i in setting[3]:
        comps.append(i.upper())
    comps.append('VA')
    phases = list(dbf.phases.keys())
    compositions_list = []
    indepcomp_list = []
    all_comp_tot = []
    df = pd.read_excel(f'{path}/composition_for_feasibilityMap.xlsx')
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
        if len(comp_list_mole) == 0:
            comp_list_mole[v.X(nocomp_list[0])]=1E-5
            all_comp.append(nocomp_list[0])
        all_comp_tot.append(all_comp)
        compositions_list.append(comp_list_mole)
    LiquidusTemp = []
    isExist = os.path.exists(path+'/Pycalphad/Equilibrium Simulation/Result/')
    if isExist:
        f = open(path+'/Pycalphad/Equilibrium Simulation'+'/Result/data_mole.json')
        eq_results = json.load(f)
        for j in eq_results.values():
            if len(j.keys()) == 1:
                j = None
                LiquidusTemp.append(intial_temperature)
                continue;
            for n,a in enumerate(j['LIQUID']):
                if n+1 == len(j['LIQUID']):
                    print('cannot find liquid phase in eq results, start with back up temperature')
                    LiquidusTemp.append(intial_temperature)
                elif float(a) < 1 and float(j['LIQUID'][n+1]) == 1:
                    LiquidusTemp.append(j['TK'][n+1])
                    break;
    if len(LiquidusTemp) == 0:
        T_liquid = [intial_temperature]*len(compositions_list)
    else:
        T_liquid = LiquidusTemp
    # Generate points for adaptive Scheil starting points (performance)
    points_dict = {}
    for phase_name in phases:
        try:
            mod = Model(dbf, comps, phase_name)
            points_dict[phase_name] = _sample_phase_constitution(mod, point_sample, True, 2000)
        except DofError:
            pass
    liquid_name = 'LIQUID'
    step_temperature = 1.0
    eq_kwargs = {'adaptive':True, 'eq_kwargs':{'calc_opts': {'points': points_dict}}}
    stop = 0.0001
    verbose = False
    adaptive = True
    # Run simulations
    iter_args_scheil = []
    for num, composition in enumerate(compositions_list):
        print(f"{composition} ({num+1}/{len(compositions_list)})")
        for key,val in composition.items():
            composition[key] = float("{:.6f}".format(val))  
        comps_new = all_comp_tot[num]
        comps_new.append('VA')  
        iter_args_scheil.append((dbf, comps_new, phases, composition, T_liquid[num], step_temperature,liquid_name,eq_kwargs,stop,verbose, adaptive))
    # Multiprocessing step:
    cores = os.cpu_count() - 1
    scheil_results_ori = []
    with Pool(cores) as p:
        for scheil_result_ori in tqdm.tqdm(p.istarmap(simulate_scheil_solidification, iter_args_scheil),
                        total=len(iter_args_scheil)):
            scheil_results_ori.append(scheil_result_ori)
            pass
    # Chang scheil result to dict
    scheil_result = defaultdict(dict)
    for num,i in enumerate(scheil_results_ori):
        scheils = i.to_dict()
        scheil_result['Point'+str(num)]['TK'] = scheils['temperatures']
        for pha,val in i.cum_phase_amounts.items():
            if np.sum(val) > 1E-6:
                scheil_result['Point'+str(num)][pha] = val
        scheil_result['Point'+str(num)]['LIQUID'] = (1.0 - np.array(scheils['fraction_solid'])).tolist()
    
    
    isExist = os.path.exists(path+'/Pycalphad/Scheil Simulation/Result/')
    if not isExist:
        os.makedirs(path+'/Pycalphad/Scheil Simulation/Result/')
        print("The new directory is created!")
    output_File = json.dumps(scheil_result)
    f = open(path+'/Pycalphad/Scheil Simulation'+'/Result/data_mole.json','w')
    f.write(output_File)
    f.close()
    return scheil_results_ori