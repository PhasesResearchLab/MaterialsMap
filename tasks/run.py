from materialsmap.core.pycalphad_run import pycalphad_eq,pycalphad_scheil
import os
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from materialsmap.core.compositions import generateCompositions,createComposition
from materialsmap.ref_data import periodic_table,materials
from materialsmap.core.pycalphad_run import pycalphad_eq,pycalphad_scheil
from materialsmap.core.GenerateEqScript import createEqScript
from materialsmap.core.ReadEqResult import getEqdata
from materialsmap.core.GenerateScheilScript import createScheilScript
from materialsmap.core.ReadScheilResult import getScheilSolidPhase
from materialsmap.plot.FeasibilityMap import plotMaps
from pycalphad import Database, equilibrium, Model, variables as v
from materialsmap.ref_data import eleweight
from pycalphad.core.calculate import _sample_phase_constitution
from pycalphad.core.errors import DofError
from pycalphad.core.utils import point_sample
from scheil import simulate_scheil_solidification
import json
# Create Compositions
comps = ['Ag','Al','Cu']
database = './Ag-Al-Cu.TDB'
eleAmountType = 'massFraction'
pressure = 101325
ngridpts = 41  # number of points along each dimension of the composition grid
TemperatureRange = (300, 1500,10) #(lower limit, upper limit, temperature step)
indep_comps = [comps[1], comps[2]]  # choose them automatically
for i in comps:
    if i in periodic_table:
        materials[i] = {i:1}
    elif i not in materials.keys():
        materials['SS304L'] = {'Ni':0.09611451943, 'Cr':0.1993865031, 'Fe':0.7044989775}         # the composition of this element/alloys(in weight fractions)
maxNumSim = 250  # maximum number of simulations in each TCM file
# Equilibrium simulation settings

 # 'C-Cr-Cu-Fe-Mo-Nb-Ni-03-08.tdb' #'C-Cr-Cu-Fe-Mo-Nb-Ni-10-05.tdb'#'Cr-Fe-Ni-Ti-V_04-05.tdb' # <userDatabase>.TDB or TCFE8
eleAmountType = 'massFraction' # Candidates: massFraction moleFraction
output_Eq = f'{TemperatureRange[0]}-{TemperatureRange[1]}-{TemperatureRange[2]}-{comps[0]}-{comps[1]}-{comps[2]}-Eq'
from datetime import datetime
current_dateTime = datetime.now()
if '.tdb' in database or '.TDB' in database:
    database_name = database.split('/')
    path = f'./Simulation/{datetime.now().strftime("%m-%d-%Y")}-{comps[0]}-{comps[1]}-{comps[2]}-database-{database_name[-1][:-4]}'
else:
    path = f'./Simulation/{datetime.now().strftime("%m-%d-%Y")}-{comps[0]}-{comps[1]}-{comps[2]}-database-{database}'
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
    print("The new directory is created!")
compositions_list = generateCompositions(indep_comps,ngridpts)
Compositions, numPoint, comp, numSimultion = createComposition(indep_comps,comps,compositions_list,materials,path)
settings = [TemperatureRange,numPoint,numSimultion,comp,comps,indep_comps,os.path.abspath(database),pressure,eleAmountType]
print(settings)
np.save(f'{path}/setting.npy',settings)
# path = './Simulation/09-07-2023-Ag-Al-Cu-database-Ag-Al-Cu'
pycalphad_eq(path)
#pycalphad_scheil(path,1500)



