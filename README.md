
[![Cross-Platform Tests](https://github.com/HUISUN24/MaterialsMap/actions/workflows/lastCommit.yml/badge.svg)](https://github.com/HUISUN24/MaterialsMap/actions/workflows/lastCommit.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/materialsmap?label=PyPI&color=green)](https://pypi.org/project/materialsmap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**2024 Paper:**
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.mtla.2024.102153-blue)](https://doi.org/10.1016/j.mtla.2024.102153)
[![Arxiv](https://img.shields.io/badge/arXiv-2403.19084-8F1515?style=flat&logo=arxiv&logoColor=red)](https://doi.org/10.48550/arXiv.2403.19084)

`MaterialsMap` is Python package for mapping properties, manufacturing feasibility, and desirability. We focus on guiding materials design graphically while proving an API to underlying methods, so that others can utlize it as an engine behind their tools, like machine learning (ML) based alloy design.

It can utilize both free and open source (FOSS) [pycalphad](http://pycalphad.org) and industry-trusted [Thermo-Calc](https://thermocalc.com) for performing underlying thermodynamic calculations.

## Theory

 uses several different methods out-of-the-box. These include thermodynamic equlibrium calculations, , and 5 different models for predicting cracking susceptibility.

Our upcoming (March 2024) article discusses in detail all implemented methods, underlying thermodynamics, and their applications.

`MaterialsMap` uses several different methods out-of-the-box. These include thermodynamic equilibrium calculations, [Scheil-Gulliver solidification](https://en.wikipedia.org/wiki/Scheil_equation), and 5 different models for predicting cracking susceptibility. The thermodynamic equilibrium calculations predict phases formed at specific conditions by minimizing the Gibbs energy of the system and are typically closer to experimnental observations of slow-colled materials, while the Scheil-Gulliver simulations capture non-equilibrium rapid solidification effects (e.g. in welding or additive manufacturing) by assuming no diffusion in solid phases, equilibrium at the solid/liquid interface, and complete mixing in the liquid phase. `MaterialsMap` combines these approaches to provide comprehensive phase formation predictions at the two extrema of cooling rates to determine feasibility. The five implemented crack susceptibility criteria include Freezing Range (FR), Crack Susceptibility Coefficient (CSC), Kou Criteria, Improved Crack Susceptibility Coefficient (iCSC), and Simplified Rappaz-Drezet-Gramaud (sRDG), enabling users to assess hot cracking risks from multiple perspectives / mechanisms. 

Please refer to [our 2024 Materialia article](https://doi.org/10.1016/j.mtla.2024.102153) which discusses in detail all implemented methods, underlying thermodynamics, and their applications.

## Use Example (SS304L-NiCr-V)

![results of SS304L-NiCr-V](https://raw.githubusercontent.com/HUISUN24/MaterialsMap/main/docs/demo-results.png)

```python
import os
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from materialsmap.core.compositions import generateCompositions, createComposition
from materialsmap.ref_data import periodic_table, materials
from materialsmap.core.pycalphad_run import pycalphad_eq, pycalphad_scheil
from materialsmap.core.GenerateEqScript import createEqScript
from materialsmap.core.ReadEqResult import getEqdata
from materialsmap.core.GenerateScheilScript import createScheilScript
from materialsmap.core.ReadScheilResult import getScheilSolidPhase
from materialsmap.plot.FeasibilityMap import plotMaps

# Create Compositions
comps = ['SS304L', 'NiCr', 'V']
eleAmountType = 'massFraction'
pressure = 101325
ngridpts = 41  # number of points along each dimension of the composition grid
TemperatureRange = (900, 2300, 10)  #(lower limit, upper limit, temperature step)
indep_comps = [comps[1], comps[2]]  # choose them automatically
for i in comps:
    if i in periodic_table:
        materials[i] = {i: 1}
    elif i not in materials.keys():
        materials['SS304L'] = {'Ni': 0.09611451943, 'Cr': 0.1993865031,
                               'Fe': 0.7044989775}  # the composition of this element/alloys(in weight fractions)
maxNumSim = 250  # maximum number of simulations in each TCM file

# Equilibrium simulation settings
pressure = 101325
database = 'TCFE8'  # <userDatabase>.TDB or TCFE8
eleAmountType = 'massFraction'  # Candidates: massFraction massPercent moleFraction molePercent
output_Eq = f'{TemperatureRange[0]}-{TemperatureRange[1]}-{TemperatureRange[2]}-{comps[0]}-{comps[1]}-{comps[2]}-Eq'

# Create folder in curent path to store simulation results
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

# Save compostion results
compositions_list = generateCompositions(indep_comps, ngridpts)
Compositions, numPoint, comp, numSimultion = createComposition(indep_comps, comps, compositions_list, materials, path)
settings = [TemperatureRange, numPoint, numSimultion, comp, comps, indep_comps, os.path.abspath(database), pressure, eleAmountType]
np.save(f'{path}/setting.npy', settings)

# Running with PyCalphad
pycalphad_eq(path)
pycalphad_scheil(path, 2000)  # temperature to start scheil if not eq results

# Running with Thermo_Calc
# Create TCM files with path
createEqScript(path)
createScheilScript(path, 2000)  # temperature to start scheil if not eq results

# Open TCM files with Thermo_Calc

# Collect results from Thermo_Calc
getEqdata(path)
getScheilSolidPhase(path)

# Plot deleterious phase diagram and crack susceptibility map 
plotMaps(path, 'pycalphad')
```

## Installation

### PyPI (recommended)

`MaterialsMap` can be quickly installed from PyPI with a simple:
```cmd
pip install materialsmap
```

### Development Versions

To install an editable development version with pip:

```cmd
git clone https://github.com/HUISUN24/materialsmap.git
cd materialsmap
pip install -e .
```
    
Upgrading scheil later requires you to run ``git pull`` in this directory.

## Testing

Automated testing is performed on every commit to the repository, as defined in `.github/workflows/lastCommit.yml` workflow. On your system, you can also run it with a simple:
```cmd
pytest
```

