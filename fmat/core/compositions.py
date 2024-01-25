import os
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm
from materialsmap.ref_data import relatedEle

def generateCompositions(indep_comps,ngridpts):
    """Generate composition list based on the independetn components and grid density

    Args:
        indep_comps (list): List of independet components 
        ngridpts (int): Int of the grids

    Returns:
        list: List combinations of independent components
    """
    axisList = [item for item in np.arange(0,1+1/ngridpts, 1/ngridpts)]
    comps = []
    for x in axisList:
        for y in axisList:
            if x+y <= 1:
                comps.append({indep_comps[0]:x, indep_comps[1]:y})
    return comps

def getCompostion(indep_comps,alloyCompostion,comps,materials_update): # get element composition from alloy composition
    """transfer the independent components to the list of independent elements

    Args:
        indep_comps (list): List of independet components 
        alloyCompostion (str): str of alloy part in the components
        comps (list): List combinations of independent components
        materials_update (dict): Dicts that contains the compositions of each compotents

    Returns:
        dict: dict of compostion with key as element
    """
    wtPercent = alloyCompostion
    totalAlloy = 0
    restAlloy = ''
    for alloy in comps:
        if alloy in indep_comps:
            totalAlloy += wtPercent[alloy]
        else:
            restAlloy = alloy
    wtPercent[restAlloy] = 1-totalAlloy
    Composition = {}
    for ele in relatedEle:
        Composition[ele] = 0
    for item in relatedEle:
        for alloy in wtPercent.keys():
            if item in materials_update[alloy].keys():
                Composition[item] += (wtPercent[alloy] * materials_update[alloy][item])
            else:
                Composition[item] += 0
    composition = dict()

    for item in Composition.keys():
        if Composition[item] != 0:
            composition[item.upper()] = Composition[item]
    return composition

def createComposition(indep_comps,comps,compositions_list,materials_update,path): #create compositions
    """Create a compositions list of element based on the input of components

    Args:
        indep_comps (list): List of independet components 
        comps (list): List of all components 
        compositions_list (list): List combinations of independent components
        materials_update (dict): Dicts that contains the compositions of each compotents
        path (str): path to store the results

    Returns:
        set: 
        Compositions(list): List combinations of independent elements
        numPoint(int): Number of calculations in each tcm files
        comp(list): the related element
        newCompositions(int): the length of the generate composition list
    """
    Compositions = dict()
    index = 0
    for num, alloys in enumerate(compositions_list): 
        composition = getCompostion(indep_comps,alloys,comps,materials_update) #alloys updated to 3 alloys in this function
        if index == 0: #initialize Compositions
            comp = []
            for item in comps:
                for item_ele in materials_update[item].keys():
                    comp.append(item_ele)
            comp = set(comp)
            comp = [item.upper() for item in comp]
            eleNum = len(comp)
            for item in comp:
                Compositions[item] = []
            for item in comps:
                Compositions['alloy_'+item] = []
            Compositions['Index'] = []
        for key in alloys.keys():
            Compositions['alloy_'+key].append(alloys[key])
        for ele in comp:
            if ele in composition.keys():
                Compositions[ele].append(composition[ele])
            else:
                Compositions[ele].append(0)
        Compositions['Index'].append(index)
        index += 1
    newCompositions = pd.DataFrame(Compositions)
    print(f'Equilibrium simulation #: {len(newCompositions)}')
    numPoint =  len(newCompositions) 
    print(f'Point #: {numPoint}')
    print(f'Element #: {eleNum}')
    print(f'relatedElement #: {comp}')
    if len(path) == 0:
        pass;
    else:
        newCompositions.to_excel(f'{path}/composition_for_feasibilityMap.xlsx')
    return Compositions, int(numPoint), comp, len(newCompositions)