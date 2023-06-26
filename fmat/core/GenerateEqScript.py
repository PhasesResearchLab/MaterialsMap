import pandas as pd
from tqdm import tqdm
import os
import numpy as np
##############################################
# Temperature in the setting file should be in Celcius
# default composition excel file should be in weight fraction
##############################################

def sortCompositions(Compositions): 
    """split the feasibility compositions into different groups with different lack elements
    Args:
        Compositions (list): List combinations of independent elements

    Returns:
        dict: Dicts of List combinations of independent elements with key as differeent element groups
    """
    output = []
    lackEle = []
    for index in tqdm(range(len(Compositions['Index']))):
        tempt = []
        keys = [item for item in Compositions.keys()]
        for key in keys:
            if key != 'Index' and 'Unnamed' not in key and 'Temperature' not in key and 'alloy' not in key:
                if round(Compositions[key][index],5) == 0:
                    tempt.append(key)
        if tuple(tempt) not in lackEle:
            lackEle.append(tuple(tempt))
    # lackEle = [item for item in set(lackEle)]
    print(lackEle)
    for index_lackEle in range(len(lackEle)):
        Elelist = {}
        for key in keys:
            if 'Unnamed' not in key:
                Elelist[key] = []
        output.append(Elelist)
    for index in tqdm(range(len(Compositions['Index']))):
        tempt = []
        for key in keys:
            if key != 'Index' and 'Unnamed' not in key and 'Temperature' not in key  and 'alloy' not in key:
                if round(Compositions[key][index],5) == 0:
                    tempt.append(key)
        tempt = tuple(tempt)   
        ii = lackEle.index(tempt)
        for key in keys:
            if 'Unnamed' not in key :
                output[ii][key].append(Compositions[key][index])
    return output

def findMainElement(compositions, elements, index):
    """when generating script, use the main element to be the dependent element

    Args:
        compositions (list): List combinations of independent elements
        elements (list): List of related elements
        index (int): the index of the main element in composition list

    Returns:
        list: add main element into element list
    """
    fraction = []
    for ele in elements:
        fraction.append(round(compositions[ele][index],5))
    return elements[fraction.index(max(fraction))]

def getSettings(path): 
    """get settings when generating the compositions of the feasibility map

    Args:
        path (str): path to the stored setting results

    Returns:
        set: rearrange the parameters in settings
    """
    ##############################get settings######################################
    settings = np.load(f'{path}/setting.npy',allow_pickle=True) #settings = [TemperatureRange,numPoint,numSimultion,relatedEles,terminalAlloys,indep_terminalAlloys,database,pressure]
    comp = settings[3]
    pressure = settings[7]
    database = settings[6]
    TRange = settings[0]
    numFile = int(settings[1])
    comp1 = settings[5][0]
    comp2 = settings[5][1]
    comps = settings[4]
    path = os.path.abspath(path)
    folder_Eq = f'{path}/Thermo-calc/Equilibrium Simulation'
    isExist = os.path.exists(folder_Eq + '/Result')
    if not isExist:
        os.makedirs(folder_Eq + '/Result')
        print("The new directory for Eq is created!")
    folder_Scheil = f'{path}/Thermo-calc/Scheil Simulation'
    isExist = os.path.exists(folder_Scheil + '/Result')
    if not isExist:
        os.makedirs(folder_Scheil + '/Result')
        print("The new directory for Scheil is created!")
    #get composition
    data = pd.read_excel(f'{path}/composition_for_feasibilityMap.xlsx')
    composition_data = dict()
    for item in data.columns:
        if 'alloy' in item:
            composition_data[item[6:]] = data[item].values.tolist()[:numFile]
    composition_data = pd.DataFrame(composition_data)
    Compositions = {}
    for item in data.keys():
        if 'Unnamed' not in item:
            Compositions[item] = data[item].values
    return [TRange, numFile, comp1, comp2, comp, folder_Eq, folder_Scheil, composition_data, Compositions, comps, pressure, database]
    

def createEqScript(path, maxNumSim = 999, database = None, eleAmountType = 'massFraction'): 
    """generate the TCM files for eq calculations based on settings

    Args:
        path (str): path to the stored setting results
        maxNumSim (int, optional): the maximum number of simulations in each TCM file. Defaults to 999.
        database (str, optional): the path to database or database name that indeside thermo_calc. Defaults to None.
        eleAmountType (str, optional): the element amount type. Defaults to 'massFraction'.

    Returns:
        TCM files: numScript(related to different comps)_comp(related elements in this script)_numFile(if exceed the maxNumSim, the script will be splited)
    """
    ################ get settings ################
    settings = getSettings(path) #[TRange, numFile, comp1, comp2, comp, folder_Eq, folder_Scheil, composition_data, Compositions, comps, pressure, database]
    path = os.path.abspath(path)
    Compositions = settings[8]
    comp = settings[4]
    if database == None:
        database = settings[11]
    output_Eq = settings[5]
    TemperatureRange = settings[0]
    pressure = settings[10]
    output = sortCompositions(Compositions)
    index_compositionList = 0
    for Compositions in output:
        compositions = {}
        relatedEles = []
        for ele in comp:
            compositions[ele] = Compositions[ele]
            if Compositions[ele][0] != 0:
                relatedEles.append(ele)
        fileName = ''
        for ele in relatedEles:
            fileName += f'{ele}-'
        Indexs = Compositions['Index']
        print(f'{len(Indexs)}th simulation in {index_compositionList}th file')
        f = open(output_Eq + '/' + f'{index_compositionList}_{fileName[:-1]}_0.TCM', 'w')
        typeNum = ['massFraction', 'moleFraction'].index(eleAmountType)
        if typeNum == 0:
            letter = 'W'
        else:
            letter = 'N'
        numSimulation = 0
        numFile = 0
        for index in tqdm(range(len(Indexs))):
            if len(Indexs) < maxNumSim:
                progress = round((numSimulation+1)/len(Indexs)*100,3)
            else:
                if len(Indexs) - maxNumSim*(numFile) < maxNumSim:
                    progress = round((numSimulation+1)/(len(Indexs) - maxNumSim*(numFile))*100,3)
                else:
                    progress = round((numSimulation+1)/maxNumSim*100,3)
            f.write(f'@@ {index}th simulation, {progress}% ##############################\n')
            # get all the elements will be used in this simulation
            allElements = ''
            for ele in comp:
                if round(compositions[ele][index],4) != 0:
                    allElements += ele + ' '
            if index == 0:
                f.write('set-echo\n\n')
            f.write('go data\n')
            if database[-4:].upper() != '.TDB':
                f.write(f'{database}\n')
            else:
                f.write(f'sw user "{database}"\n')
            f.write(f'def-sys {allElements}\n')
            f.write('get\n')
            f.write('go p-3\n')
            # reinitialize the conditions
            f.write('reinit\n')
            f.write(f's-c P={pressure} T={TemperatureRange[1]+273.15}\n')
            mainEle = findMainElement(compositions, comp, index)
            for ele in comp:
                if ele != mainEle and round(compositions[ele][index],5) != 0:
                    if typeNum in [1,3]:
                        f.write(f's-c {letter}({ele}) = {round(compositions[ele][index],5) / 100}\n')
                    else:
                        f.write(f's-c {letter}({ele}) = {round(compositions[ele][index],5)}\n')
            f.write('s-c N=1\n')
            f.write('c-e\n')
            # save the equilibrium result
            f.write('l-e\n')
            f.write('SCREEN\n')
            f.write('VWCS\n')
            f.write(f's-a-v 1 t {TemperatureRange[0]+273.15} {TemperatureRange[1]+273.15}\n')
            f.write(f'{TemperatureRange[2]}\n')
            f.write('step\n')
            f.write('NORMAL\n')
            f.write('enter function TC=T-273.15\n\n')
            f.write('ent tab\n')
            f.write('PhaseMole\n')
            f.write('TC NP(*) x(*,*)\n\n')
            f.write('ent tab\n')
            f.write('PhaseWt\n')
            f.write('TC BPW(*) w(*,*)\n\n')
            f.write('tabulate\n')
            f.write('PhaseMole\n')
            f.write('"'+path+'/Thermo-calc/Equilibrium Simulation/'+f'{Indexs[index]}_mole.exp"'+'\n')
            f.write('tabulate\n')
            f.write('PhaseWt\n')
            f.write('"'+path+'/Thermo-calc/Equilibrium Simulation/'+f'{Indexs[index]}_wt.exp"'+'\n')
            f.write('delete_sym\n')
            f.write('PhaseMole\n')
            f.write('delete_sym\n')
            f.write('PhaseWt\n')
            f.write('back\n')  
            f.write('back\n') 
            numSimulation += 1
            if numSimulation >= maxNumSim:
                numFile += 1
                f.close()
                numSimulation = 0
                f = open(output_Eq + '/' + f'{index_compositionList}_{fileName[:-1]}_{numFile}.TCM', 'w')
        f.close()
        index_compositionList += 1
    print('###################################################')
    print('All the equilibrium simulation scripts are created!')
    return None