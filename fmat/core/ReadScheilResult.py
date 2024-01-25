import json
from materialsmap.core.GenerateEqScript import getSettings
from sklearn import neighbors
import numpy as np
from tqdm import tqdm
import pandas as pd

def getPhaseNamesInSequence(fileName,liquidPhase):
    """rearrange the scheil results based on phase names

    Args:
        fileName (str): name of files that contains the scheil results
        liquidPhase (str): name of liquid phase

    Returns:
        list: list of scheil phase name
    """
    PhaseNames = []
    f = open(fileName,'r')
    lines = f.readlines()
    for i in range(len(lines)):
        words = lines[i].split()
        if words != [] and words[0] == '$E' and len(words) > 1:
            phase = words[1:]
            if liquidPhase not in phase[0] and phase[0] not in PhaseNames:
                PhaseNames.append(phase[0])
    f.close()
    return PhaseNames

def linkPhaseAndTemp(currentPhases, currentTemp):
    """create dicts that have both phase name and T

    Args:
        currentPhases (dict)): dict of all realted phases and fraction
        currentTemp (dict): dict of temperature

    Returns:
        dict: dicts of all realted phases, fraction and temperature
    """
    Models = []
    n_neighbors = 2
    weights = 'distance'
    phases = [item for item in currentPhases.keys()]
    for phase in phases:
        input = np.array(currentTemp[phase]).reshape(-1, 1)
        output = currentPhases[phase]
        Models.append(neighbors.KNeighborsRegressor(n_neighbors, weights=weights).fit(input, output))
    result = dict()
    result['TC'] = currentTemp[phase]
    for phase in phases:
        result[phase] = Models[phases.index(phase)].predict(np.array(result['TC']).reshape(-1, 1)).tolist()
    return result

def getAllPhases(ScheilResult):
    """read scheil results from exp files

    Args:
        ScheilResult (dict): dicts that contains the scheil results

    Returns:
        set: a set of unique phases
    """
    phases = []
    # i = 0
    for item in ScheilResult.keys():
        item = ScheilResult[item]
        if item != None:
            temptList = [thing for thing in item.keys() if thing != 'TC']
            phases += temptList
    phases = set(phases)
    return phases

def getFinalScheilResult(ScheilResult,folder_Scheil,numFile):
    """get scheil results from exp

    Args:
        ScheilResult (dict): dict of scheil results
        folder_Scheil (str): path to open/store the results
        numFile (int): number of files in the scheil folder

    Returns:
        xlsx: an excel sheet that contains the scheil results
    """
    print('####################################################################')
    print('Getting final Scheil result...')
    phaseNames = getAllPhases(ScheilResult)
    finalPhases = dict()
    finalPhases['Point'] = []
    for phase in phaseNames:
        finalPhases[phase] = []
    list = []
    for index in tqdm(range(numFile)):
        phaseAmount = ScheilResult[f'Point{index}']
        if phaseAmount != None:
            finalPhases['Point'].append(index)
            for phase in phaseNames:
                if phase in phaseAmount.keys():
                    try:
                        if phaseAmount[phase] != []:
                            finalPhases[phase].append(phaseAmount[phase][-1])
                        else:
                            finalPhases[phase].append(0)
                    except Exception as e:
                        print(f'error in reading {phase} in Point {index}', e)
                        list.append(index)
                else:
                    finalPhases[phase].append(0)
            sum = 0 
            for phase in phaseNames:
                sum += finalPhases[phase][-1]
            # sometimes the Scheil result will end up with a sum of 1.0000000000000002
            # in this case, we will use the second last value
            if sum > 1:
                for phase in phaseNames:
                    if phase in phaseAmount.keys():
                        try:
                            if phaseAmount[phase] != []:
                                finalPhases[phase][-1] = phaseAmount[phase][-2]
                            else:
                                finalPhases[phase][-1] = 0
                        except Exception as e:
                            print(index, e)
                            list.append(index)
                    else:
                        finalPhases[phase][-1] = 0
        else:
            finalPhases['Point'].append(index)
            for phase in phaseNames:
                finalPhases[phase].append('')
    # print('FailList:',list)
    # print(f'{len(list)} points failed to read')
    result = pd.DataFrame(finalPhases)
    result.to_excel(f'{folder_Scheil}/Result/ScheilResults.xlsx')
    print('##########################end reading Scheil#############################')
    return list

def readLiqAndSolT(folder_Scheil, numFile):
    """read the liquidius and solidius temperature from exp files

    Args:
        folder_Scheil (str): path to open/store the results
        numFile (int): number of files in the scheil folder

    Returns:
        dict: dict of liquidius and solidius temperature
        list: list of fail to read the results
    """
    print('Read Scheil Temperature vs Liquid fraction (mole):')
    result = dict()
    temperature = []
    liquidFraction = []
    x = []
    index_fail = []
    for index in range(numFile):
        temp = []
        liquid = []
        fileName = f'{folder_Scheil}/{index}_liquid_mol%.exp'
        try:
            data = open(fileName,'r')
        except:
            result[f'Point{index}'] = None
            index_fail.append(index)
            print(f'cannot find {index}th file')
            continue
        lines = data.readlines()
        startRead = False
        for i in range(len(lines)):
            content = lines[i]
            words = content.split()
            if i < len(lines) - 1:
                nextContent = lines[i + 1]
                # print(words[-1],',', nextContent,'\n')
            if words!= [] and words[-1] == 'M':
                startRead = True
            if words == ['BLOCKEND'] or words == ['CLIP', 'OFF'] or words[:2] == ['$', 'Y-AXIS:'] and startRead:
                startRead = False
            if startRead:
                temp.append(float(words[0]))
                liquid.append(float(words[1]))  
        if liquid != [] and temp != []: 
            result[f'Point{index}'] = dict() 
            result[f'Point{index}']['TC'] = temp
            result[f'Point{index}']['LIQUID'] = liquid
        else:
            print(f'No result in {index}th simulation')
            index_fail.append(index)
    return result, index_fail

def combineLiqAndSolT(ScheilLiquidResult, ScheilResult, numFile):
    """combined scheil liquid phase fraction and liquidius and solidius temperature

    Args:
        ScheilLiquidResult (dict): dcit of scheil liquidius and solidius temperature
        ScheilResult (dict): dict of scheil liquid phase fraction
        numFile (int): number of files in the scheil folder

    Returns:
        dict: dict of combined results
    """
    finalResult = dict()
    for index in range(numFile):
        solid = ScheilResult[f'Point{index}']
        liquid = ScheilLiquidResult[f'Point{index}']
        if solid == None and liquid == None:
            finalResult[f'Point{index}'] = None
        if solid == None and liquid != None:
            finalResult[f'Point{index}'] = liquid
        if solid != None and liquid == None:
            finalResult[f'Point{index}'] = solid
        if solid != None and liquid != None:
            finalResult[f'Point{index}'] = solid
            input = np.array(liquid['TC']).reshape(-1, 1)
            output = liquid['LIQUID']
            model = neighbors.KNeighborsRegressor(2, weights='distance').fit(input, output)
            finalResult[f'Point{index}']['LIQUID'] = model.predict(np.array(solid['TC']).reshape(-1, 1)).tolist()
    return finalResult

def transferTempToKelvin(data,numFile):
    for index in range(numFile):
        if data[f'Point{index}'] != None:
            data[f'Point{index}']['TK'] = [item + 273.15 for item in data[f'Point{index}']['TC']]
            del data[f'Point{index}']['TC']
    return data

def getScheilSolidPhase(path, liquidPhase = 'LIQUID'):
    """read scheil results from exp file

    Args:
        path (str): path to store the results
        liquidPhase (str, optional): name of liquid phase. Defaults to 'LIQUID'.

    Returns:
        json: data_more.json that contains all scheil results
    """
    print('#####################################################')
    print('#############start reading Scheil Result#############')
    #####################read settings######################
    settings = getSettings(path) #[TRange, numFile, comp1, comp2, comps, folder_Eq, folder_Scheil, composition_data, Compositions, comp, pressure, database]
    folder_Scheil = settings[6]
    numFile = settings[1]
    #####################read data###########################
    ScheilResult = dict()
    index_noFile = []
    index_failRead = []
    for index in tqdm(range(numFile)):
        ########################################################################
        ###########################read solid phase#############################
        fileName = f'{folder_Scheil}/{index}_solid_mol%.exp'
        try:
            phaseNames = getPhaseNamesInSequence(fileName, liquidPhase)
            data = open(fileName,'r')
        except:
            print(f'\n {fileName} does not exist, skip')
            index_noFile.append(index)
            ScheilResult[f'Point{index}'] = None
            continue
        try:
            currentPhases = dict()
            currentTemp = dict()
            for phase in phaseNames:
                currentPhases[phase] = []
                currentTemp[phase] = []
            startBlock = False
            startRead = False
            BlockNum = 0
            lines = data.readlines()
            for i in range(len(lines)):
                words = lines[i].split()
                if words != [] and words[:2] == ['$','BLOCK']:
                    BlockNum += 1
                    startBlock = True
                    phaseCount = 0
                if startBlock and words != [] and words[-1] == 'M':
                    startRead = True
                if startRead and words != [] and words == ['CLIP', 'OFF']:
                    startRead = False
                    phaseCount += 1
                if startRead and words != [] and words[:2] == ['$', 'Y-AXIS:']:
                    startRead = False
                    phaseCount += 1            
                if startRead and startBlock and words != [] and words[0] == 'BLOCKEND':
                    startRead = False
                    startBlock = False
                if startRead:
                    phase = phaseNames[phaseCount]
                    currentTemp[phase].append(float(words[0]))
                    currentPhases[phase].append(float(words[1]))
            #######################uniform temperatures for each phase##########
            ScheilResult[f'Point{index}'] = linkPhaseAndTemp(currentPhases, currentTemp)
            ####################################################################
        except Exception as e:
            print(f'fail to read {fileName}, skip',e)
            index_failRead.append(index)
            ScheilResult[f'Point{index}'] = None
            continue
    # print('Failed Index List: ',index_failRead, len(index_failRead))
    ############################################################################
    ########################end of reading solid phase##########################
    ############################################################################
    #################read the final Scheil phase result#########################
    failList2 = getFinalScheilResult(ScheilResult,folder_Scheil, numFile)
    ScheilLiquidResult, failList3 = readLiqAndSolT(folder_Scheil, numFile)
    finalScheilResult = combineLiqAndSolT(ScheilLiquidResult, ScheilResult, numFile)
    finalScheilResult = transferTempToKelvin(finalScheilResult,numFile)
    outputFile = json.dumps(finalScheilResult)
    f = open(f'{folder_Scheil}/Result/data_mole.json','w')
    f.write(outputFile)
    f.close()
    totalfailList = list(set(index_noFile + index_failRead + failList2 + failList3))
    print(f'Failed Index List: {totalfailList}')
    print(f'{len(totalfailList)} Files Failed to Read')
    return None
