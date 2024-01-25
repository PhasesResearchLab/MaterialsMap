from materialsmap.core.GenerateEqScript import getSettings
from tqdm import tqdm
import numpy as np
import json

def getColumn(words):
    columns = []
    for word in words:
        try:
            index1 = word.index('(') + 1
            index2 = word.index(')')
            columns.append(word[index1:index2])
        except:
            index1 = word.index('=') + 1
            columns.append(word[index1:-1])
    return columns

def sortOutputbyT(output,T):
    n = len(output[T])
    for i in range(n):
        # Last i elements are already in place
        for j in range(0, n-i-1):
            # Traverse the array from 0 to n-i-1
            # Swap if the element found is greater than the next element
            if output[T][j] > output[T][j+1] :
                for item in output.keys():
                    output[item][j], output[item][j+1] = output[item][j+1], output[item][j]
    return output

def readEqFromFile(folder_Eq, file_index, Mole = True): 
    """if Mole is True, then read mole fraction for phase, else read weight fraction for phase

    Args:
        folder_Eq (path): path to folder that contains eq results
        file_index (int): index of the TCM files
        Mole (bool, optional): Defaults to True.

    Returns:
        dict: eq results sorted by T
    """
    if Mole:
        fileName = f'{folder_Eq}/{file_index}_mole.exp'
    else:
        fileName = f'{folder_Eq}/{file_index}_wt.exp'
    f = open(fileName,'r+')
    lines = f.readlines()
    #############get all col###########
    columns = []
    for index in range(len(lines)):
        content = lines[index]
        words = content.split()
        if len(words) > 0 and 'col-1=' in words[0]:
            columns += getColumn(words)
    finalcolumns = set(columns)
    ############initialize output######
    output = dict()
    for item in finalcolumns:
        output[item] = []
    ###################################
    startRead = False
    for index in range(len(lines)):
        content = lines[index]
        words = content.split()
        if len(words) > 0 and 'col-1=' in words[0]:
            columns = getColumn(words)
            startRead = True
        if startRead and len(words) > 0 and 'col-1=' not in words[0]:
            if len(words) == len(columns) and 'NONE' not in words:
                for col in finalcolumns:
                    if col in columns:
                        colIndex = columns.index(col)
                        output[col].append(float(words[colIndex]))
                    else:
                        output[col].append(0)
            else:
                print(f'error in {file_index}th file, line {index}')
        if startRead and len(words) == 0:
            startRead = False
        output = sortOutputbyT(output, columns[0])
    return output

def transferTempToKelvin(data,numFile):
    for index in range(numFile):
        if data[f'Point{index}'] != None:
            data[f'Point{index}']['TK'] = [item + 273.15 for item in data[f'Point{index}']['TC']]
            del data[f'Point{index}']['TC']
    return data

def getEqdata(path, readMole = True):
    """read eq results from exp files
    if readMole is True, then read mole fraction for phase, else read weight fraction for phase

    Args:
        path (str): path to store the results
        readMole (bool, optional): _description_. Defaults to True.

    Returns:
        json: data_more.json that contains all eq results
    """
    print('###################################################################')
    print('####################### Reading Eq Result #########################')

    settings = getSettings(path) #[TRange, numFile, comp1, comp2, comps, folder_Eq, folder_Scheil, composition_data, Compositions, comp, pressure, database]
    folder_Eq = settings[5]
    numFile = settings[1]
    failList = []
    data_original = {}
    data = {}
    for index in tqdm(range(numFile)):
        # data_original[f'Point{index}'] = {}
        # data[f'Point{index}'] = {}
        try:
            Result_original = readEqFromFile(folder_Eq,index, readMole)
            Result = {}
            for key in Result_original.keys():
                if ',' not in key:
                    Result[key] = Result_original[key]
            if Result == {} or Result_original == {}:
                print(f'\n error in {index}th file, no result!')
                data_original[f'Point{index}'] = None
                data[f'Point{index}'] = None
                failList.append(index)
            else:
                data_original[f'Point{index}'] = {}
                data[f'Point{index}'] = {}
                for key in Result_original.keys():
                    data_original[f'Point{index}'][key] = Result_original[key]
                for key in Result.keys():
                    data[f'Point{index}'][key] = Result[key]
        except Exception as e:
            print(f'\n fail to read {index}th file, {e}!')
            data_original[f'Point{index}'] = None
            data[f'Point{index}'] = None
            failList.append(index)
    data = transferTempToKelvin(data,numFile)
    outputFile = json.dumps(data)
    if readMole:
        f = open(f'{folder_Eq}/Result/data_mole.json','w')
    else:
        f = open(f'{folder_Eq}/Result/data_wt.json','w')
    f.write(outputFile)
    f.close()
    print('####################### Reading Eq Result Done ####################')
    print(f'{len(failList)} files failed. FailList:')
    print(failList)
    return None
    # data = {'Point0':{'TC':[],'Phase1':[],'Phase2':[]},'Point1':{}......}

