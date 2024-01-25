import numpy as np
import json
import os
from materialsmap.core.GenerateEqScript import getSettings, sortCompositions, findMainElement
################################################################################
#==============================================================================#
############################Scheil script#######################################
# database = settings[-2]
# output_Scheil = f'{comp[0]}-{comps[1]}-{comps[2]}-Scheil'
# backupStartTemp = 1600
# temperatureStep = 1
# iterationNum = 1000
# finalLiquidFraction = 0.001
# maxNumSim = 250  # maximum number of simulations in each TCM file
# GlobalMinimization = True
# miscibilityGap = False # not ready for use yet
# miscibilityGapPhase = 'LIQUID'
# fastDiffusing = False
# fastDiffusingComp = 'C'
# retainPhase = False #False = consider all phase within the elements system; True = only consider the following phases
# if retainPhase:
#     Phases = 'LIQUID BCC_A2 FCC_A1 HCP_A3 SIGMA LAVES_PHASE_C14 NBNI3'
# rejectPhase = False
# if rejectPhase:
#     Phases = 'DELTA'
# rejectPhase and retainPhase, only one of them can be true
###############################################################
def getLiquidusTempforPoints(EqResult, backupStartTemp, numFile, liquidName):
    """Getting the liquidus temperature from eq results to start scheil simulations

    Args:
        EqResult (dict): Dict that contains the eq results
        backupStartTemp (int): the temperature that user define 
        numFile (int): the number of files
        liquidName (str): the name of liquid phase

    Returns:
        list: list of liquidus temperature based on eq results
    """
    LiquidusTemp = []
    for index in range(numFile):
        if EqResult == None or EqResult['Point'+str(index)] == None:
            LiquidusTemp.append(backupStartTemp)
        else:
            try:
                Liquid = EqResult['Point'+str(index)][liquidName]
                T = EqResult['Point'+str(index)]['TK']
                for i in range(len(Liquid)):
                    if Liquid[i] == 1:
                        LiquidusTemp.append(T[i]-273)
                        break
            except Exception as e:
                print(f'\n fail to read {index}th file, {e}!')
                LiquidusTemp.append(backupStartTemp)      
    return LiquidusTemp

def createScheilScript(path, backupStartTemp,liquidName = 'LIQUID', temperatureStep = 1,iterationNum = 2000,finalLiquidFraction = 0.001,GlobalMinimization = True, retainPhase = False,rejectPhase = False,fastDiffusing = False,maxNumSim = 999, database = None, eleAmountType = 'massFraction'):
    """generate the TCM files for scheil calculations based on settings

    Args:
        path (str): path to the stored setting results
        backupStartTemp (int): the temperature that user define 
        liquidName (str, optional): liquid phase name. Defaults to 'LIQUID'.
        temperatureStep (int, optional): Defaults to 1.
        iterationNum (int, optional): Defaults to 2000.
        finalLiquidFraction (float, optional): Defaults to 0.001.
        GlobalMinimization (bool, optional): Defaults to True.
        retainPhase (bool, optional): Defaults to False.
        rejectPhase (bool, optional): Defaults to False.
        fastDiffusing (bool, optional): Defaults to False.
        maxNumSim (int, optional): Defaults to 999.
        database (_type_, optional):  Defaults to the one in the settings.
        eleAmountType (str, optional): Defaults to 'massFraction'.

    Returns:
        TCM files: numScript(related to different comps)_comp(related elements in this script)_numFile(if exceed the maxNumSim, the script will be splited)
    """
    settings = getSettings(path) #[TRange, numFile, comp1, comp2, comp, folder_Eq, folder_Scheil, composition_data, Compositions, comps, pressure, database]
    if database == None:
        database = settings[11]
    path = os.path.abspath(path)
    comp = settings[4]
    Compositions = settings[8]
    ScheilFolder = settings[6]
    EqFolder = settings[5]
    ScheilFolder = settings[6]
    numFile = settings[1]
    try:
        try:
            f = open(f'{EqFolder}/Result/data_mole.json')
        except:
            f = open(f'{EqFolder}/Result/data_wt.json')
        EqResult = json.load(f)
    except:
        print('##################################################')
        print('No equilibrium result found! Will use the backUpStartTemp as liquidus temperature!')
        EqResult = None
    
    # print('##################Please input Scheil parameters#######################')
    # print('Liquid phase name (default: LIQUID):')
    # print('backUpStartTemp (in Celcius) (Scheil start temperature if liquidus temperature reading failed):')
    backupStartTemp = float(backupStartTemp)
    # print('temperatureStep (Scheil temperature step):')
    temperatureStep = float(temperatureStep)
    # print('iterationNum (maximum Scheil iteration limitation):')
    # while iterationNum.isdigit() == False:
    #     print('Please enter a integer !')
    #     print('iterationNum:')
    iterationNum = int(iterationNum)
    # print('finalLiquidFraction (Scheil final liquid fraction):')
    # while finalLiquidFraction.replace('.', '').isdigit() == False or float(finalLiquidFraction) > 1 or float(finalLiquidFraction) < 0:
    #     print('Please enter a number between 0 and 1!')
    #     print('finalLiquidFraction:')
    finalLiquidFraction = float(finalLiquidFraction)

    # print('retainPhase (False means considering all phases within the elements system) (y or n):')
    if retainPhase:
        print('Phases (seperate by space):')
        print('Example: LIQUID BCC_A2 FCC_A1 HCP_A3 SIGMA LAVES_PHASE_C14 NBNI3 (full phase name from TDB file)')
        print('Input:')
        Phases = input()
    if retainPhase == False:
        #print('rejectPhase (Any phase you want to exclude? False means considering all phases within the elements system) (y or n):')
        if rejectPhase:
            print('Phases (seperate by space):')
            print('Example: DELTA (full phase name from TDB file)')
            print('Input:')
            Phases = input()
    # print('fastDiffusing (Do you have fast diffusing elements?) (y or n):')

    if fastDiffusing:
        print('fastDiffusingComp (seperate by space):')
        print('Example: C ')
        print('Input:')
        fastDiffusingComp = input()
    print('##################Finish input Scheil parameters#######################')
    print('###############################################################################')
    print('Start to read liquidus temperature for each point from equilibrium results ...')
    LiquidusTemp = getLiquidusTempforPoints(EqResult, backupStartTemp, numFile, liquidName)
    print('###############################################################################')
    print('Finish reading liquidus temperature for each point from equilibrium results ...')
    print('##################Start to write Scheil script#######################')
    output = sortCompositions(Compositions)
    comp = [item.upper() for item in comp]
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
        f = open(ScheilFolder + '/' + f'{index_compositionList}_{fileName}_0.TCM', 'w')
        numSimulation = 0
        numFile = 0
        for index in range(len(Indexs)):
            startTemp = LiquidusTemp[Indexs[index]] + 25
            f.write('SET_ECHO\n\n\n')
            f.write('GO SCHEIL\n')
            if GlobalMinimization:
                f.write('global_minimization\n')
                f.write('Y\n')
            f.write('TERMINATION_CRITERIA\n')
            f.write('F\n')
            f.write(f'{finalLiquidFraction}\n')
            f.write('set-numerical-limits\n')
            f.write(f'{iterationNum}\n\n\n\n')
            f.write('TEMPERATURE-STEP\n')
            #set temperature step
            f.write(f'{temperatureStep}\n')
            f.write('START-WIZARD\n')
            if database[-3:].upper() != 'TDB':
                f.write(database)
            else:
                f.write('user ' + '"'+database+'"')
            mainEle = findMainElement(compositions, comp, index)
            f.write('\n')
            f.write(f'{mainEle}\n')
            #use mass percent
            f.write('y\n')
            includedEles = [mainEle]
            num_ele = 0
            for ele in comp:
                if ele != mainEle and round(compositions[ele][index] * 100, 5) != 0:
                    includedEles.append(ele)
                    f.write(f'{ele}\n')
                    f.write(f'{round(compositions[ele][index] * 100, 5)}\n') #should be in mass percent
            if num_ele == 0:
                for i in comp:
                    if i != mainEle:
                        f.write(f'{i}\n')
                        f.write('1E-5\n') #should be in mass percent         
                        break;       
            f.write('\n')
            f.write(f'{startTemp}\n')
            if rejectPhase:
                f.write(f'{Phases}\n')
                f.write(f'NONE\n')
            if not rejectPhase:
                f.write('*\n')
                if not retainPhase:
                    f.write('*\n')
                else:
                    f.write(f'{Phases}\n')
            f.write('NONE\n')
            f.write('y\n')
            f.write('N\n')
            if fastDiffusing and fastDiffusingComp in includedEles:
                f.write(f'{fastDiffusingComp}\n')
            else:
                f.write('NONE\n')
            f.write('s-d-a x t\n')
            f.write('s-d-a y NL\n')
            f.write('plot\n')
            f.write('m-e\n')
            f.write('"'+path+'/Thermo-calc/Scheil Simulation/'+f'{Indexs[index]}_liquid_mol%'+'"'+'\n')
            f.write('s-d-a y NS(*)\n')
            f.write('plot\n')
            f.write('m-e\n')
            f.write('"'+path+'/Thermo-calc/Scheil Simulation/'+f'{Indexs[index]}_solid_mol%'+'"'+'\n')
            f.write('s-d-a y BS(*)\n')
            f.write('plot\n')
            f.write('m-e\n')
            f.write('"'+path+'/Thermo-calc/Scheil Simulation/'+f'{Indexs[index]}_solid_wt%'+'"'+'\n')
            f.write('s-d-a y BL\n')
            f.write('plot\n')
            f.write('m-e\n')
            f.write('"'+path+'/Thermo-calc/Scheil Simulation/'+f'{Indexs[index]}_liquid_wt%'+'"'+'\n')
            f.write('s-d-a y W(*,*)\n')
            f.write('\n')
            f.write('plot\n')
            f.write('m-e\n')
            f.write('"'+path+'/Thermo-calc/Scheil Simulation/'+f'{Indexs[index]}_composition_wt%'+'"'+'\n')
            if index == numFile -1:
                f.write('set-inter\n')
            else:
                f.write('back\n')
                f.write('back\n')
            numSimulation += 1
            if numSimulation >= maxNumSim:
                numFile += 1
                f.close()
                numSimulation = 0
                f = open(ScheilFolder + '/' + f'{index_compositionList}_{fileName}_{numFile}.TCM', 'w')
        f.close()
        index_compositionList += 1
    return None
