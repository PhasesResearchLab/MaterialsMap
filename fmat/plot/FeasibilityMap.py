import json
import os
from sklearn import neighbors
from materialsmap.core.GenerateEqScript import getSettings
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import materialsmap.plot.feasibility_helpers
import math


def readResult(folder_Scheil,folder_Eq, readMole = True):
    """read the previously saved Scheil and Eq result

    Args:
        folder_Scheil (str): path to folder that contains scheil results
        folder_Eq (str): path to folder that contains eq results
        readMole (bool, optional): _description_. Defaults to True.

    Returns:
        dcit: dict of eq and scheil results
    """
    if readMole:
        f_Scheil = open(f'{folder_Scheil}/Result/data_mole.json')
        f_Eq = open(f'{folder_Eq}/Result/data_mole.json')
    else:
        f_Scheil = open(f'{folder_Scheil}/Result/data_wt.json')
        f_Eq = open(f'{folder_Eq}/Result/data_wt.json')
    ScheilResult = json.load(f_Scheil)
    EqResult = json.load(f_Eq)
    return ScheilResult, EqResult

def readDynamicFeasibility(ScheilResult, EqResult, ratio = 2/3):
    """get dynamic Eq result from ratio*ScheilSolidusT to ScheilSolidusT

    Args:
        ScheilResult (dcit): dict of scheil results
        EqResult (dcit): dict of eq results
        ratio (float): ratio of T/Tmelt. Defaults to 2/3.

    Returns:
        dict: eq results only consider custmoize T and higher
    """

    numFile = len([item for item in ScheilResult.keys()])
    newEqResult = {}
    for i in range(numFile):
        Scheil = ScheilResult[f'Point{i}']
        Eq = EqResult[f'Point{i}']
        if Scheil == None:
            newEqResult[f'Point{i}'] = 'No Scheil Result'
        elif Eq == None:
            newEqResult[f'Point{i}'] = 'No Eq Result'
        else:
            T = Scheil['TK']
            solidusT = min(T)
            T_Eq = Eq['TK']
            if solidusT * ratio < min(T_Eq):
                newEqResult[f'Point{i}'] = 'No Eq Result at low temperature'
            else:
                newEqResult[f'Point{i}'] = dict()
                for key in Eq.keys():
                    newEqResult[f'Point{i}'][key] = []
                for key in Eq.keys():
                    for index in range(len(Eq[key])):
                        if Eq['TK'][index] >= solidusT * ratio and Eq['TK'][index] <= solidusT:
                            newEqResult[f'Point{i}'][key].append(Eq[key][index])
    return newEqResult

def getPhases(result): 
    """get related phases in the simulations for the currect materials system

    Args:
        result (dict): dict of eq/scheil results

    Returns:
        list: list of phase names
    """
    numFile = len([item for item in result.keys()])
    phases = []
    for i in range(numFile):
        if isinstance(result[f'Point{i}'],dict):
            for key in result[f'Point{i}'].keys():
                if key != 'TK' and key not in phases:
                    phases.append(key)
    return phases

def findMaxUnallowedPhaseEq(EqResult, allowedPhases): 
    """find the max amount of unallowed phase for eq results

    Args:
        EqResult (dict): dict of eq results
        allowedPhases (list): list of allowed phases

    Returns:
        list: a list of max unallowed phase amount for each point (Equilibrium)
    """
    #return a list of max unallowed phase amount for each point (Equilibrium)
    numFile = len([item for item in EqResult.keys()])
    EqMaxBadPhaseAmount = []
    for i in range(numFile):
        Eq = EqResult[f'Point{i}']
        if not isinstance(Eq, dict):
            EqMaxBadPhaseAmount.append(Eq)
        else:
            allowedPhaseAmount = []
            for key_index in range(len(Eq['TK'])):
                total = 0
                total_a = 0
                for key in Eq.keys():
                    if key in allowedPhases:
                        total_a += float(Eq[key][key_index])
                    elif key != 'TK':
                        total += float(Eq[key][key_index])
                allowedPhaseAmount.append(total)
            EqMaxBadPhaseAmount.append(max(allowedPhaseAmount))
    return EqMaxBadPhaseAmount

def getFinalScheilResult(ScheilResult): 
    """get the final Scheil result for each point

    Args:
        ScheilResult (dict): dict of scheil results

    Returns:
        dcit: get the final Scheil result for each point
    """

    numFile = len([item for item in ScheilResult.keys()])
    phaseNames = getPhases(ScheilResult)
    finalPhases = dict()
    finalPhases['Point'] = []
    for phase in phaseNames:
        finalPhases[phase] = []
    list = []
    for index in range(numFile):
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
            # somemonoserif the Scheil result will end up with a sum of 1.0000000000000002
            # in this case, we will use the second last value
            if sum > 1:
                for phase in phaseNames:
                    if phase in phaseAmount.keys():
                        try:
                            if phaseAmount[phase] != []:
                                finalPhases[phase][-1] = phaseAmount[phase][-1]/sum
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
    return finalPhases

def findMaxUnallowedPhaseScheil(finalScheilResults, allowedPhases): 
    """find the max amount of unallowed phase for scheil results

    Args:
        finalScheilResults (dict): dict of scheil results
        allowedPhases (list): list of allowed phases

    Returns:
        list: a list of max unallowed phase amount for each point (Scheil)
    """
    #return a list of max unallowed phase amount for each point (Scheil)
    ScheilMaxBadPhaseAmount = []
    for i in range(len(finalScheilResults['Point'])):
        hasScheilResult = True
        for key in finalScheilResults.keys():
            if finalScheilResults[key][i] == '':
                hasScheilResult = False
                break
        if not hasScheilResult:
            ScheilMaxBadPhaseAmount.append('No Scheil Result')
            continue
        total = 0
        total_a = 0
        for key in finalScheilResults.keys():
            if key in allowedPhases:
                total_a += finalScheilResults[key][i]
            elif key != 'TK' and key != 'Point':
                total += float(finalScheilResults[key][i])
        ScheilMaxBadPhaseAmount.append(total)
    return ScheilMaxBadPhaseAmount

def getSolidLiquidTFromScheil(ScheilResult,solidCriterion):
    """read solidius and liquidius temperature from scheil results

    Args:
        ScheilResult (dict: dict of scheil results
        solidCriterion (float): fraction for solid phase fraction

    Returns:
        list: list of solidus and liquidius temperature
    """
    solidusT = []
    liquidusT = []
    numFile = len([item for item in ScheilResult.keys()])
    for index in range(numFile):
        Scheil = ScheilResult[f'Point{index}']
        if Scheil == None:
            solidusT.append(None)
            liquidusT.append(None)
        else:
            LIQUID = Scheil['LIQUID']
            if max(LIQUID) < 0.999:
                liquidusT.append(None)
            if min(LIQUID) > solidCriterion:
                solidusT.append(None)
            if max(LIQUID) >= 0.999:
                if LIQUID[0] > 0.5:
                    i = 0
                    for item in LIQUID:
                        if item >= 0.999:
                            i += 1
                        else:
                            break
                    if i == len(LIQUID):
                        liquidusT.append(Scheil['TK'][i-1])
                    else:
                        liquidusT.append(Scheil['TK'][i])
                else:
                    i = 0
                    for item in LIQUID:
                        if item >= 0.999:
                            i += 1
                        else:
                            break
                    if i == len(LIQUID):
                        liquidusT.append(Scheil['TK'][i-1])
                    else:
                        liquidusT.append(Scheil['TK'][i])
            if min(LIQUID) <= solidCriterion:
                if LIQUID[0] > 0.5:
                    i = 0
                    for item in LIQUID:
                        if item > 0.001:
                            i += 1
                        else:
                            break
                    if i == len(LIQUID):
                        solidusT.append(Scheil['TK'][i-1])
                    else:
                        solidusT.append(Scheil['TK'][i])
                else:
                    i = 0
                    for item in LIQUID:
                        if item <= 0.001:
                            i += 1
                        else:
                            break
                    if i == len(LIQUID):
                        solidusT.append(Scheil['TK'][i-1])
                    else:
                        solidusT.append(Scheil['TK'][i])
    return solidusT, liquidusT

def plotMaps(path,engine,dynamicTRange = True, dynamicRatio = 2/3, ScheilThreshold = 0.05, EqThrshold = 0.1, allowPhase = ['FCC','BCC','HCP','LIQUID'],solidCriterion = 0.001, hotTeartSettings = {'numDataThreshold':10,'CSCPoints':[0.4,0.9,0.99], 'KouPoints':[0.93,0.98], 'CDPoints':[0.7,0.98]}): 
    """plot all realted figures

    Args:
        path (str): path to open the setting and store the results
        engine (str): computational engine, 'pycalphad' or 'thermo_calc'
        dynamicTRange (bool, optional): Defaults to True.
        dynamicRatio (float, optional): Defaults to 2/3.
        ScheilThreshold (float, optional): Defaults to 0.05.
        EqThrshold (float, optional): Defaults to 0.1.
        allowPhase (list, optional): Defaults to ['FCC','BCC','HCP','LIQUID'].
        solidCriterion (float, optional): Defaults to 0.001.
        hotTeartSettings (dict, optional): Defaults to {'numDataThreshold':10,'CSCPoints':[0.4,0.9,0.99], 'KouPoints':[0.93,0.98], 'CDPoints':[0.7,0.98]}.

    Returns:
        TIF: all figures store in path
    """
    #input path(path to simulation result), dynamicTRange (should use Eq T range according to Scheil result?), dynamicRatio (ScheilSolidT*dynamicRatio to ScheilSolidT), ScheilThreshold = 0.05, EqThrshold = 0.1, allowPhase = ['FCC','BCC','HCP','LIQUID']
    ###############################get settings######################################
    path = os.path.abspath(path)
    settings = getSettings(path)
    xComp = settings[2]
    yComp = settings[3]
    comps = settings[9]
    composition_data = settings[7]
    ###############################which source to use################################
    if engine.lower() == 'pycalphad':
        path = path + '/Pycalphad'
    elif engine.lower() == 'thermo-calc' or engine.lower() == 'thermo_calc' or engine.lower() == 'thermocalc':
        path = path + '/Thermo-calc'
    else:
        raise Exception('Please choose the calculations engine')
    folder_Scheil = path + '/Scheil Simulation'
    folder_Eq = path + '/Equilibrium Simulation'
    #######################read data from the simulation result#######################
    ScheilResult, EqResult = readResult(folder_Scheil,folder_Eq)
    if dynamicTRange:
        EqResult = readDynamicFeasibility(ScheilResult, EqResult, dynamicRatio)
    finalScheilResults = getFinalScheilResult(ScheilResult)
    EqPhases = getPhases(EqResult)
    ScheilPhases = getPhases(ScheilResult)
    allPhases = list(set(EqPhases + ScheilPhases))
    allowedPhases = []
    for phase in allPhases:
        for item in allowPhase:
            if item in phase:
                allowedPhases.append(phase)
    ScheilMaxBadPhaseAmount = findMaxUnallowedPhaseScheil(finalScheilResults, allowedPhases)
    EqMaxBadPhaseAmount = findMaxUnallowedPhaseEq(EqResult, allowedPhases)
    ###############################plotting##########################################
    coord = []
    for index in range(len(composition_data)):
        x = composition_data[xComp].values[index]
        y = composition_data[yComp].values[index]
        coord.append((x, y))
    
    plotScheilEqFeasibilityMap(path,coord,EqMaxBadPhaseAmount,ScheilMaxBadPhaseAmount,EqThrshold,ScheilThreshold,xComp,yComp,dynamicTRange,comps,dynamicRatio)
    solidusT,liquidusT = getSolidLiquidTFromScheil(ScheilResult,solidCriterion)
    plotScheilTemperature(path,solidusT,liquidusT,coord,xComp,yComp,solidCriterion)
    plotScheilPhase(path, finalScheilResults,allowedPhases, coord,xComp,yComp)
    ##############################hot tearing criteria###############################
    FR, CSC, Kou, iCSC, sRDG = getCriteria(ScheilResult,solidusT,liquidusT,numDataThreshold = hotTeartSettings['numDataThreshold']
                                           ,CSCPoints = hotTeartSettings['CSCPoints'],KouPoints = hotTeartSettings['KouPoints'],CDPoints = hotTeartSettings['CDPoints'])
    plotHotTearingSusceptibilityMap(path, coord, xComp, yComp, FR, CSC, Kou, iCSC, sRDG)
    return None

def plotScheilEqFeasibilityMap(path,coord,EqMaxBadPhaseAmount,ScheilMaxBadPhaseAmount,EqThrshold,ScheilThreshold,xComp,yComp,dynamicTRange,comps,dynamicRatio):
    print('####################################################################')
    print('Plotting Scheil-Eq Feasibility Map')
    if not dynamicTRange:
        handles = [
            mpl.patches.Patch(facecolor='purple'),
            mpl.patches.Patch(facecolor='red'),
            mpl.patches.Patch(facecolor='blue'),
            mpl.patches.Patch(facecolor='green'),
            mpl.patches.Patch(facecolor='black')
        ]
        labels = [
            'Equilibrium and Scheil infeasible',
            'Equilibrium infeasible, Scheil feasible',
            'Equilibrium feasible, Scheil infeasible',
            'Both feasible',
            'No Scheil/Eq data'
        ]
    else:
        handles = [
            mpl.patches.Patch(facecolor='purple'),
            mpl.patches.Patch(facecolor='red'),
            mpl.patches.Patch(facecolor='blue'),
            mpl.patches.Patch(facecolor='green'),
            mpl.patches.Patch(facecolor='black'),
            mpl.patches.Patch(facecolor='yellow')
        ]
        labels = [
            'Equilibrium and Scheil infeasible',
            'Equilibrium infeasible, Scheil feasible',
            'Equilibrium feasible, Scheil infeasible',
            'Both feasible',
            'No Scheil/Eq data',
            'No Eq Result at low temperature'
        ]
    ###############################plot the figure##################################
    dotSize = 22
    fig = plt.figure(figsize = (4,4), dpi = 300)
    ax = fig.add_subplot(projection='triangular')
    for index in range(len(coord)):
        (x_plot, y_plot) = coord[index]
        ScheilPointResult = ScheilMaxBadPhaseAmount[index]
        EqPointResult = EqMaxBadPhaseAmount[index]
        if (isinstance(ScheilPointResult,float) or isinstance(ScheilPointResult,int)) and (isinstance(EqPointResult,float) or isinstance(EqPointResult,int)):
            eq_is_feasible = EqPointResult < EqThrshold
            scheil_is_feasible = ScheilPointResult < ScheilThreshold
            if not eq_is_feasible and not scheil_is_feasible:
                ax.scatter(x_plot, y_plot,s=dotSize, label='Equilibrium and Scheil infeasible', c='purple', marker='h')
                continue
            if not eq_is_feasible and scheil_is_feasible:
                ax.scatter(x_plot, y_plot,s=dotSize, label='Equilibrium infeasible, Scheil feasible', c='red', marker='h')
                continue
            if not scheil_is_feasible and eq_is_feasible:
                ax.scatter(x_plot, y_plot,s=dotSize, label='Equilibrium feasible, Scheil infeasible', c='blue', marker='h')
                continue
            ax.scatter(x_plot, y_plot, label='Both Feasible',s=dotSize, c='green', marker='h')
        else:
            if ScheilPointResult == 'No Scheil Result' or EqPointResult == 'No Eq Result' or EqPointResult == 'No Scheil Result' or EqPointResult == None:
                ax.scatter(x_plot, y_plot,s=dotSize, label='No Scheil/Eq data', c='black')
                continue
            if dynamicTRange and EqPointResult == 'No Eq Result at low temperature':
                ax.scatter(x_plot, y_plot,s=dotSize, label='No Eq Result at low temperature', c='yellow')
                continue

    fmtted_comps = '-'.join(sorted(set(comps)))
    if dynamicTRange:
        ax.set_title(f"{fmtted_comps}_DynamicRatio{round(dynamicRatio,3)}\nEq_Tolerance: {EqThrshold}\nScheil_Tolerance: {ScheilThreshold}", fontname = 'monoserif', fontsize = 10)
    else:
        ax.set_title(f"{fmtted_comps}_nonDynamic\nEq_Tolerance: {EqThrshold}\nScheil_Tolerance: {ScheilThreshold}", fontname = 'monoserif', fontsize = 10)
    ax.set_xlabel(f'W({xComp})', fontname = 'monoserif', fontsize = 10)
    ax.set_ylabel(f'W({yComp})', labelpad=-50, fontname = 'monoserif', fontsize = 10)
    fig.legend(handles=handles, labels=labels, loc='lower left',prop={'family':'monoserif', 'size':10}, bbox_to_anchor=(0.8, 0.6))
    ax.tick_params(labelsize=10)
    if dynamicTRange:
        fig.savefig(f'{path}/{fmtted_comps}-EqScheil-dynamic_{round(dynamicRatio,3)}.tif', bbox_inches='tight')
    else:
        fig.savefig(f'{path}/{fmtted_comps}-EqScheil-nonDynamic.tif', bbox_inches='tight')
    plt.close()
    print(f'Plotting {fmtted_comps} done!')
    return None

def plotScheilTemperature(path,solidusT,liquidusT,coord,xComp,yComp,solidCriterion):
    """plot solidius and liquidus T from scheil results

    Args:
        path (str): path to open the setting and store the results
        solidusT (list): list of solidius T
        liquidusT (list): list of liquidius T
        coord (list): list of coordation 
        xComp (str): name for x label
        yComp (str): name for y label
        solidCriterion (float, optional): Defaults to 0.001.

    Returns:
        TIF: store the reusults as TIF files
    """
    print('####################################################################')
    print('Plotting Scheil-Eq Temperature Map')
    dotSize = 20
    TName = ['solidusT (K)','liquidusT (K)']
    T_index = 0
    for T in [solidusT, liquidusT]:
        plt.figure(figsize=(4, 4), dpi=400)
        subFig = plt.subplot(projection='triangular')
        T_index += 1
        values = [item for item in T if item != None]
        top = max(values)
        bottom = min(values)
        norm = mpl.colors.Normalize(vmin = bottom, vmax = top)
        cmap = 'Greys'
        for i in range(len(coord)):
            value = T[i]
            (x,y) = coord[i]
            if value != None:
                RGB1 = mpl.cm.Greys(norm(value), bytes = True)
                RGB = (RGB1[0]/255,RGB1[1]/255,RGB1[2]/255)
            else:
                RGB = 'yellow'
            subFig.scatter(x, y,s=dotSize, color = RGB, marker = "o")
        cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        orientation='horizontal',fraction=0.035, pad=0.2,aspect=20)
        for t in cbar.ax.get_xticklabels():
            t.set_fontsize(10)
            t.set_fontname('monoserif')
        if T_index == 1:
            plt.title(f'{TName[T_index-1]}\nsolidCriterion:{solidCriterion}\n{round(bottom,3)}-{round(top,3)}',fontsize = 10, fontname = 'monoserif')
        else:
            plt.title(f'{TName[T_index-1]}\n{round(bottom,3)}-{round(top,3)}',fontsize = 10, fontname = 'monoserif')
        plt.xlabel(f'X({xComp})',fontsize = 10, fontname = 'monoserif')
        plt.ylabel(f'X({yComp})', labelpad=-40,fontsize = 10, fontname = 'monoserif')

        plt.tick_params(labelsize=10)
        plt.savefig(f'{path}/{TName[T_index-1]}.tif', bbox_inches='tight')
    plt.close()
    print('Plotting Scheil-Eq Temperature Map done!')
    return None

def plotScheilPhase(path, finalScheilResults,allowedPhases, coord,xComp,yComp):
    """plot different phases from scheil simulations

    Args:
        path (str): path to open the setting and store the results
        finalScheilResults (dict): dict of scheil results
        allowedPhases (list): list of allowed phases
        coord (list): list of coordation 
        xComp (str): name for x label
        yComp (str): name for y label

    Returns:
        TIF: store the reusults as TIF files
    """
    print('####################################################################')
    print('Plotting Scheil Phase Heat Map')
    phase_data = finalScheilResults
    phaseNames = [item for item in phase_data.keys() if item != 'Point']
    phaseNum = len(phaseNames)
    fig = plt.figure(figsize = (3 * phaseNum,3),dpi=200)
    grid = plt.GridSpec(1, phaseNum,figure = fig)
    phase_index = 0
    for phase in phaseNames:
        subFig = plt.subplot(grid[0,phase_index], projection='triangular')
        phase_index += 1
        values = [item for item in phase_data[phase] if item != '']
        norm = mpl.colors.Normalize(vmin = 0, vmax = max(values))
        cmap = 'Greys'
        isBad = True
        if phase in allowedPhases:
                isBad = False
        if isBad:
            cmap = 'Reds'
        for i in range(len(coord)):
            value = phase_data[phase][i]
            (x,y) = coord[i]
            if value != '':
                if isBad:
                    RGB1 = mpl.cm.Reds(norm(value), bytes = True)
                else:
                    RGB1 = mpl.cm.Greys(norm(value), bytes = True) 
                RGB = (RGB1[0]/255,RGB1[1]/255,RGB1[2]/255)
            else:
                RGB = 'yellow'
            subFig.scatter(x, y, color = RGB, marker = "o",s = 5)

        cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            orientation='horizontal',fraction=0.035, pad=0.2,aspect=20)
        for t in cbar.ax.get_xticklabels():
            t.set_fontsize(10)
            t.set_fontname('monoserif')
        if phase == 'LIQUID':
            plt.title('Remaining Liquid After Scheil', fontsize = 10, fontname = 'monoserif')
        else:
            plt.title(phase, fontsize = 10, fontname = 'monoserif')
        plt.xlabel(f'W({xComp})',fontsize = 10, fontname = 'monoserif')
        plt.ylabel(f'W({yComp})',labelpad=-12.5,fontsize = 10, fontname = 'monoserif')

        plt.xticks([0,0.2,0.4,0.6,0.8,1],[0,0.2,0.4,0.6,0.8,1],fontsize=10, fontname = 'monoserif')
        plt.yticks([0,0.2,0.4,0.6,0.8,1],[0,0.2,0.4,0.6,0.8,1],fontsize=10, fontname = 'monoserif')

    plt.tight_layout()

    plt.savefig(f'{path}/ScheilPhaseHeatMap.tif',bbox_inches = 'tight')
    plt.close()
    print('Plotting Scheil Phase Heat Map done!')
    return None

#######################################cracking criteria########################
def getIntegral(temperature, solidFraction):
    """combine temperature and solid phase

    Args:
        temperature (list): list of temperature
        solidCriterion (float, optional): Defaults to 0.001.

    Returns:
        dict: combined results of T and phase fraction
    """
    result = 0
    if len(temperature) > 3:
        for i in range(len(temperature)):
            if i == 0:
                delT = abs((temperature[i+1] - temperature[i])/2)
            elif i == len(temperature) - 1:
                delT = abs((temperature[i] - temperature[i-1])/2)
            else:
                delT = abs((temperature[i+1]+temperature[i])/2-(temperature[i]+temperature[i-1])/2)
            result += solidFraction[i] * delT
    return result

def getCriteria(ScheilResult,solidusT,liquidusT,numDataThreshold = 10,CSCPoints = [0.4,0.9,0.99], KouPoints = [0.93,0.98], CDPoints = [0.7,0.98]):
    """get different criteria from scheil results

    Args:
        ScheilResult (dict): dict of scheil results 
        solidusT (list): list of solidius T
        liquidusT (list): list of liquidius T
        numDataThreshold (int, optional): min data for crack criteria. Defaults to 10.
        CSCPoints (list, optional): settng for CSC criteria. Defaults to [0.4,0.9,0.99].
        KouPoints (list, optional): settng for Kou criteria. Defaults to [0.93,0.98].
        CDPoints (list, optional): settng for CD criteria. Defaults to [0.7,0.98].

    Returns:
        set: lists of different results crack criteria 
    """

    numFile = len([item for item in ScheilResult.keys()])
    solidFraction = []
    temperature = []
    for index in range(numFile):
        if ScheilResult[f'Point{index}'] != None:
            solidFraction.append([1-item for item in ScheilResult[f'Point{index}']['LIQUID']])
            temperature.append(ScheilResult[f'Point{index}']['TK'])
        else:
            solidFraction.append(None)
            temperature.append(None)

    def getFR(solidT, liquidT):
        FR = []
        for i in range(len(solidT)):
            if liquidT[i] != None and solidT[i] != None:
                FR.append(liquidT[i] - solidT[i])
            else:
                FR.append(None)
        return FR

    def getCSC(temperature, solidFraction, CSCPoints):
        CSCPoints.sort()
        T1 = []
        T2 = []
        T3 = []
        for i in range(len(temperature)):
            Temperature = temperature[i]
            solidFrac = solidFraction[i]

            if Temperature != None and len(Temperature) >= numDataThreshold and max(solidFrac) > max(CSCPoints):
                n_neighbors = 2
                weights = 'distance'
                model = neighbors.KNeighborsRegressor(n_neighbors, weights=weights).fit(np.array(solidFrac).reshape(-1, 1), Temperature)
                T1.append(model.predict(np.array([CSCPoints[2]]).reshape(-1, 1))[0])
                T2.append(model.predict(np.array([CSCPoints[1]]).reshape(-1, 1))[0])
                T3.append(model.predict(np.array([CSCPoints[0]]).reshape(-1, 1))[0])
            else:
                T1.append(None)
                T2.append(None)
                T3.append(None)
        CSC = []
        for i in range(len(T1)):
            if T1[i] != None:
                try:
                    CSC.append((T1[i]-T2[i])/(T2[i]-T3[i]))
                except:
                    CSC.append(None)
            else:
                CSC.append(None)
        return CSC

    def getKou(temperature, solidFraction, KouPoints):
        KouPoints.sort()
        T1 = []
        T2 = []
        solid1 = []
        solid2 = []
        for i in range(len(temperature)):
            Temperature = temperature[i]
            solidFrac = solidFraction[i]
            if Temperature != None and len(Temperature) >= numDataThreshold and max(solidFrac) > max(KouPoints):
                n_neighbors = 2
                weights = 'distance'
                model = neighbors.KNeighborsRegressor(n_neighbors, weights=weights).fit(np.array(solidFrac).reshape(-1, 1), Temperature)
                T1.append(model.predict(np.array([KouPoints[1]]).reshape(-1, 1))[0])
                T2.append(model.predict(np.array([KouPoints[0]]).reshape(-1, 1))[0])
                solid1.append(KouPoints[1])
                solid2.append(KouPoints[0])
            else:
                T1.append(None)
                T2.append(None)
                solid1.append(None)
                solid2.append(None)
        Kou = []
        for i in range(len(T1)):
            if T1[i] != None:
                Kou.append(abs((T1[i]-T2[i])/((solid1[i]**0.5-solid2[i]**0.5))))
            else:
                Kou.append(None)
        return Kou

    def getCD(temperature, solidFraction, CDPoints): #sRDG and iCSC criteria
        CDPoints.sort()
        CD1 = []
        CD2 = []
        fs_0 = CDPoints[0]
        fs_co = CDPoints[1]
        for i in range(len(temperature)):
            Temperature = temperature[i]
            solidFrac = solidFraction[i]
            if Temperature != None and len(Temperature) >= numDataThreshold and max(solidFrac) > max(CDPoints):
                n_neighbors = 2
                weights = 'distance'
                model1 = neighbors.KNeighborsRegressor(n_neighbors, weights=weights).fit(np.array(Temperature).reshape(-1, 1), solidFrac)
                model2 = neighbors.KNeighborsRegressor(n_neighbors, weights=weights).fit(np.array(solidFrac).reshape(-1, 1), Temperature)
                T0 = model2.predict(np.array([fs_0]).reshape(-1, 1))[0]
                Tco = model2.predict(np.array([fs_co]).reshape(-1, 1))[0]
                try:
                    deltT = min((T0 - Tco)/5,0.002)
                    Trange = [item for item in np.arange(Tco,T0+deltT,deltT)]
                    solidFrac = model1.predict(np.array(Trange).reshape(-1, 1))
                    solidFrac = [item for item in solidFrac if item <= fs_co]
                    Trange = [Trange[i] for i in range(len(solidFrac))]
                    CD1.append(getIntegral(Trange, [item**2/(1-item)**2 for item in solidFrac]))
                    CD2.append(getIntegral(Trange, solidFrac))
                except:
                    CD1.append(None)
                    CD2.append(None)
            else:
                CD1.append(None)
                CD2.append(None)
        return CD1, CD2
    
    lackPointList = []
    noResultList = []
    for i in range(len(temperature)):
        if temperature[i] == None:
            noResultList.append(i)
        elif len(temperature[i]) < numDataThreshold:
            lackPointList.append(i)
    print('##############################################################')
    print(f'Warning: Lack of data points (<{numDataThreshold}) in the following samples: ', lackPointList)
    print(f'{len(lackPointList)} samples do not have enough data point.')
    print('Warning: No result in the following samples: ', noResultList)
    print(f'{len(noResultList)} samples have no result.')
    print('##############################################################')

    FR = getFR(solidusT,liquidusT)
    CSC = getCSC(temperature, solidFraction, CSCPoints)
    Kou = getKou(temperature, solidFraction, KouPoints)
    CD1, CD2 = getCD(temperature, solidFraction, CDPoints)

    return FR, CSC, Kou, CD2, CD1
    #CD1: sRDG, CD2: iCSC

def plotHotTearingSusceptibilityMap(path, coord, xComp, yComp, FR, CSC, Kou, iCSC, sRDG):
    """plot the hot crack susceptibility map

    Args:
        path (str): path to open the setting and store the results
        coord (list): list of coordation 
        xComp (str): name for x label
        yComp (str): name for y label
        FR (list): list of values of FR from scheil results
        CSC (list): list of values of CSC from scheil results
        Kou (list): list of values of Kou from scheil results
        iCSC (list): list of values of iCSC from scheil results
        sRDG (list): list of values of sRDG from scheil results

    Returns:
        TIF: store the crack reusults as TIF files
    """
    dotSize = 5
    print('##############################################################')
    print('Plotting Hot Tearing Susceptibility Map...')
    fig = plt.figure(figsize = (8,6),dpi = 400)
    criteria={'FR': FR, 'CSC': CSC, 'Kou': Kou,'iCSC': iCSC,'sRDG':sRDG}
    ii = 1
    for criterion in ['FR', 'CSC', 'Kou','iCSC','sRDG']:
        plt.subplot(2,3,ii,projection='triangular') 
        values = criteria[criterion]
        values_pure = [item for item in values if item != None and not math.isnan(item) and not math.isinf(item)]
        if len(values_pure) == 0:
            print('no data for '+criterion+' from scheil')
            continue;
        top = max(values_pure)
        bottom = min(values_pure)
        norm1 = mpl.colors.Normalize(vmin = bottom, vmax = top)
        for i in range(len(coord)):
            value = values[i]
            (x, y) = coord[i]
            if value != None and not math.isnan(value) and not math.isinf(value):
                RGB1 = mpl.cm.Greys(norm1(value), bytes = True)
                RGB = (RGB1[0]/255,RGB1[1]/255,RGB1[2]/255)
                plt.scatter(x, y,s=dotSize , color = RGB)
            else:
                plt.scatter(x, y,s=dotSize , c = 'yellow')

        cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm1, cmap='Greys'),
            orientation='horizontal',fraction=0.035, pad=0.2,aspect=20)
        for t in cbar.ax.get_xticklabels():
            t.set_fontsize(10)
            t.set_fontname('monoserif')
        # result = np.array(result)
        plt.title(f'{criterion}',fontsize=10, fontname = 'monoserif')
        plt.xlabel(f'W({xComp})',fontsize=10, fontname = 'monoserif')
        plt.ylabel(f'W({yComp})', labelpad=-15,fontsize=10, fontname = 'monoserif')
        plt.xticks([0,0.2,0.4,0.6,0.8,1],[0,0.2,0.4,0.6,0.8,1],fontsize=10, fontname = 'monoserif')
        plt.yticks([0,0.2,0.4,0.6,0.8,1],[0,0.2,0.4,0.6,0.8,1],fontsize=10, fontname = 'monoserif')
        ii += 1
    plt.subplots_adjust(wspace = 0)
    fig.savefig(f'{path}/hotTearing_normal.tif',bbox_inches='tight')
    plt.close()
    print('Plotting Hot Tearing Susceptibility Map...Done')
    return None
