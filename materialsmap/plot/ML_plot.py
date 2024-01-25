import os
from materialsmap.core.GenerateEqScript import getSettings
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import json
import numpy as np

def readResult(folder_ML,properties):
    """read the previously saved ML result

    Args:
        folder_ML (str): path to folder that contains ML results
        properties (str): properties of interested

    Returns:
        dcit: dict of ML results
    """
    f_ML = open(f'{folder_ML}/{properties}/Result/data_ML.json')
    MLResult = json.load(f_ML)
    return MLResult

def getMLResult(result_ML,properties):
    """_summary_

    Args:
        result_ML (dict): ML results in dict format
        properties (str): properties of interested

    Returns:
        list: list of value of properties
    """
    properties_val = []
    if properties == 'melting_temperature':
        key = 'TK'
    else:
        raise "Properties not implement yet"
    for val in result_ML.values():
        properties_val.append(val[key])
    return properties_val


def ML_plot(path,properties='melting_temperature'):
    """_summary_

    Args:
        path (str): path to open the setting and store the results
        properties (str): properties of interested. Defaults to melting_temperature.
    """
    path = os.path.abspath(path)
    settings = np.load(f'{path}/setting.npy',allow_pickle=True) #settings = [TemperatureRange,numPoint,numSimultion,relatedEles,terminalAlloys,indep_terminalAlloys,database,pressure]
    xComp = settings[5][0]
    yComp = settings[5][1]
    comps = settings[4]
    numFile = int(settings[1])
    data = pd.read_excel(f'{path}/composition_for_feasibilityMap.xlsx')
    composition_data = dict()
    
    for item in data.columns:
        if 'alloy' in item:
            composition_data[item[6:]] = data[item].values.tolist()[:numFile]
    composition_data = pd.DataFrame(composition_data)
    coord = []
    for index in range(len(composition_data)):
        x = composition_data[xComp].values[index]
        y = composition_data[yComp].values[index]
        coord.append((x, y))
    folder_ML = path + '/ML'
    result_ML = readResult(folder_ML,properties)
    properties_val = getMLResult(result_ML,properties)
    # plot figure
    dotSize = 12
    plt.figure(figsize=(4, 4), dpi=400)
    subFig = plt.subplot(projection='triangular')
    top = max(properties_val)
    bottom = min(properties_val)
    norm = mpl.colors.Normalize(vmin = bottom, vmax = top)
    cmap = 'coolwarm'
    for i in range(len(coord)):
        value = properties_val[i]
        (x,y) = coord[i]
        if value != None:
            RGB1 = mpl.cm.coolwarm(norm(value), bytes = True)
            RGB = (RGB1[0]/255,RGB1[1]/255,RGB1[2]/255)
        else:
            RGB = 'yellow'
        subFig.scatter(x, y,s=dotSize, color = RGB, marker = "o")
    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    orientation='horizontal',fraction=0.035, pad=0.2,aspect=20)
    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(10)
        t.set_fontname('times')
    plt.title(f'{properties} from ML\n{round(bottom,3)}-{round(top,3)}',fontsize = 10, fontname = 'times')
    plt.xlabel(f'X({xComp})',fontsize = 10, fontname = 'times')
    plt.ylabel(f'X({yComp})', labelpad=-40,fontsize = 10, fontname = 'times')

    plt.tick_params(labelsize=10)
    plt.savefig(f'{path}/ML_results.tif', bbox_inches='tight')