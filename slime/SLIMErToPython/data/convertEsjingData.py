import numpy as np
import pandas as pd
from cobra import Model

from slime.SLIMErToPython.data import matchToModel
from slime.SLIMErToPython.models import getBackboneName, getMWfromFormula

def convertEsjingData(data_old, model, condense):

    """
    Reads data from a CSV, extracts name, abundunce, st. dev and removes lipid that contain 10, 12, 14, 20, 22
    carbon chains OR contain :2 in their names OR don't belong to lipid classes, then creates a dictionary.
    """

    MWs = np.zeros(len(data_old['metNames']))
    backbones = [None] * len(data_old['metNames'])
    
    for i, met_name in enumerate(data_old['metNames']):
        pos = matchToModel(model, met_name)
        if np.sum(pos) > 0:
            MWi = getMWfromFormula([model.metabolites[i].formula for i in pos])
            MWs[i] = np.mean(MWi)
            backbone = model.metabolites[pos[0]].name
            backbones[i] = getBackboneName(backbone).split('[')[0].strip()
    
    mask = MWs > 0
    data_old = {key: np.array(value)[mask] for key, value in data_old.items()}
    data_old['abundance'] /= 100
    data_old['std'] /= 100
    MWs = MWs[mask] * 1000
    
    if condense:
        data = {'lipidData': change_units(data_old, MWs)}
        data['lipidData'] = squash_lipids(data['lipidData'], backbones)
        
        chain_data = pd.read_csv('chainData_Lahtvee2017.csv')
        chain_names = chain_data.iloc[:, 0].str.replace('C', '').str.replace(' chain', '')
        chain_formulas = chain_data.iloc[:, 1]
        chain_MWs = getMWfromFormula(chain_formulas.tolist())
        
        data['chainData'] = squash_lipids(data_old, chain_names.tolist())
        data['chainData']['formulas'] = chain_formulas.tolist()
        data['chainData']['metNames'] = chain_names.tolist()
        data['chainData'] = change_units(data['chainData'], chain_MWs)
    else:
        data = change_units(data_old, MWs)
        data['MWs'] = MWs
    
    data['otherData'] = {
        'metIDs': ['protein', 'RNA'],
        'abundance': [0.5, 0.06]
    }
    data['fluxData'] = {
        'rxnIDs': ['r_1714', 'r_2111'],
        'averages': [-20.4, 0.41],
        'stdevs': [20.4 * 0.01, 0.41 * 0.01]
    }
    
    return data

def change_units(data, MWs):
    lipid_content = 0.08 * np.sum(data['abundance'])
    data['std'] *= MWs
    data['abundance'] *= MWs
    data['std'] /= np.sum(data['abundance'])
    data['abundance'] /= np.sum(data['abundance'])
    data['std'] *= lipid_content
    data['abundance'] *= lipid_content
    return data

def squash_lipids(data_old, met_names):
    unique_met_names = list(set(met_names))
    abundance = np.zeros(len(unique_met_names))
    std = np.zeros(len(unique_met_names))
    
    is_tail = len(unique_met_names) == len(met_names)
    
    for i, met_name in enumerate(unique_met_names):
        if is_tail:
            hits = np.array([met_name in name for name in data_old['metNames']])
        else:
            hits = np.array(met_names) == met_name
        
        abundance[i] = np.sum(data_old['abundance'] * hits)
        std[i] = np.mean(data_old['std'] * hits)
    
    return {'metNames': unique_met_names, 'abundance': abundance, 'std': std}
