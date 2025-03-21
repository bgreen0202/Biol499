import pandas as pd

def readLahtveeData(i):
    """
    Reads data from a CSV and stores it in a dict.
    """
    data = {}
    
    # Lipid data
    lipid_df = pd.read_csv('lipidData_Lahtvee2017.csv')
    data['lipidData'] = {
        'metAbbrev': lipid_df.iloc[:, 0].tolist(),
        'metNames': lipid_df.iloc[:, 1].tolist(),
        'metIDs': lipid_df.iloc[:, 2].tolist(),
        'abundance': lipid_df.iloc[:, 2 + i].tolist()
    }
    
    # Chain data
    chain_df = pd.read_csv('chainData_Lahtvee2017.csv')
    data['chainData'] = {
        'metNames': chain_df.iloc[:, 0].tolist(),
        'formulas': chain_df.iloc[:, 1].tolist(),
        'abundance': chain_df.iloc[:, 1 + 2 * i].tolist(),
        'std': chain_df.iloc[:, 2 + 2 * i].tolist()
    }
    
    # Other composition data
    comp_df = pd.read_csv('compData_Lahtvee2017.csv')
    data['otherData'] = {
        'metIDs': comp_df.iloc[:, 1].tolist(),
        'abundance': comp_df.iloc[:, 1 + i].tolist()
    }
    
    # Flux data
    flux_df = pd.read_csv('fluxData_Lahtvee2017.csv')
    data['fluxData'] = {
        'rxnIDs': flux_df.iloc[:, 1].tolist(),
        'averages': flux_df.iloc[:, 2 * i].tolist(),
        'stdevs': flux_df.iloc[:, 2 * i + 1].tolist()
    }
    
    return data
