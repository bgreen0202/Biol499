import numpy as np
from scipy.optimize import minimize
from cobra import Model

from SLIMErToPython.models import adjustModel

def scaleAbundancesInModel(model, data, scaling):
    """
    Direct translation of scaleAbundancesInModel from MATLAB to Python (Benjamín J. Sánchez, 2018-09-04).

    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """
    # Initial guess for k
    k0 = 1

    # Find optimal scaling factor using fminsearch equivalent (scipy's minimize)
    result = minimize(lambda k: unusedLipid(k[0], model, data, scaling), [k0], method='Nelder-Mead')
    #Find the optimal scaling factor kOpt. This ensures minimal unused lipid accumulation in the model.
    kOpt = result.x[0]

    # Find optimality range
    krange = [0, 0]
    result1 = minimize(lambda k: +minScaling(k[0], model, data, scaling, kOpt), [kOpt], method='Nelder-Mead')
    krange[0] = result1.x[0]
    
    result2 = minimize(lambda k: -minScaling(k[0], model, data, scaling, kOpt), [kOpt], method='Nelder-Mead')
    krange[1] = result2.x[0]

    print(f"Optimality range: k = [ {krange[0]} , {krange[1]} ]")

    if krange[0] == krange[1]:
        raise ValueError('Could not find an optimality range!')

    # Scale using the average of the optimality range
    k = np.mean(krange)
    model = adjustModel(model, k, True, scaling)
    print(f"Scaled {scaling[:-1]} data in model: k = {k}")

    return model, k


def unusedLipid(k, model, data, scaling):
    """
    Objective function to minimize unused lipid abundance.
    """

    # Adjust stoichiometry using scaling factor
    model = adjustModel(model, k, False, scaling)

    # Simulate growth
    try:
        sol = simulateGrowth(model, data['fluxData'])
    except:
        sol = {'x': np.ones(len(model.reactions))}

    # Objective function: unused tails or backbones
    exchange_tails = getReactionFlux(model, sol, 'lipid - tails exchange')
    exchange_backs = getReactionFlux(model, sol, 'lipid - backbones exchange')
    exchange = exchange_tails + exchange_backs

    print(f"Scaling abundance data: k = {k} -> exchange = {exchange}")

    return exchange


def minScaling(k, model, data, scaling, kOpt):
    """
    Objective function to find the minimum k within feasible range.
    """

    # Adjust stoichiometry with k
    model = adjustModel(model, k, True, scaling)

    # Simulate growth and get NGAM 
    try:
        sol = simulateGrowth(model, data['fluxData'])
        posNGAM = getReactionIndex(model, 'non-growth associated maintenance reaction')
        maintenance = sol['x'][posNGAM]
        print(f"Finding scaling range: k = {k} -> Maintenance = {maintenance}")
    except:
        print(f"Finding scaling range: k = {k} -> Maintenance = 0")
        k = kOpt  # Any unfeasible simulation returns the original value

    return k


def getReactionFlux(model, sol, rxn_name):
    """
    Helper to fetch the flux of a reaction by its name.
    """
    try:
        rxn_index = getReactionIndex(model, rxn_name)
        return sol['x'][rxn_index]
    except:
        return 0


def getReactionIndex(model, rxn_name):
    """
    Helper to get the index of a reaction by name.
    """
    for i, rxn in enumerate(model.reactions):
        if rxn.name == rxn_name:
            return i
    raise ValueError(f"Reaction '{rxn_name}' not found in model.")


def simulateGrowth(model, fluxData):
    solution = model.optimize()  
    return {'x': solution.fluxes.values} 
    pass
