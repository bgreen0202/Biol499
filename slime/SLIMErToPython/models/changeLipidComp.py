from cobra import Reaction
from SLIMErToPython.models import getNewIndex

def changeLipidComp(model, lipidData):
    """
    Python version of changeLipidComp that updates the lipid backbone pseudo-reaction.
    Uses classic for-loops for clarity, no zip or comprehension.
    """

    # Detect if the model has a metabolite for backbones
    metID = None
    for met in model.metabolites:
        if met.name == 'lipid - backbones':
            metID = met.id
            break

    if metID:
        # Create new pseudo-reaction for backbones
        newID = getNewIndex(model.reactions)
        rxnID = f'r_{newID}'
        rxnName = 'lipid pseudoreaction - backbone'
    else:
        # Use fallback 
        rxnID = 'r_2108'
        rxnName = 'lipid pseudoreaction'
        metID = 's_1096[c]'  # Assume this ID exists in the model

    # Create lipid pseudo-reaction for backbones (or modify normal one)
    metaboliteList = lipidData['metIDs'] + [metID]
    stoichCoeffList = [-x for x in lipidData['abundance']] + [1]

    stoichDict = {}
    for i in range(len(metaboliteList)):
        met = model.metabolites.get_by_id(metaboliteList[i])
        coeff = stoichCoeffList[i]
        stoichDict[met] = coeff

    # Create and add the reaction
    rxn = Reaction(id=rxnID, name=rxnName, lower_bound=0, upper_bound=1000)
    rxn.add_metabolites(stoichDict)
    model.add_reactions([rxn])

    return model

