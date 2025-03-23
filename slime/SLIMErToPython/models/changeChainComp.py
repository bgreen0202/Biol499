import numpy as np
import cobra
from SLIMErToPython.models import getNewIndex

def changeChainComp(model, chainData):
    """
    Updates the metabolic model by adding a lipid pseudoreaction for tails.
    
    Parameters:
        model (cobra.Model): The metabolic model.
        chainData (dict): Dictionary with keys:
            'metNames': List of tail metabolite base names (without compartment suffix).
            'abundance': List of abundances for each tail metabolite.
    

    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """
    # Identify metabolite IDs for the tails
    tailIDs = []
    for metName in chainData['metNames']:
        tailName = metName + ' [cytoplasm]'
        try:
            tailMet = next(met for met in model.metabolites if met.name == tailName)
            tailIDs.append(tailMet)
        except StopIteration:
            raise ValueError(f"Metabolite {tailName} not found in model.")
    
    # Generate new reaction ID and name
    newID = getNewIndex([r.id for r in model.reactions])
    rxnID = 'r_' + newID
    rxnName = 'lipid pseudoreaction - tail'
    
    # Find the main 'lipid - tails' metabolite
    try:
        tailMet = next(met for met in model.metabolites if met.name == 'lipid - tails [cytoplasm]')
    except StopIteration:
        raise ValueError("Metabolite 'lipid - tails [cytoplasm]' not found in model.")
    
    # Create the new reaction
    reaction = cobra.Reaction(rxnID)
    reaction.name = rxnName
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    
    # Define metabolite stoichiometry
    # IE: tailIDs = ['C16:0', 'C18:0']
    # chainData['abundance'] = [0.3, 0.7]
    # [('C16:0', 0.3), ('C18:0', 0.7)]
    metabolites_dict = {met: -abundance for met, abundance in zip(tailIDs, chainData['abundance'])}
    metabolites_dict[tailMet] = 1  #T his ensures that lipid - tails is produced with a coefficient of +1.
    
    reaction.add_metabolites(metabolites_dict)
    model.add_reactions([reaction])

    #printRxnFormula(model, rxnID, True, True, True)

    # Set confidence score 
    if hasattr(model, 'rxnConfidenceScores'):
        rxnIndex = model.rxns.index(rxnID)
        model.rxnConfidenceScores[rxnIndex] = 1

    return model
