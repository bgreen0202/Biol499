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
    metabolite_lookup = {met.name: met for met in model.metabolites}  # only one pass through model.metabolites

    tailIDs = []
    for metName in chainData['metNames']:
        tailName = metName
        if tailName in metabolite_lookup:
            tailIDs.append(metabolite_lookup[tailName])
        else:
            raise ValueError(f"Metabolite {tailName} not found in model.")
    
    # Generate new reaction ID and name
    newID = getNewIndex([r.id for r in model.reactions])
    rxnID = 'r_' + newID
    rxnName = 'lipid pseudoreaction - tail'
    
    # Find the main 'lipid - tails' metabolite
    if 'lipid - tails' not in metabolite_lookup:
        raise ValueError("Metabolite 'lipid - tails' not found in model.")
    tailMet = metabolite_lookup['lipid - tails']
    
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
    metabolites_dict[tailMet] = 1  # This ensures that lipid - tails is produced with a coefficient of +1.
    
    reaction.add_metabolites(metabolites_dict)
    model.add_reactions([reaction])

    return model
