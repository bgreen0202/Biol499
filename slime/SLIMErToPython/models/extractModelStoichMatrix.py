import numpy as np
from cobra.util.array import create_stoichiometric_matrix

def extractModelStoichMatrix(model):
    """
    Extracts stoichiometric matrix (S), reaction names (rxnNames), and reaction IDs (rxns) from a COBRApy model.

    Returns:
        tuple: (S, rxnNames, rxns) where:
            - S is the stoichiometric matrix (numpy array).
            - rxnNames is a list of reaction names.
            - rxns is a list of reaction IDs.

    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """
    S = create_stoichiometric_matrix(model, array_type='dense')
    rxnNames = [rxn.name for rxn in model.reactions]  
    rxns = [rxn.id for rxn in model.reactions] 
    
    return S, rxnNames, rxns