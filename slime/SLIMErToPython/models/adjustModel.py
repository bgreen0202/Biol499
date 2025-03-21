import numpy as np
from SLIMErToPython.models import extractModelStoichMatrix

def adjustModel(model, k, block, scaling):
   
    """
    Adjust the stoichiometry matrix (S) of a metabolic model.

    Parameters:
        model (cobra.Model): The metabolic model.
        k (float): Scaling factor for stoichiometry.
        block (bool): Whether to block exchange reactions for tails and backbones.
        scaling (str): Either 'backbones' or 'tails', indicating which pseudoreaction to rescale.

    Returns:
        cobra.Model: Modified metabolic model.
    """
    # Extract model data
    S, rxnNames, rxns = extractModelStoichMatrix(model)

    # Block exchange of tails and backbones
    if block:
        posT = np.where(np.array(rxnNames) == 'lipid - tails exchange')[0]
        posB = np.where(np.array(rxnNames) == 'lipid - backbones exchange')[0]

        if posT.size > 0:
            rxn = model.reactions.get_by_id(rxns[posT[0]])
            rxn.bounds = (0, 0)

        if posB.size > 0:
            rxn = model.reactions.get_by_id(rxns[posB[0]])
            rxn.bounds = (0, 0)

    # Determine which pseudoreaction to rescale
    if scaling == 'backbones':
        rxnName = 'lipid pseudoreaction - backbone'
    elif scaling == 'tails':
        rxnName = 'lipid pseudoreaction - tail'
    else:
        raise ValueError(f"Invalid scaling type: {scaling}")

    # Find positions
    scaleRxnPos = np.where(np.array(rxnNames) == rxnName)[0]
    if scaleRxnPos.size == 0:
        raise ValueError(f"Reaction '{rxnName}' not found in model.")

    scaleRxnPos = scaleRxnPos[0]

    # Find metabolites associated with the pseudoreaction (negative stoichiometry = reactants)
    scaleMets = np.where(S[:, scaleRxnPos] < 0)[0]

    # Rescale stoichiometric coefficients
    for metPos in scaleMets:
        S[metPos, scaleRxnPos] *= k

    return model  # The model itself is modified, but S is not automatically updated in COBRApy
