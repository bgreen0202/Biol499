import numpy as np
from SLIMErToPython.models import extractModelStoichMatrix

def sumBioMass(model, comps):
    """
    Calculates the breakdown of biomass into fractions.

    Parameters:
        model (cobra.Model): The metabolic model.
        comps (list): List of tuples, where each tuple contains:
            (metabolite ID, molecular weight, component type).

    Returns:
        tuple: (X, P, C, R, D, L)
            X: Biomass fraction without lipids [g/gDW]
            P: Protein fraction [g/gDW]
            C: Carbohydrate fraction [g/gDW]
            R: RNA fraction [g/gDW]
            D: DNA fraction [g/gDW]
            L: Lipid fraction [g/gDW]
    """

    # Extract stoichiometric matrix
    S, met_ids, rxn_ids = extractModelStoichMatrix(model)  

    # Initialize biomass fractions
    X = 0

    # Get each fraction
    P, X = getFraction(model, comps, 'P', X, S, met_ids, rxn_ids, water_correction=True)
    C, X = getFraction(model, comps, 'C', X, S, met_ids, rxn_ids, water_correction=True)
    R, X = getFraction(model, comps, 'R', X, S, met_ids, rxn_ids, water_correction=True)
    D, X = getFraction(model, comps, 'D', X, S, met_ids, rxn_ids, water_correction=True)
    L, X = getFraction(model, comps, 'L', X, S, met_ids, rxn_ids, water_correction=False)  

    # Find biomass reaction index
    try:
        bio_index = rxn_ids.index('r_4041')  # Need our biomass reaction ID
    except ValueError:
        raise ValueError("Biomass reaction 'r_4041' not found in model.")

    # Add up any remaining components
    for i, met in enumerate(met_ids):
        matching_comp = next((c for c in comps if c[0] == met), None)
        if matching_comp:
            MW = matching_comp[1]
            abundance = -S[i, bio_index] * MW / 1000  # Convert to g/gDW
            X += abundance

    return X, P, C, R, D, L


def getFraction(model, comps, compType, X, S, met_ids, rxn_ids, water_correction):
    """
    Helper function to calculate fraction of a component type.
    """

    # Define pseudoreaction name
    rxnNames = {
        'P': 'protein pseudoreaction',
        'C': 'carbohydrate pseudoreaction',
        'R': 'RNA pseudoreaction',
        'D': 'DNA pseudoreaction',
        'L': 'lipid pseudoreaction'
    }
    
    rxnName = rxnNames.get(compType)
    if not rxnName:
        raise ValueError(f"Invalid component type: {compType}")

    # Find reaction index
    try:
        fraction_index = rxn_ids.index(rxnName)
    except ValueError:
        raise ValueError(f"Reaction '{rxnName}' not found in model.")

    # Add up fraction
    F = 0
    if compType == 'L':
        # Lipid backbone special case
        try:
            lipid_index = rxn_ids.index('lipid pseudoreaction - backbone')
        except ValueError:
            raise ValueError("Reaction 'lipid pseudoreaction - backbone' not found in model.")

        # Sum up negative stoichiometry coefficients
        F = -np.sum(abs(S[:, lipid_index][S[:, lipid_index] < 0]))  # g/gDW
    else:
        # Other biomolecules
        componentMets = [c for c in comps if c[2] == compType]

        for i, met in enumerate(met_ids):
            matchingComp = next((c for c in componentMets if c[0] == met), None)
            if matchingComp:
                MW = matchingComp[1]
                MW_adjusted = MW - 18 if water_correction else MW  # Apply water loss correction
                abundance = -S[i, fraction_index] * MW_adjusted / 1000  # Convert to g/gDW
                F += abundance

    X += F
    return F, X
