import numpy as np
from cobra import Model, Reaction, Metabolite

from SLIMErToPython.models import getNewIndex
from SLIMErToPython.models import addLipidSpecies, addSLIMErxn, changeChainComp, changeLipidComp, getBackboneName
from SLIMErToPython.models import extractModelStoichMatrix

def SLIMEr(model, data, includeTails):
    """
    Python translation of SLIMEr.
    Adapted from the MATLAB function by Benjamín J. Sánchez (2018-05-20).
    """

    # Add new backbones
    for met in model.metabolites:
        backName = getBackboneName(met.name)
        if backName:
            existing_met = model.metabolites.get_by_id(backName) if backName in model.metabolites else None
            
            model = addLipidSpecies(model, backName, '', False)

            # Add transport reaction to cytoplasm for non-cytoplasmic backbones
            if 'cytoplasm' not in backName:
                cytoName = f"{backName.split('[')[0]}[cytoplasm]"
                existing_cyto_met = model.metabolites.get_by_id(cytoName) if cytoName in model.metabolites else None

                if existing_cyto_met:
                    existing_cyto_met.formula = ''
                else:
                    model = addLipidSpecies(model, cytoName, '', False)

                # Avoid duplicating transport reactions
                transport_rxn_id = f'r_{getNewIndex(model.reactions)}'
                if transport_rxn_id not in model.reactions:
                    trans_rxn = Reaction(transport_rxn_id)
                    trans_rxn.name = f"{backName.split('[')[0]} transport"
                    trans_rxn.lower_bound = 0
                    trans_rxn.upper_bound = 1000
                    trans_rxn.add_metabolites({existing_met: -1, existing_cyto_met: 1})
                    model.add_reactions([trans_rxn])

    # Get backbone metabolite IDs
    metIDs = []
    for metName in data['lipidData']['metNames']:
        cytoName = f"{metName} [cytoplasm]"
        met = model.metabolites.get_by_id(cytoName) if cytoName in model.metabolites else None
        if met:
            metIDs.append(met.id)
    data['lipidData']['metIDs'] = metIDs

    # Add chains
    for i, metName in enumerate(data['chainData']['metNames']):
        fullName = f"{metName} [cytoplasm]"
        formula = data['chainData']['formulas'][i]
        model = addLipidSpecies(model, fullName, formula, not includeTails)

    # Create lipid pseudo-reaction for backbones
    model = addLipidSpecies(model, 'lipid - backbones [cytoplasm]', '', includeTails)
    model = changeLipidComp(model, data['lipidData'])

    # Create lipid pseudo-reaction for tails
    if includeTails:
        model = addLipidSpecies(model, 'lipid - tails [cytoplasm]', '', includeTails)
        model = changeChainComp(model, data['chainData'])

    # Replace ISA reactions with SLIME reactions
    toDelete = []
    for rxn in model.reactions:
        if 'isa ' in rxn.name:
            result = addSLIMErxn(model, rxn.id, None)
            if result is not None:
                model, delete = result
                if delete:
                    toDelete.append(rxn)
        elif rxn.name == 'complex sphingolipid transport':
            toDelete.append(rxn)

    # Add new SLIME reactions for backbones
    for met in model.metabolites:
        backName = getBackboneName(met.name)
        if backName:
            model, _ = addSLIMErxn(model, None, met.id)

    # Change overall lipid pseudo-reaction
    back_met = model.metabolites.get_by_id('lipid - backbones [cytoplasm]')
    lipid_met = model.metabolites.get_by_id('lipid [cytoplasm]')

    if includeTails:
        tail_met = model.metabolites.get_by_id('lipid - tails [cytoplasm]')
        mets = [tail_met, back_met, lipid_met]
        coeffs = [-1, -1, 1]
    else:
        mets = [back_met, lipid_met]
        coeffs = [-1, 1]

    lipid_rxn = Reaction('r_2108')
    lipid_rxn.name = 'lipid pseudoreaction - merge'
    lipid_rxn.lower_bound = 0
    lipid_rxn.upper_bound = 1000
    lipid_rxn.add_metabolites({mets[i]: coeffs[i] for i in range(len(mets))})
    model.add_reactions([lipid_rxn])

    # Remove unused ISA reactions properly
    model.remove_reactions(toDelete)

    # Remove GAM requirement
    GAM = 0.0
    bio_rxn = model.reactions.get_by_id('biomass pseudoreaction')
    ATP = model.metabolites.get_by_id('ATP [cytoplasm]')
    H2O = model.metabolites.get_by_id('H2O [cytoplasm]')
    ADP = model.metabolites.get_by_id('ADP [cytoplasm]')
    H = model.metabolites.get_by_id('H+ [cytoplasm]')
    P = model.metabolites.get_by_id('phosphate [cytoplasm]')

    bio_rxn.add_metabolites({ATP: -GAM, H2O: -GAM, ADP: GAM, H: GAM, P: GAM})

    S = extractModelStoichMatrix(model)

    return model, S 
