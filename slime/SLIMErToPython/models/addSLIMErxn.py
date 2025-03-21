import numpy as np
import cobra
from cobra import Model, Reaction, Metabolite
from SLIMErToPython.models import getStoichFromFormula, getBackboneName, getMWfromFormula, getNewIndex

def addSLIMErxn(model, rxnID, specID):
    """
    Python translation of addSLIMErxn from MATLAB to COBRApy.
    Adds SLIME reactions to model.
    """
    toDelete = False
    
    if specID == '':  # Case 1: Existing ISA reaction
        reaction = model.reactions.get_by_id(rxnID)
        specID = list(reaction.metabolites.keys())[0].id
    
    if rxnID == '':  # Case 2: New ISA reaction
        rxnID = 'r_' + getNewIndex([r.id for r in model.reactions])
    
    specMet = model.metabolites.get_by_id(specID)
    specName = specMet.name
    rxnName = specName + ' SLIME rxn'
    
    backName = getBackboneName(specName)
    if backName is None:
        print('Removing unrecognized ISA reaction')
        return model, True
    
    backMet = next((m for m in model.metabolites if m.name.startswith(backName + ' [')), None)
    if not backMet:
        print(f'Backbone metabolite {backName} not found')
        return model, True
    
    specFormula = specMet.formula
    specMW = getMWfromFormula(specFormula)
    
    # Determine fatty acid tails for different lipid classes
    tailsRxn = []
    if backName in ['phosphatidylserine', 'phosphatidylethanolamine', 'phosphatidate', 'phosphatidylglycerol', 'cardiolipin']:
        tailsRxn = ['C' + part for part in specName.split('(')[-1].replace(')', '').split(', ')]

    # Need to decide our fatty acid names, ie: palmitate C16:0? myristoyl acid C14:0? iso-hepta-1-methyl-decanoyl as C15:0 ?? 
    # Recall that according to FAMES on BG8 has [C10:0, C12:0, C14:0, C15:0. C16:0, C16:ln9?, C16:ln7?, C16:ln6?, C16:ln5?, C16unknown1, C16:2, C16:3, C18:0, C18:ln9, C18:ln7 ]
    elif backName in ['fatty acid', 'ergosterol ester', 'long-chain base', 'long-chain base phosphate']:
        species_map = {
            'myristate': 'C14:0', 'palmitate': 'C16:0', 'stearate': 'C18:0',
            'oleate': 'C18:1', 'sphinganine': 'C18:0', 'phytosphingosine': 'C18:0'
        }
        tailsRxn = [species_map.get(specName)]
    
    tailIDs, tailsFormulas, tailsMWs = [], [], []
    for met in model.metabolites:
        if ' chain [cytoplasm]' in met.name:
            tailIDs.append(met.id)
            tailsFormulas.append(met.formula)
            tailsMWs.append(getMWfromFormula(met.formula))
    
    tailCoeffs = np.zeros(len(tailIDs))
    prodFormulas = []
    
    for tail in tailsRxn:
        tailName = tail + ' chain [cytoplasm]'
        try:
            idx = next(i for i, m in enumerate(model.metabolites) if m.name == tailName)
            tailCoeffs[idx] += 1
            prodFormulas.append(tailsFormulas[idx])
        except StopIteration:
            raise ValueError(f"Tail metabolite '{tailName}' not found in the model.")
    
    # Assign molecular formula to backbone (balance SLIME rxn)
    backFormula = ''
    for elem in ['C', 'H', 'N', 'O', 'P', 'S']:
        Nin = getStoichFromFormula(specFormula, elem)
        Nout = sum(getStoichFromFormula(f, elem) for f in prodFormulas)
        diff = Nin - Nout
        if diff > 0:
            backFormula += f'{elem}{diff}' if diff > 1 else elem
    
    backMet.formula = backFormula  # Assign to backbone metabolite
    
    # Add SLIME reaction
    reaction = Reaction(rxnID)
    reaction.name = rxnName
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    
    reaction.add_metabolites({
        specMet: -1,
        backMet:specMW
    })
    
    for tailID, coeff in zip(tailIDs, tailCoeffs):
        if coeff != 0:
            tailMet = model.metabolites.get_by_id(tailID)
            reaction.add_metabolites({tailMet: coeff})
    
    model.add_reactions([reaction])
    
    #below calculates confidence scores but ??
    try:
        print(f"Added reaction: {reaction.build_reaction_string(use_metabolite_names=True)}")
        if hasattr(model, 'rxnConfidenceScores'):
            model.rxnConfidenceScores[rxnID] = 1
    except Exception as e:
        print(f"Repeated: {rxnName} ({e})")
    
    return model, toDelete
