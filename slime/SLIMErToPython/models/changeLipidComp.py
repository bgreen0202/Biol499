from SLIMErToPython.models import getNewIndex
import cobra 

def changeLipidComp(model, lipidData):
    """
    Direct translation of changeLipidComp from MATLAB to Python (Benjamín J. Sánchez, 2018-05-20).
    Changes stoichiometric data for backbone

    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """

    # Detect if the model has a metabolite for backbones
    metPos = [i for i, name in enumerate(model.metNames) if name == 'lipid - backbones [cytoplasm]']
    if len(metPos) > 0:
        # Create new pseudo-reaction for backbones
        newID = getNewIndex(model.rxns) 
        rxnID = 'r_' + newID
        rxnName = 'lipid pseudoreaction - backbone'
        metID = model.mets[metPos[0]]
    else:
        # Modify normal lipid pseudo-reaction --> OBVIOUSLY NEED TO MAKE ONE FOR US WHAT ID?
        rxnID = 'r_2108'
        rxnName = 'lipid pseudoreaction'
        metID = 's_1096[c]'

    # Create lipid pseudo-reaction for backbones (or modify normal one)
    metaboliteList = lipidData['metIDs'] + [metID]
    stoichCoeffList = [-x for x in lipidData['abundance']] + [1]

    model.add_reaction(model, rxnID,
                        reactionName=rxnName,
                        metaboliteList=metaboliteList,
                        stoichCoeffList=stoichCoeffList,
                        reversible=False,
                        lowerBound=0,
                        upperBound=1000)

    # Print reaction formula
    #printRxnFormula(model, rxnID, True, True, True)

    # Set confidence score
    if hasattr(model, 'rxnConfidenceScores'):
        rxnIndex = model.rxns.index(rxnID)
        model.rxnConfidenceScores[rxnIndex] = 1

    return model
