from cobra import Metabolite, Reaction
from SLIMErToPython.models import getNewIndex

def addLipidSpecies(model, metName, metFormula, exchange):
    """
    Direct translation of the MATLAB function addLipidSpecies into Python using COBRApy.
    Add lipid species and exchange reactions if needed to the model.

    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """
    newID = getNewIndex(model.metabolites)

    compName = metName[metName.find('[') + 1 : metName.find(']')]

    compPos = model.compartments.get(compName, None)
    if compPos is None:
        raise ValueError(f"Compartment '{compName}' not found in model.")

    metID = f's_{newID}[{compPos}]'

    new_met = Metabolite(
        id=metID,
        name=metName,
        formula=metFormula,
        compartment=compPos
    )
    model.add_metabolites([new_met])

    # Add exchange reaction if needed
    if exchange:
        newID = getNewIndex(model.reactions)
        baseMetName = metName[:metName.find('[')].strip()
        rxnName = f'{baseMetName} exchange'

        exchange_rxn = Reaction(
            id=f'r_{newID}',
            name=rxnName,
            lower_bound=0,     
            upper_bound=1000
        )
        exchange_rxn.add_metabolites({new_met: -1})

        model.add_reactions([exchange_rxn])

        print(f'{exchange_rxn.id}: {rxnName} - {metID} ->')

    return model
