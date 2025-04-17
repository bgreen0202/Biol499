from cobra import Metabolite, Reaction
from SLIMErToPython.models import getNewIndex

def addLipidSpecies(model, metName, metFormula, exchange):
    """
    Simplified version of addLipidSpecies that ignores compartments.
    Adds lipid species and optionally creates an exchange reaction.
    """

    newID = getNewIndex(model.metabolites)
    metID = f's_{newID}'

    new_met = Metabolite(
        id=metID,
        name=metName,
        formula=metFormula
    )
    model.add_metabolites([new_met])

    # Add exchange reaction if needed
    if exchange:
        newID = getNewIndex(model.reactions)
        rxnName = f'{metName} exchange'

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

