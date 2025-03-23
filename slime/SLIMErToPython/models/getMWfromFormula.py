import numpy as np
from SLIMErToPython.models import getStoichFromFormula

def getMWfromFormula(metFormulas):
    """
    Direct translation of getMWfromFormula from MATLAB to Python (Benjamín J. Sánchez, 2018-09-04).
    Gets MW based on checical formula in g/mmol.

    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """

    # Atomic weights (g/mol, same order as elements)
    elements = ['C', 'H', 'O', 'N', 'P', 'S']
    atomic_weights = np.array([12.011, 1.00794, 15.9994, 14.00674, 30.973762, 32.066])

    MWs = np.zeros(len(metFormulas))

    # For each element, calculate contribution to MW
    for idx, element in enumerate(elements):
        # Use getStoichFromFormula to get count of each element
        stoich = np.array([getStoichFromFormula(formula, element) for formula in metFormulas])

        # Add contribution to MW (convert from mg to g if needed, hence /1000 if necessary)
        MWs += stoich * atomic_weights[idx] / 1000  # g/mmol (g per mmol)

    return MWs
