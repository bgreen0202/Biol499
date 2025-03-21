import re
import numpy as np

def getStoichFromFormula(metFormulas, element):
    """
    Direct translation of getStoichFromFormula from MATLAB to Python (Benjamín J. Sánchez, 2018-09-04).
    Gets stoich coeff for a certain element in all metabolites formulas 
    ie: metFormulas = ['C6H12O6', 'H2O', 'C2H5OH'] element = "C", printed output = [6. 0. 2.]
    """
    stoich = np.zeros(len(metFormulas), dtype=float)

    element_pattern = re.compile(rf'{element}(\d*)')

    for i, formula in enumerate(metFormulas):
        matches = element_pattern.findall(formula)

        if len(matches) > 1:
            raise ValueError(f'Non-standard formula: {formula}')

        elif len(matches) == 0:
            stoich[i] = 0

        elif matches[0] == '':
            # Element exists with no number, e.g., "C" counts as 1
            stoich[i] = 1

        else:
            stoich[i] = int(matches[0])

    return stoich
