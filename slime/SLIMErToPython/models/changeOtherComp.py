from SLIMErToPython.models import sumBioMass


def changeOtherComp(model, data):
    """
    Direct translation of changeOtherComp from MATLAB to Python (Benjamín J. Sánchez, 2018-09-04).
    Changes composition for other stuff in BIOMASS objective function

    most data we need to change to reflect ours
    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """

    otherData = data['otherData']

    # Define the components list for biomass NEED THIS FOR OURS
    comps = [
        ('s_0404[c]', 89.09, 'P'), ('s_0542[c]', 121.16, 'P'), ('s_0432[c]', 133.11, 'P'),
        ('s_0748[c]', 147.13, 'P'), ('s_1314[c]', 165.19, 'P'), ('s_0757[c]', 75.07, 'P'),
        ('s_0832[c]', 155.15, 'P'), ('s_0847[c]', 131.17, 'P'), ('s_1099[c]', 146.19, 'P'),
        ('s_1077[c]', 131.17, 'P'), ('s_1148[c]', 149.21, 'P'), ('s_0430[c]', 132.12, 'P'),
        ('s_1379[c]', 115.13, 'P'), ('s_0747[c]', 146.14, 'P'), ('s_0428[c]', 174.2, 'P'),
        ('s_1428[c]', 105.09, 'P'), ('s_1491[c]', 119.12, 'P'), ('s_1561[c]', 117.15, 'P'),
        ('s_1527[c]', 204.23, 'P'), ('s_1533[c]', 181.19, 'P'), ('s_0001[ce]', 180.16, 'C'),
        ('s_0004[ce]', 180.16, 'C'), ('s_0509[c]', 221.21, 'C'), ('s_0773[c]', 180.16, 'C'),
        ('s_1107[c]', 180.16, 'C'), ('s_1520[c]', 342.296, 'C'), ('s_0423[c]', 347.22, 'R'),
        ('s_0526[c]', 323.2, 'R'), ('s_0782[c]', 363.22, 'R'), ('s_1545[c]', 324.18, 'R'),
        ('s_0584[c]', 331.22, 'D'), ('s_0589[c]', 307.2, 'D'), ('s_0615[c]', 345.21, 'D'),
        ('s_0649[c]', 322.21, 'D'), ('s_3714[c]', 852.83, 'N'), ('s_1405[c]', 376.36, 'N'),
        ('s_1467[c]', 96.06, 'N')
    ]

    # Get initial component masses
    X, P, C, R, D, _ = sumBioMass(model, comps)

    protPos = model.rxnNames.index('protein pseudoreaction')
    rnaPos = model.rxnNames.index('RNA pseudoreaction')

    for i, metID in enumerate(otherData['metIDs']):
        abundance = otherData['abundance'][i]

        if metID == 'protein':
            fP = abundance / P
            isAA = ['tRNA' in name for name in model.metNames]
            for j, aa in enumerate(isAA):
                if aa:
                    model.S[j, protPos] *= fP

        elif metID == 'RNA':
            fR = abundance / R
            nucs = [c[0] for c in comps if c[2] == 'R']
            for nuc in nucs:
                modelPos = model.mets.index(nuc)
                model.S[modelPos, rnaPos] *= fR

        else:
            modelPos = model.mets.index(metID)
            comp = next(c for c in comps if c[0] == metID)

            if comp[2] == 'C':
                rxnPos = model.rxnNames.index('carbohydrate pseudoreaction')
                MW = comp[1] - 18
            elif comp[2] == 'D':
                rxnPos = model.rxnNames.index('DNA pseudoreaction')
                MW = comp[1] - 18
            elif comp[2] == 'N':
                rxnPos = model.rxnNames.index('biomass pseudoreaction')
                MW = comp[1]

            model.S[modelPos, rxnPos] = -abundance / MW * 1000

    # Recompute biomass and balance with sugars
    X, _, C, _, _, _ = sumBioMass(model, comps)
    delta = X - 1  # Difference to balance

    # Balance with all carbohydrate components
    mets = [c[0] for c in comps if c[2] == 'C']
    massPre = C
    massPost = massPre - delta
    fC = massPost / massPre

    carbPos = model.rxnNames.index('carbohydrate pseudoreaction')
    for met in mets:
        modelPos = model.mets.index(met)
        model.S[modelPos, carbPos] *= fC

    # Estimate polymerization-associated maintenance energy
    _, P, C, R, D, _ = sumBioMass(model, comps)
    GAMpol = P * 37.7 + C * 12.8 + R * 26.0 + D * 26.0 

    return model, GAMpol
