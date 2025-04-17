import re

def getBackboneName(spec_name):
    """
    Direct translation of the MATLAB function getBackboneName into Python using COBRApy.
    Works with a single metabolite ID and parses out the backbone name.
    Handles parsing out the compartment from name.

    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """

    back_name = ""

     # Group 1: keep generic name if conditions match
    group1 = [
        "phosphatidylserine",
        "phosphatidylethanolamine",
        "phosphatidate",
        "phosphatidylglycerol",
        "cardiolipin"
    ]

    for name in group1:
        if spec_name.lower().startswith(name.lower()) and spec_name.lower() != name.lower() and "phosphate" not in spec_name.lower():
            back_name = name 
            break

    # Group 2: replace specific name with generic backbone name
    if not back_name:
        group2 = {
            "palmitate": "fatty acid",
            "palmitoleate": "fatty acid",
            "stearate": "fatty acid",
            "oleate": "fatty acid",
        }

        if spec_name in group2:
            back_name = group2[spec_name]

    return back_name