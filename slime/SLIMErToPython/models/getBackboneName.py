import re

def getBackboneName(spec_name):
    """
    Direct translation of the MATLAB function getBackboneName into Python using COBRApy.
    Works with a single metabolite ID and parses out the backbone name.
    Handles parsing out the compartment from name.
    """

    back_name = ""

    # Extract compartment (compName) and core species name
    match = re.search(r'\[(.*?)\]', spec_name)
    if not match:
        return back_name  # If no compartment, return empty

    comp_name = match.group(1)
    spec_name = spec_name[:match.start()].strip()

    # Define the groups GROUPS WE WILL DECIDE
    # ALSO MAY NOT NECESSARILY NEED THE COMPARTMENT
    group1 = [
        ("phosphatidylserine", "endoplasmic reticulum membrane"),
        ("phosphatidylethanolamine", "endoplasmic reticulum membrane"),
        ("phosphatidate", "endoplasmic reticulum membrane"),
        ("phosphatidylglycerol", "mitochondrial membrane"),
        ("cardiolipin", "mitochondrial membrane")
    ]

    # Group 1: keep generic name if conditions match
    for name, location in group1:
        if spec_name.startswith(name) and spec_name != name and comp_name == location and "phosphate" not in spec_name:
            back_name = f"{name} [{comp_name}]"
            return back_name

    # Group 2: replace specific name with generic backbone name
    # GROUPS WE WILL DECIDE
    # ALSO MAY NOT NECESSARILY NEED THE COMPARTMENT
    group2 = [
        ("palmitate", "fatty acid", "cytoplasm"),
        ("palmitoleate", "fatty acid", "cytoplasm"),
        ("stearate", "fatty acid", "cytoplasm"),
        ("oleate", "fatty acid", "cytoplasm"),
        ("ergosterol", "ergosterol", "cytoplasm"),
        ("ergosteryl palmitoleate", "ergosterol ester", "endoplasmic reticulum membrane"),
        ("ergosteryl oleate", "ergosterol ester", "endoplasmic reticulum membrane"),
        ("sphinganine", "long-chain base", "endoplasmic reticulum"),
        ("phytosphingosine", "long-chain base", "endoplasmic reticulum"),
        ("sphinganine 1-phosphate", "long-chain base phosphate", "endoplasmic reticulum"),
        ("phytosphingosine 1-phosphate", "long-chain base phosphate", "endoplasmic reticulum")
    ]

    for name, generic, location in group2:
        if spec_name == name and comp_name == location:
            back_name = f"{generic} [{comp_name}]"
            return back_name

    return back_name
