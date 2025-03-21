import re

def getNewIndex(ids):
    """
    Direct translation of getNewIndex from MATLAB to Python (Benjamín J. Sánchez, 2018-09-04).
    Numbers from each id is extracted, works for both reactions (r_###) and metabolites (s_###)
    """
    # Extract numeric parts from the IDs
    numeric_ids = []
    for id_str in ids:
        # Extract all numeric substrings (assuming IDs follow the "s_###" or "r_###" format)
        numbers = re.findall(r'\d+', id_str)
        if numbers:
            numeric_ids.extend([int(n) for n in numbers])

    # Find the maximum number and increment it
    if numeric_ids:
        new_id = max(numeric_ids) + 1
    else:
        new_id = 1  # Start from 1 if no numeric IDs exist

    return str(new_id)
