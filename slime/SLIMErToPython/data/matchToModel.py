import numpy as np
from itertools import permutations

def match_to_model(model, met_name):
    """
    Checks metabolite name against known mappings, processes tail group variations and returns matching positions in model.
    """

    codes = {
        'TAG': ('triglyceride', 'endoplasmic reticulum membrane'),
        'DAG': ('diglyceride', 'endoplasmic reticulum membrane'),
        'PA': ('phosphatidate', 'endoplasmic reticulum membrane'),
        'PS': ('phosphatidyl-L-serine', 'endoplasmic reticulum membrane'),
        'PE': ('phosphatidylethanolamine', 'endoplasmic reticulum membrane'),
        'PC': ('phosphatidylcholine', 'endoplasmic reticulum membrane'),
        'LPI': ('sn-2-acyl-1-lysophosphatidylinositol', 'endoplasmic reticulum membrane'),
        'PI': ('1-phosphatidyl-1D-myo-inositol', 'cytoplasm'),
        'PG': ('phosphatidylglycerol', 'mitochondrial membrane'),
        'CL': ('cardiolipin', 'mitochondrial membrane'),
        'LCB 18:0;2': ('sphinganine', 'endoplasmic reticulum'),
        'LCB 18:0;3': ('phytosphingosine', 'endoplasmic reticulum'),
        'LCBP 18:0;2': ('sphinganine 1-phosphate', 'endoplasmic reticulum'),
        'LCBP 18:0;3': ('phytosphingosine 1-phosphate', 'endoplasmic reticulum'),
        'Cer': ('ceramide', 'Golgi'),
        'IPC': ('inositol-P-ceramide', 'Golgi'),
        'MIPC': ('mannosylinositol phosphorylceramide', 'Golgi'),
        'M(IP)2C': ('inositol phosphomannosylinositol phosphoceramide', 'Golgi'),
        'Ergosterol': ('ergosterol', 'cytoplasm')
    }
    
    pos = np.zeros(len(model.mets), dtype=bool)
    
    if 'LCB' in met_name or 'Ergost' in met_name:
        if met_name in codes:
            backbone, compartment = codes[met_name]
            formatted_name = f"{backbone} [{compartment}]"
            pos = np.array([name == formatted_name for name in model.metNames])
    else:
        parts = met_name.split(' ')
        if len(parts) < 2:
            return pos
        
        back_code, tail_code = parts[0], parts[1]
        if back_code in codes:
            backbone, compartment = codes[back_code]
            tails = tail_code.split('-')
            
            if back_code in {'TAG', 'DAG', 'PA', 'PS', 'PE', 'PC', 'LPI', 'PI', 'PG', 'CL'}:
                for tail_combo in permutations(tails):
                    if back_code == 'LPI':
                        tail_str = tail_combo[0]
                    else:
                        tail_str = ', '.join([f"{i+1}-{tail}" for i, tail in enumerate(tail_combo)])
                    
                    formatted_name = f"{backbone} ({tail_str}) [{compartment}]"
                    pos |= np.array([name == formatted_name for name in model.metNames])
            elif back_code in {'Cer', 'IPC', 'MIPC', 'M(IP)2C'}:
                sphingos = [
                    ('1', 'A', '2', '0'),
                    ('2', 'B', '3', '0'),
                    ('3', 'C', '3', '1'),
                    ('4', 'D', '3', '2'),
                    ("2'", "B'", '2', '1')
                ]
                
                short_tail, long_tail = tails[0][-1], tails[1][-1]
                for s in sphingos:
                    if s[2] == short_tail and s[3] == long_tail:
                        formatted_name = f"{backbone}-{s[0]} (C{tails[1][:2]}) [{compartment}]"
                        if formatted_name not in model.metNames:
                            formatted_name = f"{backbone} {s[1]} (C{tails[1][:2]}) [{compartment}]"
                        pos |= np.array([name == formatted_name for name in model.metNames])
                        break
    
    return pos
