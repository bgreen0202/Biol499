import pandas as pd

def readEsjingData(i, filename='fullData_Ejsing2009.csv'):
    """
    Reads lipid data from a csv, selects appropriate abundunces and st. dev. columns and filters out unwanted species.
    
    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """


    lipid_data = pd.read_csv(filename)
    
    # Extract relevant columns
    met_names = lipid_data.iloc[:, 0]
    abundance = lipid_data.iloc[:, 2 * i]
    std_dev = lipid_data.iloc[:, 2 * i + 1]
    
    # Filter out unwanted lipid species
    filter_out = met_names.str.contains('10|12|14|20|22|:2|LPA|LPS|LPE|LPC|LIPC|Ergostadienol', regex=True)
    data = {
        'metNames': met_names[~filter_out].tolist(),
        'abundance': abundance[~filter_out].tolist(),
        'std': std_dev[~filter_out].tolist()
    }
    
    return data
