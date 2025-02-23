import json

# List of JSON file paths
json_files = [
    "Anteisocardiolipin to Even lipids.json",
    "Isoheptadecanoylcardiolipin to odd lipids.json",
    "OtherNamedFAs.json",
    "Phosphatidylserine_didodecanoyl and Stearoylcardiolipin.json",
    "FA New w_Biomass.json"
]

unique_rxn_values = set()

for file in json_files:
    with open(file, "r") as f:
        data = json.load(f)  # Load JSON content
        
        # Loop through the list (since JSON is a list of dictionaries)
        for entry in data:
            if isinstance(entry, dict) and "reactions" in entry:
                for reaction in entry["reactions"].values():
                    # Extract and check "name"
                    if "name" in reaction and isinstance(reaction["name"], str) and reaction["name"].startswith("rxn"):
                        unique_rxn_values.add(reaction["name"])
                    
                    # Extract and check "bigg_id"
                    if "bigg_id" in reaction and isinstance(reaction["bigg_id"], str) and reaction["bigg_id"].startswith("rxn"):
                        unique_rxn_values.add(reaction["bigg_id"])

# Convert set to sorted list
result = sorted(unique_rxn_values)

# Write the output to a JSON file
with open("output_reactions.json", "w") as output_file:
    json.dump(result, output_file, indent=4)

print("Results have been written to 'output_reactions.json'.")
