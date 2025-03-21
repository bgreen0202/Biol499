import cobra
import pandas as pd

"""
Direct translation of listaGem a mixture of snippets from "reconstructionProtocol" and "addSLIMEreactions"
"""

# Load model (assumes it's in JSON format, replace with actual path)
model = cobra.io.load_json_model("model.json")

# Function to remove existing SLIME reactions
slime_rxns = [rxn.id for rxn in model.reactions if "SLIME rxn" in rxn.name]
model.remove_reactions(slime_rxns)

# Load SLIME template
slime_template = pd.read_csv("SLIMERtemplates.txt", sep="\t")
slime_template.columns = ["metName", "bbID", "bbMW", "comps"] + [f"chain{i}" for i in range(len(slime_template.columns) - 4)]
slime_template = slime_template.fillna("")

# Define known fatty acid data
dat = {
    "chain": ["16:0", "16:1", "18:0", "18:1", "18:2", "18:3"],
    "name": ["palmitate", "palmitoleate", "stearate", "oleate", "linoleate", "linolenate"],
    "MW": [256.42, 254.4, 284.48, 282.46, 280.44, 278.43],
    "ID": ["s_3740", "s_3741", "s_3742", "s_3743", "m_0102", "m_0103"]
}

# Iterate through the SLIME template
for _, row in slime_template.iterrows():
    chains = [row[f"chain{i}"] for i in range(len(dat["chain"])) if row[f"chain{i}"]]
    for chain_set in chains:
        chain_list = chain_set.split(',')
        lipid = row["metName"]
        
        if lipid.startswith("fatty acid"):
            idx = dat["chain"].index(chain_list[0])
            lipid = dat["name"][idx]
        else:
            for i, chain in enumerate(chain_list, 1):
                lipid = lipid.replace(f"CHAIN{i}", chain)
        
        unique_chains = list(set(chain_list))
        coeffs = [dat["MW"][dat["chain"].index(c)] for c in unique_chains]
        totalMW = row["bbMW"] + sum(coeffs)
        
        substrates = [f"{lipid}[{row['comps']}]"]
        products = [row["bbID"]] + [dat["ID"][dat["chain"].index(c)] for c in unique_chains]
        coeffs = [-1] + [mw / 1000 for mw in [totalMW] + coeffs]
        
        reaction = cobra.Reaction(f"t_{len(model.reactions) + 1}")
        reaction.name = f"{lipid} [{row['comps']}] SLIME rxn"
        reaction.lower_bound = 0
        reaction.upper_bound = 1000
        
        for met, coef in zip(substrates + products, coeffs):
            if met not in model.metabolites:
                model.add_metabolites([cobra.Metabolite(met)])
            reaction.add_metabolites({model.metabolites.get_by_id(met): coef})
        
        model.add_reactions([reaction])

# Save the modified model
cobra.io.save_json_model(model, "model_r3.json")

# Export model statistics
print(f"Number of genes: {len(model.genes)}, Reactions: {len(model.reactions)}, Metabolites: {len(model.metabolites)}")
