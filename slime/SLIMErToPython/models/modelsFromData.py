from SLIMErToPython.models import SLIMEr
from SLIMErToPython.models import adjustModel, changeOtherComp, scaleAbundancesInModel


def modelsFromData(model, data, scaling):
    """
    Direct translation of modelsFromData from MATLAB to Python (Benjamín J. Sánchez, 2018-09-04).
    It takes the original model and experimental data, then outputs new models with correction.

    Parameters:
        data (dict): Experimental data, including lipid composition and other components.
        scaling (float): Scaling factor for adjusting abundances.

    Returns:
        tuple: (model_corrComp, model_SLIMEr, k, GAMpol)
            model_corrComp: Model with corrected lipid composition.
            model_SLIMEr: Model with lipid and chain length constraints.
            k: Scaling factor applied to abundances.
            GAMpol: Polymerization-associated growth-associated maintenance energy.


    ***note that this is a direct translation rn and we'll need to make adjustments for our model
    """

    # Model with lipid composition corrected
    model_corrComp = SLIMEr(model, data, constrain_chain_lengths=False)

    # Model with both lipid composition and chain length constrained
    model_SLIMEr = SLIMEr(model, data, constrain_chain_lengths=True)

    # Correct the biomass composition with data for both models
    model_corrComp, _ = changeOtherComp(model_corrComp, data)
    model_SLIMEr, _ = changeOtherComp(model_SLIMEr, data)

    # Make abundances of backbones & chains consistent
    model_SLIMEr, k = scaleAbundancesInModel(model_SLIMEr, data, scaling)

    # Adjust model_corrComp using the same scaling factor k
    model_corrComp = adjustModel(model_corrComp, k, constrain_chain_lengths=False, scaling=scaling)

    # Correct biomass composition again to account for updated lipid content
    model_corrComp, _ = changeOtherComp(model_corrComp, data)
    model_SLIMEr, GAMpol = changeOtherComp(model_SLIMEr, data)

    return model_corrComp, model_SLIMEr, k, GAMpol
