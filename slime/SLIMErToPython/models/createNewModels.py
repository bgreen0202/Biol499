import os
import numpy as np
from cobra import Model, io

from slime.SLIMErToPython.data import convertEsjingData, readEsjingData, readLahtveeData
from slime.SLIMErToPython.models import modelsFromData

"""
Direct translation of createNewModels from MATLAB to Python (Benjamín J. Sánchez, 2018-09-04).
Loads the original model, creates various new models, reads the data for the associated model.

"""

# WE NEED TO MAKE THIS BASED ON OUR DATA BUT I DONT UNDERSTAND WHAT THEY ARE DOING EXACTLY
# ONE SET OF DATA SEEMS TO BE FOR SLIME, ONE SET SEEMS TO BE FOR STUDIES THAT VALIDATE ??
# WOULD WE DO THIS BECAUSE CATHERINES PAPER HAS 4 CONDITIONS THEN TECHNICALLY WE HAVE 6 WITH OUR VALIDATIONS ??

model = io.read_sbml_model('yeast780.xml')  # or 'yeast780.json', depending on format

# Create 2 models for each of the 10 conditions in the stress dataset
# One model is a "corrected" (corrComp) model and the other is the model adjusted for slime (SLIMEr)
model_corrComp = [None] * 10
model_SLIMEr = [None] * 10
GAMpol = np.zeros(10)
k = np.zeros(10)

for i in range(10):
    data = readLahtveeData(i + 1)  
    model_corrComp[i], model_SLIMEr[i], k[i], GAMpol[i] = modelsFromData(model, data, 'backbones')

# Create 2 models for each of the validation studies
# One model is a "corrected" (corrComp) model and the other is the model adjusted for slime (SLIMEr)
model_corrComp_val = [None] * 8
model_SLIMEr_val = [None] * 8

for i in range(8):
    data = readEsjingData(i + 1)  
    data = convertEsjingData(data, model, True)
    model_corrComp_val[i], model_SLIMEr_val[i], _, _ = modelsFromData(model, data, 'tails')


