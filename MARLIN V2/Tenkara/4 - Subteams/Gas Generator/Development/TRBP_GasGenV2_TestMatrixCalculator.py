# TRBP Gas Generator V2.0 Test Matrix Calculator

import numpy as np
import pandas as pd
import yaml
import os

# Inputs

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)  # Go up one level
def_path = os.path.join(parent_dir, 'Inputs', 'GG_hardware_definition.yaml')
print(def_path)
config = yaml.safe_load(open(def_path, 'r'))

inputs = pd.read_csv('TestMatrix_Inputs.csv')



# Calculations

# Outputs
