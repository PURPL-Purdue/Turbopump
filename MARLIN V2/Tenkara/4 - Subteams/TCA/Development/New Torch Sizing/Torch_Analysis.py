# Torch_Analysis
# Torch_Analysis
from New_TCA_Igniter_Sizing import choked_backpressure

## Define torch parameters
At = 1     # Throat area
S_f = 1     # Fuel injector size
S_o = 1     # Oxidizer injector size
k = 1       #Specific heat ratio
p_c = 1     # Torch Chamber Pressure
stiff = 0.30 # Stiffness of the torch

line_pressure = choked_backpressure(k, stiff, p_c)


