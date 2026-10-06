import numpy as np
from pyfluids import Fluid, FluidsList, Input

g_f = 1.31  # ch4 specific heat ratio
g_ox = 1.40 # gox specific heat ratio  

pc = 192.4    #psia

p_crit_ox = pc * ((g_ox + 1) / 2) ** (g_ox/(g_ox -1))

p_crit_f = pc * ((g_f + 1) / 2) ** (g_f/(g_f -1))

crit_stiff_ox = ((p_crit_ox - pc) / pc) * 100

crit_stiff_f = ((p_crit_f - pc) / pc) * 100

print(f"The minimum feed pressure for oxygen is {p_crit_ox : .3f} psia")
print(f"The minimum stiffness for oxygen is {crit_stiff_ox : .2f}%")
print(f"The minimum feed pressure for fuel is {p_crit_f : .3f} psia")
print(f"The minimum stiffness for fuel is {crit_stiff_f : .2f}%")