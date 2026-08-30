"""
This file contains functions pertaining to compressible flow.

Index:
 - mach_to_area

Future implementations:
 - area_to_mach
 - total_properties
 - normal_shock

Author: Louis DeSano
Created: 08/23/2026
Last Modified: 08/23/2026
"""

import numpy as np

def mach_to_area(chamberDiameter, mach, gam=1.4, supersonic=True):
    """
    Description: Determine area from mach.

    Methods: Isentropic Flow Area-Mach relation

    Arguments:
        chamberDiameter (float) Plenum Diameter [length^2] (unit impartial)
        mach (float) Local mach number
        gam (float) Specific heat ratio of working fluid

    Returns: 
        area (float) Local area [length^2]
        
    """
    # determine area from diameter
    chamberArea = (np.pi/4) * chamberDiameter**2
    
    # isentropic flow area ratio
    area_ratio = (1/mach) * ((2/(gam+1)) * (1 + ((gam-1)/2)*(mach**2)))**(gam+1)/(2*(gam-1))
    
    # 
    
    """if supersonic:

    else:"""


    return area




