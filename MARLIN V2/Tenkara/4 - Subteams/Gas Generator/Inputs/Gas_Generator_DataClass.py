# Gas Generator Hardware Definition
# Last Updated: Nov 10, 2026

from pint import Quantity as Q_

class hardware:

    # Injector Parameters
    LOx_Cda = Q_(0.00000190, 'm^2')     # Assumed Cd of 0.65 w/ #68 drill bit
                                        # Last updated: Nov 5th, 2026 by Miles Cohen

    Fuel_Cda = Q_(0.00000435, 'm^2')    # Assumed Cd of 0.65 w/ 3/64 drill bit
                                        # Last updated: Nov 5th, 2026 by Miles Cohen

    Film_Cda = Q_(0.00001028, 'm^2')    # Assumed Cd of 0.65
                                        # Last updated: Oct 28th, 2026 by Arvin Peng

    # Film Cooling Inserts
    Insert_Cda = {                      # (m^2) Assumed Cd = 0.65
        '0%':  Q_(0.0, 'm^2'),
        '5%':  Q_(0.000000319, 'm^2'),
        '10%': Q_(0.000000711, 'm^2'),
        '15%': Q_(0.00000115, 'm^2'),
        '20%': Q_(0.00000162, 'm^2'),
        '25%': Q_(0.00000204, 'm^2'),
        '30%': Q_(0.00000270, 'm^2'),
        '35%': Q_(0.00000345, 'm^2'),
        '40%': Q_(0.00000459, 'm^2')
    }                                   # Last updated: Nov 5th, 2026 by Arvin Peng

    # Stator Parameters
    Stator_Area = Q_(98, 'mm^2')    # Assumed to be the same as the throat area
                                    # Last updated: Sept 27th, 2026 by Miles Cohen

    Stator_throat_num = Q_(4, 'dimensionless')  # Number of nozzles in the stator
                                                # Last updated: May 25th, 2026 by Miles Cohen

    Stator_eff_combustion_volume = Q_(0, 'm^3') # Volume of the HGM that is effectivly used for combustion.
                                                # Last updated: May 25th, 2026 by Miles Cohen

    Cstar_eff = Q_(0.7, 'dimensionless')    # Estimated C* efficiency
                                            # Last updated October 9th by Miles Cohen

    # Combustion Chamber Hardware Definition
    combustion_chamber_volume = Q_(0.0001393, 'm^3')    # Total volume of the combustion chamber.
                                                        # Last updated: Sept 28th, 2026 by Miles Cohen

    # Propellant Definitions
    GG_fuel_type = "Isopropanol"
    IPA_rho = Q_(785, 'kg/m^3')
    
    cp_IPA_l = Q_(0.1933, 'kJ/(mol*K)')
    h_vap_IPA = Q_(45, 'kJ/mol')
    Hvap = Q_(45, 'kJ/mol')
    cp_fuel_g_m = Q_(0.0001331, 'kJ/(mol*K^2)') # Linear model slope for gaseous IPA heat capacity
    cp_fuel_g_b = Q_(0.06588, 'kJ/(mol*K)') # Heat capacity intercept

    GG_ox_type = "LO2"
    LOx_rho = Q_(1141, 'kg/m^3')

    # Torch Hardware Definitions
    CH4_Cda = Q_(0.000000438, 'm^2')    # Assumed Cd of 0.9 
                                        # Last updated: Nov 2nd, 2026 by Arvin Peng

    GOx_Cda = Q_(0.000000438, 'm^2')    # Assumed Cd of 0.9
                                        # Last updated: Nov 2nd, 2026 by Arvin Peng

    torch_cstar_eff = Q_(0.9, 'dimensionless')  # Estimated torch C* efficiency
                                                # Last updated: May 25th, 2026 by Miles Cohen

    torch_combustion_volume = Q_(0, 'm^3')      # This is the volume of the torch that is used for combustion.
                                                # Last updated: May 25th, 2026 by Miles Cohen

    torch_throat_area = Q_(0.000007917, 'm^2')  # This is the area of the torch throat.
                                                # Last updated: Sept 28th, 2026 by Miles Cohen

    torch_eff_throat_area = Q_(0.000007917, 'm^2')  # This is the effective area of the torch throat, accounting for any losses.
                                                    # Last updated: Sept 28th, 2026 by Miles Cohen

    # Torch Propellant Definitions
    torch_fuel_type = "CH4"         # NASA CEA input for torch propellant
    CH4_gamma = Q_(1.304, 'dimensionless') # Methane specific heat ratio
    CH4_R = Q_(518.28, 'J/(kg*K)') # Methane specific gas constant (R/M)

    torch_ox_type = "O2"            # NASA CEA input for torch oxidizer
    GOx_gamma = Q_(1.395, 'dimensionless') # Gaseous oxygen specific heat ratio
    GOx_R = Q_(259.84, 'J/(kg*K)') # Oxygen specific gas constant (R/M)
        