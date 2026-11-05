from pint import Quantity as Q_



class hardware:

    LOx_Cda = Q_(0.00000175, 'm^2')
    Fuel_Cda = Q_(0.00000307, 'm^2')
    Film_Cda = Q_(0.00000791, 'm^2')   
    Insert_Cda = {
        '0%':  Q_(0.0, 'm^2'),
        '5%':  Q_(0.0000002678, 'm^2'),
        '10%': Q_(0.0000006026, 'm^2'),
        '15%': Q_(0.0000009719, 'm^2'),
        '20%': Q_(0.000001366, 'm^2'),
        '25%': Q_(0.000001728, 'm^2'),
        '30%': Q_(0.000002285, 'm^2'),
        '35%': Q_(0.000002920, 'm^2'),
        '40%': Q_(0.000003888, 'm^2')
    }
            
    Stator_Area = Q_(84, 'mm^2') # Assumed to be the same as the injector area #ORIGINALY WAS 0.000077 m^2
                                # Last updated: Sept 27th, 2026 by Miles Cohen
                                #TENTATIVE/RANDOM

    stator_throat_num = Q_(4, 'dimensionless') # Number of nozzles in the stator
                                # Last updated: May 25th, 2026 by Miles Cohen

    stator_eff_combustion_volume = Q_(0, 'm^3') # This is the volume of the HGM that is effectivly used for combustion.
                                    # Last updated: May 25th, 2026 by Miles Cohen

    CSTAR_eff = Q_(0.7, 'dimensionless') # Last updated October 9th by Miles Cohen

    # Combustion Chamber Hardware Definition
    combustion_chamber_volume = Q_(0.0001393, 'm^3') # This is the total volume of the combustion chamber.
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
    CH4_Cda = Q_(0.000000438, 'm^2') # Assumed Cd of 0.9 
                                # Last updated: Nov 2nd, 2026 by Arvin Peng

    GOx_Cda = Q_(0.000000438, 'm^2') # Assumed Cd of 0.9
                                # Last updated: Nov 2nd, 2026 by Arvin Peng

    torch_cstar_eff = Q_(0.9, 'dimensionless') # Assumed to be 1.0
                                # Last updated: May 25th, 2026 by Miles Cohen

    torch_combustion_volume = Q_(0, 'm^3') # This is the volume of the torch that is used for combustion.
                                # Last updated: May 25th, 2026 by Miles Cohen

    torch_throat_area = Q_(0.000007917, 'm^2') # This is the area of the torch throat.
                                # Last updated: Sept 28th, 2026 by Miles Cohen

    torch_eff_throat_area = Q_(0.000007917, 'm^2') # This is the effective area of the torch throat, accounting for any losses.
                                            
                                # Efficiency is assumed to be 1.0
                                # Last updated: Sept 28th, 2026 by Miles Cohen

    torch_fuel_type = "CH4"         # NASA CEA input for torch propellant
    CH4_gamma = Q_(1.304, 'dimensionless') # Methane specific heat ratio
    CH4_R = Q_(518.28, 'J/(kg*K)') # Methane specific gas constant (R/M)

    torch_ox_type = "O2"            # NASA CEA input for torch oxidizer
    GOx_gamma = Q_(1.395, 'dimensionless') # Gaseous oxygen specific heat ratio
    GOx_R = Q_(259.84, 'J/(kg*K)') # Oxygen specific gas constant (R/M)
        
