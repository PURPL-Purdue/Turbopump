##################################
#Importing Necessary Libraries
##################################
import numpy as np
import csv
import yaml
from rocketcea.cea_obj import CEA_Obj

np.set_printoptions(legacy='1.25')    #Fix for some numpy float printing isssues that happened

##################################
#Define NASA Chemical Equilibrium with Applications (CEA) Object
##################################

#Define CEA Object as C
C = CEA_Obj( oxName='LOX', fuelName='RP1')

##################################
#Define Global Conversion Factors
##################################

#psia to pascals conversion
psi_to_pa = 6894.76
#feet to meters conversion
ft_to_m = 0.3048
#m to cm
mcm = 100
#Rankine to Kelvin conversion factor
rankineToKelvin = 5.0 / 9.0
#Rankine to Fahrenheit conversion
rankineToF = -459.67
#lbf to N conversion, N = kg-m/s^2
lbf_N = 4.4482216
#universal gas constant (J / kmol-K)
Ru_m = 8314.462618
#m to in
m_in = 39.3701
#gravity in m/s^2
g = 9.81
#kg to lbm conversion
kg_to_lbm = 2.20462
# Specific Heat Capacity SI Conversion (to W / m-K)
# Default assumes input from CEA is in kcal/kg-K or cal/g-K basis that requires
# multiplying by ~4184 to get J/(kg-K). This can be overridden in `TCA_params.yaml`
# by setting `cp_conversion` to the desired factor.
c_siconv = 4186.8
#Thermal Conductivity SI Conversion (to J / kg-K)
tc_siconv = 418.4
#Dynamic Viscosity SI Conversion  (to Pa-s)
vis_siconv = 0.0001

#MAIN FUNCTION

pe = 14.7   #Exit pressure

with open(r'C:\Users\RAFA\Desktop\TCA_params.yaml') as file:
        tca_params = yaml.safe_load(file)

# Allow overriding conversion factors from the yaml config
c_siconv = tca_params.get('cp_conversion', c_siconv)
tc_siconv = tca_params.get('tc_conversion', tc_siconv)
vis_siconv = tca_params.get('visc_conversion', vis_siconv)

print(f"Using conversion factors: cp_conversion={c_siconv}, tc_conversion={tc_siconv}, visc_conversion={vis_siconv}")

#Contraction Ratio
con_r = tca_params['tca_contraction_ratio']  #contraction ratio chamber area/throat area

mrLo = float(input('What is your lower bound for O/F ratio?: '))
mrHi = float(input('What is your upper bound for O/F ratio?: '))
mrInt = float(input('What interval would you like between your O/F ratios?: '))
pc = float(input('What chamber pressure are you operating at (psia)?: '))

MRs = np.arange(mrLo, (mrHi+ mrInt), mrInt)

idx = 0

# Define the number of rows and columns for the matrix
rows = len(MRs)
cols = 19

data = []

row1 = ['O/F', 'Chamber Heat Capacity (W/m-K)', 'Throat Heat Capacity (W/m-K)', 'Exit Heat Capacity (W/m-K)', 'Chamber Viscosity (Pa-s)', 
        'Throat Viscosity (Pa-s)', 'Exit Viscosity (Pa-s)', 'Chamber Thermal Conductivity (J/kg-K)', 'Throat Thermal Conductivity (J/kg-K)', 
        'Exit Thermal Conductivity(J/kg-K)', 'Chamber Prandtl Number', 'Throat Prandtl Number', 'Exit Prandtl Number', 'Chamber Mach', 'Throat Mach', 'Exit Mach', 
        'Chamber Molecular Weight (kg/kmol)', 'Throat Molecular Weight (kg/kmol)', 'Exit Molecular Weight (kg/kmol)', 
        'Chamber Temperature (F)', 'Throat Temperature (F)', 'Exit Temperature (F)']

data.append(row1)

for ROW in range(rows):
    row = []
    for mr in MRs:
        Eps = C.get_eps_at_PcOvPe(Pc = pc, MR = mr, PcOvPe = (pc / pe), frozen = 0)

        heat_cap_c, viscosity_c, therm_con_c, p_num_c = C.get_Chamber_Transport(Pc = pc, MR = mr, eps = Eps, frozen=0)
        CM = C.get_Chamber_MachNumber(Pc = pc, MR = mr, fac_CR = con_r)
        MW_C = C.get_Chamber_MolWt_gamma(Pc = pc, MR = mr, eps = Eps)

        hcc = heat_cap_c * c_siconv
        vc = viscosity_c * vis_siconv
        tcc = therm_con_c * tc_siconv
        # Print raw and converted cp once for debugging/verification
        if idx == 0:
            try:
                print(f"Raw CEA chamber heat capacity (raw) = {heat_cap_c}")
                print(f"Converted chamber cp (J/kg-K) = {hcc}")
            except Exception:
                pass
            idx = 1
    
        heat_cap_t, viscosity_t, therm_con_t, p_num_t = C.get_Throat_Transport(Pc = pc, MR = mr, eps = Eps, frozen=0)
        TM = 1
        MW_T = C.get_Throat_MolWt_gamma(Pc = pc, MR = mr, eps = Eps)

        hct = heat_cap_t * c_siconv
        vt = viscosity_t * vis_siconv
        tct = therm_con_t * tc_siconv        
    
        heat_cap_e, viscosity_e, therm_con_e, p_num_e = C.get_Exit_Transport(Pc = pc, MR = mr, eps = Eps, frozen=0)
        EM = C.get_MachNumber(Pc = pc, MR = mr, eps = Eps, frozen=0)
        MW_E = C.get_exit_MolWt_gamma(Pc = pc, MR = mr, eps = Eps, frozen=0)

        hce = heat_cap_e * c_siconv
        ve = viscosity_e * vis_siconv
        tce = therm_con_e * tc_siconv
        
        # Get temperatures at chamber, throat, and exit (in Rankine from CEA)
        Tc_R, Tt_R, Te_R = C.get_Temperatures(Pc = pc, MR = mr, eps = Eps, frozen=0, frozenAtThroat=0)
        # Convert from Rankine to Fahrenheit
        Tc_F = Tc_R + rankineToF
        Tt_F = Tt_R + rankineToF
        Te_F = Te_R + rankineToF
        
        # Convert temperatures from Rankine to Kelvin for polynomial fitting
        Tc_K = Tc_R * rankineToKelvin
        Tt_K = Tt_R * rankineToKelvin
        Te_K = Te_R * rankineToKelvin
        

        row.append([mr, hcc, hct, hce, vc, vt, ve, tcc, tct, tce, CM, TM, EM, MW_C[0], MW_T[0], MW_E[0], Tc_F, Tt_F, Te_F])

        data.append(row)
    
        # Exit the loop after first iteration since we only need data once
        break

# ==================== COLLECT DATA FOR POLYNOMIAL FITTING ====================
# After collecting all data in the CSV, now we collect temperature and property data
# Lists to hold all temperature and property values for polynomial fitting
temps_K = []
cp_values = []
thermal_k_values = []
viscosity_values = []

# Go through data again to extract temperatures and properties for polynomial fitting
for row_data in data[1:]:  # Skip header row
        if isinstance(row_data[0], list):
                # Handle nested list structure
                for item in row_data:
                        # Each item is [mr, hcc, hct, hce, vc, vt, ve, tcc, tct, tce, CM, TM, EM, MW_C[0], MW_T[0], MW_E[0], Tc_F, Tt_F, Te_F]
                        # Tc_F, Tt_F, Te_F are at indices 16, 17, 18
                        Tc_F = item[16]
                        Tt_F = item[17]
                        Te_F = item[18]
                        # Convert back to Kelvin
                        Tc_K = (Tc_F - rankineToF) * rankineToKelvin
                        Tt_K = (Tt_F - rankineToF) * rankineToKelvin
                        Te_K = (Te_F - rankineToF) * rankineToKelvin
                        # Extract properties (indices: hcc=1, hct=2, hce=3, vc=4, vt=5, ve=6, tcc=7, tct=8, tce=9)
                        temps_K.extend([Tc_K, Tt_K, Te_K])
                        cp_values.extend([item[1], item[2], item[3]])
                        thermal_k_values.extend([item[7], item[8], item[9]])
                        viscosity_values.extend([item[4], item[5], item[6]])


output_file_path = r"dimensions.csv"

with open(output_file_path, 'w', newline='') as csvfile:
  csv_writer = csv.writer(csvfile)
  csv_writer.writerows(data)

# ==================== POLYNOMIAL FITTING ====================
# Fit polynomials to the collected data
# Using cubic polynomials (degree=3) for good accuracy

print("\n" + "="*60)
print("POLYNOMIAL INTERPOLATION RESULTS")
print("Independent Variable: Temperature (Kelvin)")
print("="*60)

# Convert lists to numpy arrays for polyfit
temps_K_arr = np.array(temps_K)
cp_arr = np.array(cp_values)
thermal_k_arr = np.array(thermal_k_values)
viscosity_arr = np.array(viscosity_values)

# Fit cubic polynomials (degree 3)
poly_degree = 3

cp_coeffs = np.polyfit(temps_K_arr, cp_arr, poly_degree)
thermal_k_coeffs = np.polyfit(temps_K_arr, thermal_k_arr, poly_degree)
viscosity_coeffs = np.polyfit(temps_K_arr, viscosity_arr, poly_degree)

# Print results
print(f"\nPolynomial Degree: {poly_degree}")
print("\n1. Specific Heat Capacity (cp) - W/(m-K) vs Temperature (K)")
print(f"   Coefficients (highest to lowest degree): {cp_coeffs}")

print("\n2. Thermal Conductivity - J/(kg-K) vs Temperature (K)")
print(f"   Coefficients (highest to lowest degree): {thermal_k_coeffs}")

print("\n3. Dynamic Viscosity - Pa-s vs Temperature (K)")
print(f"   Coefficients (highest to lowest degree): {viscosity_coeffs}")

# Write coefficients to a separate output file
poly_output_file = r"polynomial_coefficients.txt"
with open(poly_output_file, 'w') as f:
        f.write("POLYNOMIAL INTERPOLATION COEFFICIENTS\n")
        f.write("Independent Variable: Temperature (Kelvin)\n")
        f.write(f"Polynomial Degree: {poly_degree}\n")
        f.write(f"Chamber Pressure: {pc} psia\n")
        f.write(f"O/F Ratio Range: {mrLo} to {mrHi}\n")
        f.write("="*60 + "\n\n")
    
        f.write("1. Specific Heat Capacity (cp) - W/(m-K) vs Temperature (K)\n")
        f.write(f"   y = {cp_coeffs[0]:.6e}*T^3 + {cp_coeffs[1]:.6e}*T^2 + {cp_coeffs[2]:.6e}*T + {cp_coeffs[3]:.6e}\n")
        f.write(f"   Coefficients: {cp_coeffs}\n\n")
    
        f.write("2. Thermal Conductivity - J/(kg-K) vs Temperature (K)\n")
        f.write(f"   y = {thermal_k_coeffs[0]:.6e}*T^3 + {thermal_k_coeffs[1]:.6e}*T^2 + {thermal_k_coeffs[2]:.6e}*T + {thermal_k_coeffs[3]:.6e}\n")
        f.write(f"   Coefficients: {thermal_k_coeffs}\n\n")
    
        f.write("3. Dynamic Viscosity - Pa-s vs Temperature (K)\n")
        f.write(f"   y = {viscosity_coeffs[0]:.6e}*T^3 + {viscosity_coeffs[1]:.6e}*T^2 + {viscosity_coeffs[2]:.6e}*T + {viscosity_coeffs[3]:.6e}\n")
        f.write(f"   Coefficients: {viscosity_coeffs}\n")

print(f"\n✓ Polynomial coefficients saved to: {poly_output_file}")
print("="*60)

