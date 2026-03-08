import math

def calc_radius_turbine(height_blade, radius_hub):
    return 0.5 * height_blade + radius_hub

def calc_mass_blade(Length_blade, width_blade, height_bmin, rho_blade):
    return Length_blade * width_blade * height_bmin * rho_blade

def calc_Force_centrifugal(mass_blade, w, radius_turbine):
    return mass_blade * (w ** 2) * radius_turbine

def calc_stress_centrifugal(radius_turbine, rho_blade, height_blade, w):
    return radius_turbine * rho_blade * height_blade * (w ** 2)

def calc_Force_tangential(m_dot, V_1, beta_1, V_2, beta_2):
    return m_dot * (V_1 * math.cos(beta_1) + V_2 * math.cos(beta_2))

def calc_Force_axial(m_dot, C_1, alpha_1, V_2, beta_2):
    return m_dot * (C_1 * math.sin(alpha_1) + V_2 * math.sin(beta_2))

def calc_torque_blade(Force_tangential, height_blade):
    return 0.5 * Force_tangential * height_blade

def calc_torque_turbine(Force_tangential, radius_turbine, Z_blade):
    return Force_tangential * radius_turbine * Z_blade

def calc_P(torque_turbine, w):
    return torque_turbine * w

def calc_Force_gas(Force_tangential, Force_axial):
    return math.sqrt(Force_tangential ** 2 + Force_axial ** 2)

def calc_Moment_Bending(height_blade, Z_blade, Force_gas):
    return (height_blade / (2 * Z_blade)) * Force_gas

def calc_I(Length_blade, height_bmin):
    return (1 / 12) * (Length_blade ** 3) * height_bmin

def calc_stress_gas(height_bmin, Force_gas, width_b, I):
    return (0.5 * height_bmin * Force_gas * width_b) / I