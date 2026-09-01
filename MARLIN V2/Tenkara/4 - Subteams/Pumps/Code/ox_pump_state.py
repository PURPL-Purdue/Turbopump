import numpy as np
import CoolProp as cp
from CoolProp.CoolProp import PropsSI




def ox_pump_state(mdot, safety_factor, shaft_speed_rpm, headrise, surrogate_flag, inlet_temperature, inlet_pressure, break_hp):
    g = 9.81

    if surrogate_flag == True:
        fluid_temp = float(input("Fluid Temperature (Celsius): "))
        density = float(input("Fluid Density (kg / m3): "))
    else:
        fluid_temp =  inlet_temperature 
        fluid_temp_kelvin = fluid_temp + 273.15 # kelvin
        density = PropsSI('D', 'T', fluid_temp_kelvin, 'P', inlet_pressure, "oxygen") # kg/m3
        vapor_pressure = PropsSI('P', 'T', fluid_temp_kelvin, 'Q', 0, 'oxygen') #Pa


    vol_flow = mdot / density
    NPSHa = (inlet_pressure - vapor_pressure) / (density * g)
    target_NPSHr = NPSHa / safety_factor

    def inducer(vol_flow, target_NPSHr, safety_factor, shaft_speed_rpm, headrise):

        #------------------------Geometric Contraints (Get rid of)---------------------------
        tip_radius = 0.0553
        tip_diameter = tip_radius * 2
        hub_radius = 0.0166
        hub_diameter = hub_radius * 2
        blade_number = 3


        tip_diameter = 0.1105
        hub_diameter = 0.0332
        tip_radius = tip_diameter / 2
        hub_radius = hub_diameter / 2
        #----------------------------------------------------------------------------------------
        g = 9.81


        angular_velocity = np.pi * 2 * shaft_speed_rpm / 60
        tip_speed = tip_radius * angular_velocity

        target_NPSHr = NPSHa / safety_factor
        suction_specific_speed = shaft_speed_rpm * np.sqrt(vol_flow) / (target_NPSHr**0.75)
        suction_specific_angular_speed = angular_velocity * np.sqrt(vol_flow) / ((target_NPSHr * g)**0.75)

        cavitation_number = target_NPSHr / ((tip_speed ** 2) / (2 * g))



        flow_coeff_optimal = np.sqrt(cavitation_number / ((2 * cavitation_number) + 3))


        meridional_velocity = flow_coeff_optimal * tip_speed
        hub_tangential_speed = angular_velocity * hub_radius
        tip_flow_angle = np.arctan(np.deg2rad(meridional_velocity / tip_speed))
        hub_flow_angle = np.arctan(np.deg2rad(meridional_velocity / hub_tangential_speed))

        incidence_angle = 0.35 * (tip_flow_angle / 0.65)
        leading_blade_angle_tip = tip_flow_angle + incidence_angle
        leading_blade_angle_hub = np.arctan((tip_diameter / hub_diameter) * np.tan(np.deg2rad(leading_blade_angle_tip)))

        head_coeff = g * headrise / (tip_speed ** 2)


        iter = 0
        wiesner_slip = 1
        while iter < 50:
            euler_head = headrise / wiesner_slip
            outlet_swirl_velocity = euler_head * g / tip_speed
            trailing_blade_angle_tip = np.arctan(meridional_velocity / (tip_speed - outlet_swirl_velocity))
            wiesner_slip = 1 - (np.sqrt(np.sin(trailing_blade_angle_tip))) / blade_number ** 0.7

            iter = iter + 1

        return(suction_specific_speed, suction_specific_angular_speed, cavitation_number, flow_coeff_optimal, head_coeff, euler_head, outlet_swirl_velocity, trailing_blade_angle_tip, wiesner_slip, hub_flow_angle)

    [suction_specific_speed_ind, suction_specific_angular_speed_ind, cavitation_number_ind, flow_coeff_optimal_ind, head_coeff_ind, euler_head_ind, outlet_swirl_velocity_ind, trailing_blade_angle_tip_ind, wiesner_slip_ind, 
     hub_flow_angle_ind] = inducer(vol_flow, target_NPSHr, safety_factor, shaft_speed_rpm, headrise)


    def lox_pump(shaft_speed_rpm, headrise, break_hp, vol_flow, target_NPSHr, density):

        #-------------- Hardware Params (Yaml imports) -------------------------------------------------------------
        blade_angle = 27
        blade_count = 6
        #circumfrential_blockage = 3
        #effective_tip_circumference = 3
        tip_radius = 3
        #exit_height = 3
        #hub_diam = 3
        #inner_blade_angle = 3
        #outer_blade_angle = 3
        #inner_depth = 3
        #outer_depth = 3
        #-----------------------------------------------------------------------------------------
    
        outer_radius = tip_radius
        #inner_radius = hub_diam / 2    
        shaft_speed_rad = shaft_speed_rpm * ((2 * np.pi) / 60)
        g = 9.81 # m / s2
    



        #headrise =  shaft_speed_rad**2 * (outer_radius) * (outer_radius - ((inner_radius**2 * inner_depth / (outer_radius * outer_depth)) * np.tan(np.deg2rad(inner_blade_angle)) * (1/np.tan(np.deg2rad(outer_blade_angle))))) / g # meters
        hydraulic_power_watts = density * g * headrise * vol_flow # Watts
 
        #---------------------------------------------------------------------------------

        specific_speed = shaft_speed_rad * np.sqrt(vol_flow) / (headrise ** (3/4))

        tip_speed = shaft_speed_rad * tip_radius
        suction_specific_speed = shaft_speed_rad * np.sqrt(vol_flow) / (target_NPSHr ** (3/4))
        hydraulic_power_hp = hydraulic_power_watts * 745.7 # 745.7 is a conversion factor from watts to hp

        #For testing efficiency is changed t0 100%
        overall_efficiency = hydraulic_power_hp / break_hp # pg. 194 Huzzle and Huang
        overall_efficiency = 1

        head_coeff = headrise * g / (tip_speed**2) # unitless Ansys https://ansyshelp.ansys.com/public/account/secured?returnurl=/views/secured/corp/v251/en/wb_bm/bm_vista_cpd_gui_componentcontrols.html
        flow_coeff = vol_flow / (0.5 * shaft_speed_rad* (outer_radius * 2)**3) # Unitless Ansys https://ansyshelp.ansys.com/public/account/secured?returnurl=/views/secured/corp/v251/en/wb_bm/bm_vista_cpd_gui_componentcontrols.html

        meridional_velocity = tip_speed * flow_coeff

        weisner_slip = 1 - (np.sqrt(np.sin(np.deg2rad(blade_angle))) / (blade_count ** 0.7))
        tangential_velocity = (meridional_velocity * (np.tan(np.deg2rad(blade_angle)))) + (tip_speed * (1 - weisner_slip))
        abs_total_velocity = np.sqrt(meridional_velocity**2 + tangential_velocity**2)

        max_head = (tip_speed ** 2) / g
        ideal_head = tangential_velocity * tip_speed / g
        actual_head = ideal_head * overall_efficiency
        theoretical_head = headrise / overall_efficiency


        discharge_tip_radius = tip_speed * 2 / shaft_speed_rad
        specific_diam = discharge_tip_radius * (actual_head ** (1/4)) / np.sqrt(vol_flow)
        #inlet_diam_m = discharge_tip_radius / specific_diam
        #inlet_diam_cm = inlet_diam_m / 100 

        return(specific_speed, tip_speed, suction_specific_speed, hydraulic_power_hp, overall_efficiency, head_coeff, flow_coeff, meridional_velocity, weisner_slip, tangential_velocity, abs_total_velocity, max_head, ideal_head, actual_head, theoretical_head, specific_diam)
    [specific_speed_imp, tip_speed_imp, suction_specific_speed_imp, hydraulic_power_hp_imp, overall_efficiency, head_coeff_imp, flow_coeff_imp, meridional_velocity_imp, weisner_slip_imp, tangential_velocity_imp, abs_total_velocity_imp, max_head_imp, ideal_head_imp, actual_head_imp, theoretical_head_imp, specific_diam_imp] = lox_pump(shaft_speed_rpm, headrise, break_hp, vol_flow, target_NPSHr, density) 


    return(vol_flow, target_NPSHr, suction_specific_speed_ind, suction_specific_angular_speed_ind, cavitation_number_ind,
  flow_coeff_optimal_ind, head_coeff_ind, euler_head_ind, outlet_swirl_velocity_ind, trailing_blade_angle_tip_ind, 
  wiesner_slip_ind, hub_flow_angle_ind, specific_speed_imp, tip_speed_imp, suction_specific_speed_imp, hydraulic_power_hp_imp, overall_efficiency, 
 head_coeff_imp, flow_coeff_imp, meridional_velocity_imp, weisner_slip_imp, tangential_velocity_imp, abs_total_velocity_imp, 
 max_head_imp, ideal_head_imp, actual_head_imp, theoretical_head_imp, 
 specific_diam_imp)
