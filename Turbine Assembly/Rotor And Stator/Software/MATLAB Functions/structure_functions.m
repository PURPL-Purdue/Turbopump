function funs = structure_functions
    % function that returns a struct of function handles
    funs.calc_radius_turbine = @calc_radius_turbine;
    funs.calc_mass_blade = @calc_mass_blade;
    funs.calc_Force_centrifugal = @calc_Force_centrifugal;
    funs.calc_stress_centrifugal = @calc_stress_centrifugal;
    funs.calc_Force_tangential = @calc_Force_tangential;
    funs.calc_Force_axial = @calc_Force_axial;
    funs.calc_torque_blade = @calc_torque_blade;
    funs.calc_torque_turbine = @calc_torque_turbine;
    funs.calc_P = @calc_P;
    funs.calc_Force_gas = @calc_Force_gas;
    funs.calc_Moment_Bending = @calc_Moment_Bending;
    funs.calc_I = @calc_I;
    funs.calc_stress_gas = @calc_stress_gas;
end

% function definitions
function radius_turbine = calc_radius_turbine(height_blade, radius_hub)
    radius_turbine = 0.5 * height_blade + radius_hub;
end

function mass_blade = calc_mass_blade(Length_blade, width_blade, height_bmin, rho_blade)
    mass_blade = Length_blade * width_blade * height_bmin * rho_blade;
end

function Force_centrifugal = calc_Force_centrifugal(mass_blade, w, radius_turbine)
    Force_centrifugal = mass_blade * (w^2) * radius_turbine;
end

function stress_centrifugal = calc_stress_centrifugal(radius_turbine, rho_blade, height_blade, w)
    stress_centrifugal = radius_turbine * rho_blade * height_blade * (w^2);
end

function Force_tangential = calc_Force_tangential(m_dot, V_1, beta_1, V_2, beta_2)
    Force_tangential = m_dot * (V_1 * cos(beta_1) + V_2 * cos(beta_2));
end

function Force_axial = calc_Force_axial(m_dot, C_1, alpha_1, V_2, beta_2)
    Force_axial = m_dot * (C_1 * sin(alpha_1) + V_2 * sin(beta_2));
end

function torque_blade = calc_torque_blade(Force_tangential, height_blade)
    torque_blade = 0.5 * Force_tangential * height_blade;
end

function torque_turbine = calc_torque_turbine(Force_tangential, radius_turbine, Z_blade)
    torque_turbine = Force_tangential * radius_turbine * Z_blade;
end

function P = calc_P(torque_turbine, w)
    P = torque_turbine * w;
end

function Force_gas = calc_Force_gas(Force_tangential, Force_axial)
    Force_gas = sqrt((Force_tangential^2) + (Force_axial^2));
end

function Moment_Bending = calc_Moment_Bending(height_blade, Z_blade, Force_gas)
    Moment_Bending = (height_blade / (2 * Z_blade)) * Force_gas;
end

function I = calc_I(Length_blade, height_bmin)
    I = (1 / 12) * (Length_blade^3) * height_bmin;
end

function stress_gas = calc_stress_gas(height_bmin, Force_gas, width_b, I)
    stress_gas = (0.5 * height_bmin * Force_gas * width_b) / I;
end
