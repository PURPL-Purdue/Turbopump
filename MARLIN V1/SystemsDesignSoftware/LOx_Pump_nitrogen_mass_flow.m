% Constants

LOx_mass = 4.2433+.9246+.187+.041;        % LOX casing mass (lbm)
h_LN2 = 200e3;            % latent heat of vaporization for LN2 (J/kg)

T_final = 90;             % target casing temperature (K)
T_initial = 300;          % initial casing temperature (K)

% Convert SS-316 mass to kg
m_LOx_kg = LOx_mass * 0.45359237;   % kg

% SS-316 specific heat formula coefficients (J/(kg-K))
a = -1879.464;
b = 3643.198;
c = 76.70125;
d = -6176.028;
e = 7437.6247;
f = -4305.7217;
g = 1382.4627;
h = -237.22704;
i = 17.05262;

% Use representative temperature for cp
T_rep = 0.5 * (T_initial + T_final);

% Calculating specific heat for SS-316
logT = log10(T_rep);
logy = a + b*logT + c*logT^2 + d*logT^3 + e*logT^4 + ...
       f*logT^5 + g*logT^6 + h*logT^7 + i*logT^8;
c_LOx = 10^logy;          % J/(kg-K)

% Required minimum LN2 mass
m_LN2_kg = m_LOx_kg * c_LOx * (T_initial - T_final) / h_LN2;

% Apply factor of safety / extra usage factor
m_LN2_kg = m_LN2_kg * 5;

% Convert back to lbm
m_LN2_lbm = m_LN2_kg * 2.20462;

fprintf('5 times the minimum required liquid nitrogen mass is %.2f lbm\n', m_LN2_lbm);