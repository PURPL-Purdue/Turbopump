%% =========================================================================
%  PRESSURE VESSEL ANALYSIS — Pipe with Custom Waterjet Plate Flanges
%
%  Applicable Codes:
%    • ASME B31.3        — Process Piping (pipe wall thickness, MAWP)
%    • ASME BPVC Sec. VIII Div. 1, Appendix 2  — Flat plate flange design
%    • ASME BPVC Sec. VIII Div. 1, UW-12/UW-18 — Weld joint requirements
%    • ASME BPVC Sec. II Part D               — Material allowables
%
%  Flange Type: CUSTOM FLAT PLATE (Loose-type, no hub taper)
%    — Waterjet cut from plate stock, then fillet-welded to pipe
%    — B16.5 P-T ratings do NOT apply
%    — Weld throat must be independently checked
%    — Hub factors: g0 = g1 = t_pipe (no taper)
% =========================================================================
clc; clear; close all;

fprintf('=================================================================\n');
fprintf('  PRESSURE VESSEL ANALYSIS: Pipe + Custom Waterjet Plate Flanges\n');
fprintf('  ASME B31.3 / ASME BPVC Sec. VIII Div.1 App.2 / UW-12/UW-18\n');
fprintf('=================================================================\n\n');

%% ── 1. INPUT PARAMETERS ─────────────────────────────────────────────────
%        *** EDIT THIS SECTION FOR YOUR SPECIFIC DESIGN ***

% -------------------------------------------------------------------------
% 1a. Pipe Geometry
% -------------------------------------------------------------------------
OD     = 6.625;     % Pipe outer diameter              [in]   (e.g. 10" NPS)
t_nom  = 0.719;     % Pipe nominal wall thickness      [in]
mill   = 0.125;     % Mill undertolerance fraction     [-]    (12.5% per ASTM)
CA     = 0.125;     % Corrosion allowance              [in]

% -------------------------------------------------------------------------
% 1b. Design Conditions
% -------------------------------------------------------------------------
P      = 500;       % Internal design pressure         [psi]
T_des  = 4040;       % Design temperature               [°F]

% -------------------------------------------------------------------------
% 1c. Pipe Material  (e.g. ASTM A106 Gr.B — from ASME Sec. II Part D)
% -------------------------------------------------------------------------
S_pipe  = 6500;    % Pipe allowable stress @ T_des    [psi]
E_joint = 0.7;      % Weld joint quality factor        [-]  (1.0 = seamless)
Y_coeff = 0.7;      % Temp coefficient (B31.3 §304)    [-]  (< 900°F steel)

% -------------------------------------------------------------------------
% 1d. Custom Plate Flange Geometry
%     (Flat plate, waterjet cut — NO tapered hub)
% -------------------------------------------------------------------------
A_f    = 9.4;     % Flange OD                        [in]
t_f    = 0.5;     % Flange plate thickness           [in]   ← design variable
C_f    = 0.4687;     % Bolt circle diameter             [in]
B_f    = OD;        % Flange bore ≈ pipe OD (slip-over) [in]
%  NOTE: For a plate flange welded to the outside of the pipe:
%        B_f = OD.  If the flange slips over and welds at the bore,
%        set B_f = OD - 2*t_nom (pipe ID).  Adjust as needed.

% Plate flange: no hub taper → g0 = g1 = pipe wall thickness
g0     = t_nom;     % Hub small-end thickness = pipe wall [in]
g1     = t_nom;     % Hub large-end thickness = pipe wall [in]
%  (g0 = g1 signals "no hub" to the App.2 formulation)

% -------------------------------------------------------------------------
% 1e. Flange Material  (e.g. ASTM A36 or A516 Gr.70 plate)
% -------------------------------------------------------------------------
S_f    = 6500;     % Flange plate allowable stress    [psi]  (from ASME Sec. II Part D)

% -------------------------------------------------------------------------
% 1f. Bolting
% -------------------------------------------------------------------------
n_b    = 12;        % Number of bolts                  [-]
d_b    = 0.34;     % Bolt nominal (root) diameter     [in]
S_b    = 25000;     % Bolt allowable stress            [psi]  (SA-193 B7)

% -------------------------------------------------------------------------
% 1g. Gasket  (e.g. Spiral-Wound ASME B16.20)
% -------------------------------------------------------------------------
b_g      = 0.25;    % Effective gasket seating width   [in]
G_mean   = 9.4;   % Gasket mean (reaction) diameter  [in]
m_g      = 2.0;     % Gasket factor m                  [-]    (ASME App. 2 Table 2-5.1)
y_g      = 10000;   % Min. seating stress y            [psi]  (ASME App. 2 Table 2-5.1)

% -------------------------------------------------------------------------
% 1h. Weld Parameters (Pipe-to-Plate Fillet Weld)
%     Category D weld per ASME VIII UW-18
% -------------------------------------------------------------------------
weld_leg   = 0.500; % Fillet weld leg size             [in]   (each side)
E_weld_cat = 0.70;  % Category D fillet weld efficiency[-]    (UW-12, Table UW-12)
%  Full-penetration groove weld → set E_weld_cat = 1.0

%% ── 2. PIPE WALL THICKNESS  (ASME B31.3 §304.1.2) ──────────────────────

fprintf('--- 2. PIPE WALL THICKNESS  (ASME B31.3 §304.1.2) ---\n');

t_req      = (P * OD) / (2 * (S_pipe * E_joint + P * Y_coeff));
t_req_CA   = t_req + CA;
t_avail    = t_nom * (1 - mill);     % Net available after mill tolerance
ID         = OD - 2 * t_nom;        % Pipe inner diameter
t_eff      = t_avail - CA;          % Effective thickness for pressure

fprintf('  Pipe OD  = %.3f in  |  ID = %.3f in  |  t_nom = %.3f in\n', OD, ID, t_nom);
fprintf('  Required thickness (pressure only):  t_calc  = %.4f in\n', t_req);
fprintf('  Required thickness (+ CA):           t_req   = %.4f in\n', t_req_CA);
fprintf('  Available thickness (nom - %.0f%% mill): t_avail = %.4f in\n', mill*100, t_avail);

pass_wall = t_avail >= t_req_CA;
fprintf('  Wall Thickness Check: %s  (t_avail=%.4f vs t_req=%.4f in)\n\n', ...
    pf(pass_wall), t_avail, t_req_CA);

%% ── 3. PIPE PRINCIPAL STRESSES ──────────────────────────────────────────

fprintf('--- 3. PIPE PRINCIPAL STRESSES ---\n');

sigma_h = (P * OD) / (2 * t_avail);          % Hoop — Barlow's formula  [psi]
sigma_L = (P * OD) / (4 * t_avail);          % Longitudinal (closed-end)[psi]
sigma_r = -P / 2;                             % Radial (avg through wall)[psi]
sigma_vm = sqrt(sigma_h^2 - sigma_h*sigma_L + sigma_L^2);  % von Mises  [psi]

fprintf('  Hoop stress          sigma_h  = %9.1f psi\n', sigma_h);
fprintf('  Longitudinal stress  sigma_L  = %9.1f psi\n', sigma_L);
fprintf('  Radial stress (avg)  sigma_r  = %9.1f psi\n', sigma_r);
fprintf('  von Mises equivalent sigma_vm = %9.1f psi\n', sigma_vm);
fprintf('  Allowable stress     S_pipe   = %9.1f psi\n', S_pipe);

pass_hoop = sigma_h <= S_pipe;
fprintf('  Hoop Stress Check: %s\n\n', pf(pass_hoop));

%% ── 4. MAXIMUM ALLOWABLE WORKING PRESSURE (MAWP) ────────────────────────

fprintf('--- 4. MAWP — PIPE (ASME B31.3) ---\n');

MAWP_pipe = (2 * S_pipe * E_joint * t_eff) / (OD - 2 * Y_coeff * t_eff);
margin_p  = 100 * (MAWP_pipe - P) / P;

fprintf('  Effective thickness   t_eff    = %.4f in  (t_avail - CA)\n', t_eff);
fprintf('  MAWP (pipe)           MAWP     = %.1f psi\n', MAWP_pipe);
fprintf('  Design Pressure       P        = %.1f psi\n', P);
fprintf('  Pressure Margin                = %.1f psi  (%.1f%%)\n', MAWP_pipe-P, margin_p);

pass_mawp = MAWP_pipe >= P;
fprintf('  MAWP Check: %s\n\n', pf(pass_mawp));

%% ── 5. CUSTOM FLAT PLATE FLANGE — SHAPE CONSTANTS (ASME App. 2) ─────────

fprintf('--- 5. FLAT PLATE FLANGE SHAPE CONSTANTS (ASME VIII App. 2) ---\n');
fprintf('  [Loose-type, no hub taper: g0 = g1 = t_pipe]\n\n');

% Diameter ratio
K = A_f / B_f;
fprintf('  K = A_f/B_f = %.4f\n', K);

% App.2 shape constants (Appendix 2, Tables 2-7.1 through 2-7.3)
T_K = (K^2*(1 + 8.55246*log10(K)) - 1) / ((1.04720 + 1.9448*K^2)*(K - 1));
U_K = (K^2*(1 + 8.55246*log10(K)) - 1) / (1.36136*(K^2 - 1)*(K - 1));
Y_K = (1.0/(K-1)) * (0.66845 + 5.71690*((K^2*log10(K))/(K^2 - 1)));
Z_K = (K^2 + 1) / (K^2 - 1);

fprintf('  T = %.5f\n', T_K);
fprintf('  U = %.5f\n', U_K);
fprintf('  Y = %.5f\n', Y_K);
fprintf('  Z = %.5f\n\n', Z_K);

% Hub factor h0 — for flat plate (no hub): h0 = sqrt(B_f * g0)
h0   = sqrt(B_f * g0);
fprintf('  Hub factor h0 = sqrt(B_f*g0) = %.4f in\n', h0);

% For flat plate (g0=g1), hub correction factor F=1, V=1, f=1 (App.2 simplified)
F_hub = 1.0;   % Hub bending correction (flat plate → 1)
V_hub = 0.550; % From App.2 chart for g1/g0=1, h/h0≈0 (conservative)
f_hub = 1.0;   % Hub stress correction factor (flat plate)

fprintf('  Hub correction factors (flat plate, g1/g0=1.0):\n');
fprintf('    F = %.3f  |  V = %.3f  |  f = %.3f\n\n', F_hub, V_hub, f_hub);

%% ── 6. GASKET & BOLT LOADS (ASME App. 2) ────────────────────────────────

fprintf('--- 6. GASKET & BOLT LOADS (ASME VIII App. 2) ---\n');

G = G_mean;     % Gasket reaction diameter [in]

% Hydrostatic end force
H    = (pi/4) * G^2 * P;            % [lbf]

% Gasket compression force (operating)
H_p  = pi * G * b_g * 2 * m_g * P; % [lbf]

% Required bolt load — operating condition
W_op = H + H_p;                     % [lbf]

% Required bolt load — seating condition
W_seat = pi * G * b_g * y_g;        % [lbf]

% Available bolt area
A_b_total = n_b * (pi/4) * d_b^2;  % [in²]

% Required bolt area (governing condition)
W_gov   = max(W_op, W_seat);
A_b_req = W_gov / S_b;

% Design bolt load (App.2: average of available and required)
W_design = (A_b_total * S_b + W_gov) / 2;

fprintf('  Hydrostatic end force        H       = %10.1f lbf\n', H);
fprintf('  Gasket compression (op.)     H_p     = %10.1f lbf\n', H_p);
fprintf('  Required bolt load (op.)     W_op    = %10.1f lbf\n', W_op);
fprintf('  Required bolt load (seat.)   W_seat  = %10.1f lbf\n', W_seat);
fprintf('  Governing bolt load          W_gov   = %10.1f lbf  [%s controls]\n', ...
    W_gov, ternary(W_op >= W_seat, 'Operating', 'Seating'));
fprintf('  Required bolt area           A_b_req = %10.4f in²\n', A_b_req);
fprintf('  Available bolt area          A_b_avail=%10.4f in²  (%d x %.3f" dia)\n', ...
    A_b_total, n_b, d_b);
fprintf('  Design bolt load (App.2)     W_design= %10.1f lbf\n', W_design);

pass_bolt = A_b_total >= A_b_req;
fprintf('  Bolt Area Check: %s\n\n', pf(pass_bolt));

%% ── 7. FLANGE MOMENTS (ASME App. 2) ─────────────────────────────────────

fprintf('--- 7. FLANGE MOMENTS (ASME VIII App. 2) ---\n');

% Moment arms
hD = (C_f - B_f) / 2;          % Arm for pressure on bore area   [in]
hG = (C_f - G)   / 2;          % Arm for gasket load             [in]
hT = (hD + hG)   / 2;          % Arm for pressure on flange face [in]

% Force components — operating
HD  = (pi/4) * B_f^2 * P;      % Pressure force on bore area [lbf]
HT  = H - HD;                  % Pressure on flange annular face [lbf]
HG  = W_op - H;                % Net gasket load             [lbf]

% Operating moment
M_op = abs(HD*hD) + abs(HT*hT) + abs(HG*hG);

% Seating moment (W_seat acts at hG)
M_seat = W_seat * hG;

% Governing moment
M_gov = max(M_op, M_seat);

fprintf('  Moment arms:\n');
fprintf('    hD = %.4f in  (bore pressure arm)\n', hD);
fprintf('    hG = %.4f in  (gasket load arm)\n',   hG);
fprintf('    hT = %.4f in  (face pressure arm)\n', hT);
fprintf('  Force components (operating):\n');
fprintf('    HD = %9.1f lbf\n', HD);
fprintf('    HT = %9.1f lbf\n', HT);
fprintf('    HG = %9.1f lbf\n', HG);
fprintf('  Operating moment  M_op   = %12.1f in·lbf\n', M_op);
fprintf('  Seating moment    M_seat = %12.1f in·lbf\n', M_seat);
fprintf('  Governing moment  M_gov  = %12.1f in·lbf  [%s controls]\n\n', ...
    M_gov, ternary(M_op >= M_seat, 'Operating', 'Seating'));

%% ── 8. FLAT PLATE FLANGE STRESSES (ASME App. 2 — Loose Type) ────────────

fprintf('--- 8. FLAT PLATE FLANGE STRESS CALCULATIONS (App. 2, Loose) ---\n');
fprintf('  [No hub: S_H treated as bore tangential stress only]\n\n');

% For a LOOSE-TYPE flat plate flange (App. 2, 2-3):
%   - S_H (hub) is not meaningful — flange is not integral with the vessel
%   - Primary checks are S_R and S_T in the plate ring
%   - S_H limit becomes the tangential bore stress check

% Radial stress in flange ring
S_R = (1.33 * t_f * E_joint + 1) * M_gov / (t_f^2 * B_f);
%  Note: For flat plate (loose), simplified as:
S_R = 6 * M_gov / (t_f^2 * pi * B_f);       % Radial bending at bore [psi]

% Tangential stress in flange ring
S_T = Y_K * M_gov / (t_f^2 * B_f) - Z_K * S_R;

% Longitudinal hub stress (flat plate: no taper → conservative estimate)
%  For loose type: S_H ≤ allowable (check as bending at weld junction)
S_H = f_hub * V_hub * M_gov / (U_K * g0^2 * B_f);

% Combined average stresses (App.2 combination rules)
S_HR = (S_H + S_R) / 2;
S_HT = (S_H + S_T) / 2;

fprintf('  Flat plate flange stresses:\n');
fprintf('    Longitudinal/hub stress   S_H = %10.1f psi\n', S_H);
fprintf('    Radial ring stress        S_R = %10.1f psi\n', S_R);
fprintf('    Tangential ring stress    S_T = %10.1f psi\n', S_T);
fprintf('    (S_H + S_R)/2            S_HR = %9.1f psi\n', S_HR);
fprintf('    (S_H + S_T)/2            S_HT = %9.1f psi\n\n', S_HT);

% ASME App. 2 stress limits (loose type)
lim_H  = 1.5 * S_f;
lim_R  = S_f;
lim_T  = S_f;
lim_HR = S_f;
lim_HT = S_f;

pass_SH  = S_H  <= lim_H;
pass_SR  = S_R  <= lim_R;
pass_ST  = S_T  <= lim_T;
pass_SHR = S_HR <= lim_HR;
pass_SHT = S_HT <= lim_HT;

fprintf('  ASME App.2 Stress Checks (Loose-Type Flat Plate):\n');
fprintf('    S_H  ≤ 1.5·Sf  : %s  (%8.1f ≤ %8.1f psi)\n', pf(pass_SH),  S_H,  lim_H);
fprintf('    S_R  ≤ Sf      : %s  (%8.1f ≤ %8.1f psi)\n', pf(pass_SR),  S_R,  lim_R);
fprintf('    S_T  ≤ Sf      : %s  (%8.1f ≤ %8.1f psi)\n', pf(pass_ST),  S_T,  lim_T);
fprintf('    S_HR ≤ Sf      : %s  (%8.1f ≤ %8.1f psi)\n', pf(pass_SHR), S_HR, lim_HR);
fprintf('    S_HT ≤ Sf      : %s  (%8.1f ≤ %8.1f psi)\n\n',pf(pass_SHT), S_HT, lim_HT);

%% ── 9. FLAT PLATE THICKNESS — ALTERNATE DIRECT CHECK ────────────────────
%       ASME VIII Div.1 UG-34: Flat heads and covers (conservative backup)

fprintf('--- 9. FLAT PLATE THICKNESS CHECK  (ASME VIII UG-34) ---\n');

% UG-34 flat plate formula:  t = d * sqrt(C*P/S)
% C = 0.3 for welded flat plate with full-pen weld (conservative)
C_plate = 0.3;
d_plate = B_f;          % Unsupported span = bore diameter [in]
t_req_ug34 = d_plate * sqrt(C_plate * P / S_f);

fprintf('  UG-34 flat plate formula:  t = d*sqrt(C*P/S)\n');
fprintf('    C (weld factor)  = %.2f  (full-pen groove weld, UG-34(c)(2))\n', C_plate);
fprintf('    d (bore dia.)    = %.3f in\n', d_plate);
fprintf('    Required t_ug34  = %.4f in\n', t_req_ug34);
fprintf('    Actual t_f       = %.4f in\n', t_f);

pass_ug34 = t_f >= t_req_ug34;
fprintf('  UG-34 Thickness Check: %s\n\n', pf(pass_ug34));

%% ── 10. PIPE-TO-FLANGE WELD CHECK (ASME VIII UW-18) ─────────────────────

fprintf('--- 10. PIPE-TO-PLATE WELD CHECK  (ASME VIII UW-18) ---\n');
fprintf('  [Category D fillet weld — two-sided]\n\n');

% Fillet weld throat
weld_throat = weld_leg / sqrt(2);       % Effective throat [in]

% Total weld throat area (both sides, circumferential)
A_weld_total = 2 * pi * OD/2 * weld_throat;   % [in²]

% Shear force on weld from internal pressure (longitudinal force)
F_long = (pi/4) * ID^2 * P;            % Longitudinal pressure force [lbf]

% Shear force from bolt load (weld must transfer flange moment reaction)
%  Simplified: weld carries the net unbalanced load = W_gov - H
F_gasket = max(0, W_gov - H);          % [lbf]

% Governing axial load on weld
F_weld = max(F_long, F_gasket);        % [lbf]

% Shear stress on weld throat
tau_weld = F_weld / A_weld_total;      % [psi]

% Allowable weld shear stress (ASME UW-18 Category D)
tau_allow = E_weld_cat * S_pipe;       % [psi]

fprintf('  Weld leg size              = %.4f in\n', weld_leg);
fprintf('  Effective weld throat      = %.4f in\n', weld_throat);
fprintf('  Total weld throat area     = %.4f in²  (both sides)\n', A_weld_total);
fprintf('  Longitudinal pressure force F_long   = %9.1f lbf\n', F_long);
fprintf('  Net gasket/bolt force       F_gasket  = %9.1f lbf\n', F_gasket);
fprintf('  Governing weld axial load   F_weld   = %9.1f lbf\n', F_weld);
fprintf('  Weld shear stress           tau_weld = %9.1f psi\n', tau_weld);
fprintf('  Allowable weld shear stress tau_allow= %9.1f psi  (E=%.2f × S_pipe)\n', ...
    tau_allow, E_weld_cat);

pass_weld = tau_weld <= tau_allow;
fprintf('  Weld Shear Check: %s\n\n', pf(pass_weld));

% Minimum weld size check (per ASME D1.1 / VIII general guidance)
weld_min = min(t_f, t_nom) * 0.7;     % Rule of thumb: 70% of thinner part
fprintf('  Minimum recommended weld leg = %.4f in  (0.7 × thinner part)\n', weld_min);
pass_weld_size = weld_leg >= weld_min;
fprintf('  Minimum Weld Size Check: %s  (%.4f in provided vs %.4f in min)\n\n', ...
    pf(pass_weld_size), weld_leg, weld_min);

%% ── 11. BOLT STRESS & GASKET SEATING VERIFICATION ───────────────────────

fprintf('--- 11. BOLT & GASKET SEATING VERIFICATION ---\n');

% Actual bolt stress (operating)
sigma_b_op   = W_op   / A_b_total;
sigma_b_seat = W_seat / A_b_total;
sigma_b_gov  = W_gov  / A_b_total;

fprintf('  Bolt stress — operating  : %8.1f psi\n', sigma_b_op);
fprintf('  Bolt stress — seating    : %8.1f psi\n', sigma_b_seat);
fprintf('  Bolt stress — governing  : %8.1f psi\n', sigma_b_gov);
fprintf('  Bolt allowable stress S_b: %8.1f psi\n', S_b);

pass_bolt_stress = sigma_b_gov <= S_b;
fprintf('  Bolt Stress Check: %s\n\n', pf(pass_bolt_stress));

% Gasket compression check
sigma_g_op   = W_op / (pi * G * b_g);
sigma_g_seat = W_seat / (pi * G * b_g);
fprintf('  Gasket stress — operating: %8.1f psi  (min required: %.1f psi = 2·m·P)\n', ...
    sigma_g_op, 2*m_g*P);
fprintf('  Gasket stress — seating  : %8.1f psi  (min required: %.1f psi = y_g)\n', ...
    sigma_g_seat, y_g);

pass_gasket_op   = sigma_g_op   >= 2*m_g*P;
pass_gasket_seat = sigma_g_seat >= y_g;
fprintf('  Gasket Seating Check (op.)  : %s\n', pf(pass_gasket_op));
fprintf('  Gasket Seating Check (seat.): %s\n\n', pf(pass_gasket_seat));

%% ── 12. SUMMARY TABLE ────────────────────────────────────────────────────

all_pass = pass_wall && pass_hoop && pass_mawp && pass_bolt && ...
           pass_SH && pass_SR && pass_ST && pass_SHR && pass_SHT && ...
           pass_ug34 && pass_weld && pass_weld_size && ...
           pass_bolt_stress && pass_gasket_op && pass_gasket_seat;

fprintf('=================================================================\n');
fprintf('                      RESULTS SUMMARY\n');
fprintf('=================================================================\n');
fprintf('  PIPE\n');
fprintf('  %-44s %8.4f in   [%s]\n', 'Wall thickness (required, incl. CA):', t_req_CA,   pf(pass_wall));
fprintf('  %-44s %8.4f in\n',        'Wall thickness (available):',          t_avail);
fprintf('  %-44s %8.1f psi  [%s]\n', 'Hoop stress:',                         sigma_h,    pf(pass_hoop));
fprintf('  %-44s %8.1f psi  [%s]\n', 'MAWP (pipe):',                         MAWP_pipe,  pf(pass_mawp));
fprintf('\n  FLAT PLATE FLANGE (App. 2 — Loose Type)\n');
fprintf('  %-44s %8.1f psi  [%s]\n', 'S_H (hub/bore stress):',               S_H,        pf(pass_SH));
fprintf('  %-44s %8.1f psi  [%s]\n', 'S_R (radial ring stress):',            S_R,        pf(pass_SR));
fprintf('  %-44s %8.1f psi  [%s]\n', 'S_T (tangential ring stress):',        S_T,        pf(pass_ST));
fprintf('  %-44s %8.1f psi  [%s]\n', '(S_H+S_R)/2:',                         S_HR,       pf(pass_SHR));
fprintf('  %-44s %8.1f psi  [%s]\n', '(S_H+S_T)/2:',                         S_HT,       pf(pass_SHT));
fprintf('  %-44s %8.4f in   [%s]\n', 'Flange thickness (UG-34 check):',      t_f,        pf(pass_ug34));
fprintf('\n  GOVERNING LOADS\n');
fprintf('  %-44s %8.1f lbf\n',       'Governing bolt load (W_gov):',         W_gov);
fprintf('  %-44s %8.1f in·lbf\n',   'Governing flange moment (M_gov):',     M_gov);
fprintf('\n  PIPE-TO-FLANGE WELD  (UW-18 Category D)\n');
fprintf('  %-44s %8.1f psi  [%s]\n', 'Weld shear stress:',                   tau_weld,   pf(pass_weld));
fprintf('  %-44s %8.4f in   [%s]\n', 'Weld leg size:',                       weld_leg,   pf(pass_weld_size));
fprintf('\n  BOLTING & GASKET\n');
fprintf('  %-44s %8.1f psi  [%s]\n', 'Governing bolt stress:',               sigma_b_gov,pf(pass_bolt_stress));
fprintf('  %-44s %8.1f psi  [%s]\n', 'Gasket stress (operating):',           sigma_g_op, pf(pass_gasket_op));
fprintf('  %-44s %8.1f psi  [%s]\n', 'Gasket stress (seating):',             sigma_g_seat,pf(pass_gasket_seat));
fprintf('=================================================================\n');
fprintf('  OVERALL RESULT: %s\n', ternary(all_pass, '*** ALL CHECKS PASSED ***', '!!! ONE OR MORE CHECKS FAILED !!!'));
fprintf('=================================================================\n\n');

%% ── 13. VISUALIZATION ───────────────────────────────────────────────────

figure('Name','Pressure Vessel Analysis — Custom Plate Flange', ...
       'Color','white','Position',[80 60 1400 820]);

theta_c = linspace(0, 2*pi, 400);

% ── Panel 1: Pipe principal stresses ─────────────────────────────────────
subplot(3,3,1);
vals_s  = [sigma_h, sigma_L, abs(sigma_r), sigma_vm];
lbls_s  = {'\sigma_h Hoop','\sigma_L Long.','\sigma_r Radial','\sigma_{VM}'};
clrs_s  = [0.20 0.50 0.85; 0.25 0.72 0.45; 0.85 0.45 0.20; 0.60 0.20 0.70];
b1 = bar(vals_s,'FaceColor','flat');
for ci = 1:4; b1.CData(ci,:) = clrs_s(ci,:); end
yline(S_pipe,'r--','LineWidth',1.8,'Label','S_{allow}','LabelHorizontalAlignment','left');
set(gca,'XTickLabel',lbls_s,'FontSize',9); ylabel('Stress [psi]');
title('Pipe Principal Stresses'); grid on;

% ── Panel 2: Wall thickness comparison ───────────────────────────────────
subplot(3,3,2);
tw = [t_req, t_req_CA, t_avail, t_nom, t_eff];
nw = {'t_{calc}','t_{req}+CA','t_{avail}','t_{nom}','t_{eff}'};
b2 = bar(tw,'FaceColor','flat');
b2.CData = [0.6 0.6 0.9; 0.9 0.6 0.2; 0.2 0.8 0.3; 0.75 0.75 0.75; 0.3 0.7 0.7];
set(gca,'XTickLabel',nw,'FontSize',9); ylabel('Thickness [in]');
title('Wall Thickness Breakdown'); grid on;

% ── Panel 3: Bolt load components ────────────────────────────────────────
subplot(3,3,3);
pie([H, H_p, max(0, W_seat-W_op)], ...
    {'Hydrostatic H', 'Gasket H_p', 'Seat. excess'});
title(sprintf('Bolt Load Components\nW_{gov} = %.0f lbf', W_gov));

% ── Panel 4: Flange stress checks vs allowables ───────────────────────────
subplot(3,3,4);
fl_actual = [S_H,  S_R,  S_T,  S_HR, S_HT];
fl_limits = [lim_H,lim_R,lim_T,lim_HR,lim_HT];
fl_labels = {'S_H','S_R','S_T','S_{HR}','S_{HT}'};
x_fl = 1:5;
bar(x_fl, fl_limits, 'FaceColor',[0.90 0.85 0.85],'EdgeColor','none'); hold on;
bar(x_fl, fl_actual, 0.5,'FaceColor','flat');
colormap(gca, cool(5));
set(gca,'XTickLabel',fl_labels,'FontSize',9);
ylabel('Stress [psi]');
legend('Limit (ASME)','Actual','Location','NorthEast');
title('Flange Stress Checks (App. 2 Loose)'); grid on; hold off;

% ── Panel 5: Weld stress check ────────────────────────────────────────────
subplot(3,3,5);
bar_vals = [tau_weld, tau_allow];
bar_labs = {'Weld \tau_{actual}','Allowable \tau_{allow}'};
bw = bar(bar_vals,'FaceColor','flat');
bw.CData(1,:) = ternary_color(pass_weld);
bw.CData(2,:) = [0.75 0.75 0.75];
set(gca,'XTickLabel',bar_labs,'FontSize',10);
ylabel('Shear Stress [psi]');
title(sprintf('Weld Shear Check (UW-18)\nLeg = %.3f in | Throat = %.3f in', ...
    weld_leg, weld_throat));
grid on;

% ── Panel 6: MAWP vs Design Pressure ─────────────────────────────────────
subplot(3,3,6);
th = linspace(0,2*pi,200);
r_mawp = MAWP_pipe; r_p = P;
fill(r_mawp*cos(th), r_mawp*sin(th), [0.85 0.93 1.0],'EdgeColor',[0.5 0.5 0.8]); hold on;
fill(r_p*cos(th),    r_p*sin(th),    [0.25 0.55 0.90],'EdgeColor','none');
axis equal; axis off;
text(0, 0, sprintf('P = %g psi\nMAWP = %g psi\nMargin: %.1f%%', ...
    P, MAWP_pipe, margin_p), 'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
title('Pressure Margin (Pipe)'); hold off;

% ── Panel 7: Gasket stress check ─────────────────────────────────────────
subplot(3,3,7);
g_act  = [sigma_g_op,  sigma_g_seat];
g_req  = [2*m_g*P,     y_g];
g_labs = {'Op. \sigma_g','Seat. \sigma_g'};
bg = bar([g_act; g_req]');
bg(1).FaceColor = [0.30 0.65 0.90];
bg(2).FaceColor = [0.90 0.40 0.40];
set(gca,'XTickLabel',g_labs,'FontSize',10);
ylabel('Stress [psi]'); legend('Actual','Required','Location','NorthWest');
title('Gasket Seating Check'); grid on;

% ── Panel 8: Flange cross-section (face view) ────────────────────────────
subplot(3,3,8);
% Flange plate ring
fill(A_f/2*cos(theta_c), A_f/2*sin(theta_c), [0.65 0.65 0.70]); hold on;
fill(B_f/2*cos(theta_c), B_f/2*sin(theta_c), [0.92 0.92 0.95]);
% Pipe wall
fill(OD/2*cos(theta_c), OD/2*sin(theta_c),   [0.80 0.80 0.85]);
fill(ID/2*cos(theta_c), ID/2*sin(theta_c),   [0.97 0.97 1.00]);
% Bolt holes on PCD
for i = 1:n_b
    ang = 2*pi*(i-1)/n_b;
    bx  = C_f/2*cos(ang); by = C_f/2*sin(ang);
    fill(bx + d_b/2*cos(theta_c), by + d_b/2*sin(theta_c),[0.20 0.20 0.20]);
end
% Weld indication (orange arc)
weld_r = OD/2 + weld_leg*0.5;
plot(weld_r*cos(theta_c), weld_r*sin(theta_c), '-', ...
    'Color',[0.95 0.55 0.10],'LineWidth',3);
axis equal; axis off;
title(sprintf('Plate Flange — Face View\nOD=%.2f"  FlangeOD=%.2f"  %d bolts  t_f=%.3f"', ...
    OD, A_f, n_b, t_f));
legend({'Flange plate','Bore','Pipe wall','Bore ID','Bolts','Fillet weld'}, ...
    'Location','SouthOutside','FontSize',7,'NumColumns',3);

% ── Panel 9: Pass/Fail summary scorecard ─────────────────────────────────
subplot(3,3,9);
axis off;
checks = { 'Wall thickness (B31.3)',     pass_wall;
           'Hoop stress',                pass_hoop;
           'MAWP (pipe)',                pass_mawp;
           'Bolt area',                  pass_bolt;
           'Flange S_H',                 pass_SH;
           'Flange S_R',                 pass_SR;
           'Flange S_T',                 pass_ST;
           'Flange (S_H+S_R)/2',         pass_SHR;
           'Flange (S_H+S_T)/2',         pass_SHT;
           'Flat plate UG-34',           pass_ug34;
           'Weld shear (UW-18)',         pass_weld;
           'Weld leg size',              pass_weld_size;
           'Bolt stress',               pass_bolt_stress;
           'Gasket stress (op.)',        pass_gasket_op;
           'Gasket stress (seat.)',      pass_gasket_seat };
n_chk = size(checks,1);
for ci = 1:n_chk
    y_pos = 1 - (ci-1)/(n_chk);
    passed = checks{ci,2};
    clr = ternary_color(passed);
    rectangle('Position',[0.65, y_pos-0.04, 0.25, 0.06],'FaceColor',clr, ...
              'EdgeColor','none','Curvature',0.3);
    text(0.02, y_pos-0.01, checks{ci,1},'FontSize',8,'VerticalAlignment','middle');
    text(0.775, y_pos-0.01, ternary(passed,'PASS','FAIL'), ...
         'FontSize',8,'FontWeight','bold','Color','white', ...
         'HorizontalAlignment','center','VerticalAlignment','middle');
end
title('Check Scorecard','FontWeight','bold');

sgtitle({'Pressure Vessel Analysis — Pipe with Custom Waterjet Plate Flange', ...
         'ASME B31.3 / BPVC Sec. VIII Div.1 App.2 / UW-18'}, ...
        'FontSize',13,'FontWeight','bold');

fprintf('Visualization complete.\n');

%% ── LOCAL HELPER FUNCTIONS ───────────────────────────────────────────────

function s = pf(flag)
    % Returns 'PASS' or 'FAIL' string
    if flag; s = 'PASS'; else; s = 'FAIL'; end
end

function s = ternary(cond, a, b)
    % Returns a if cond true, else b
    if cond; s = a; else; s = b; end
end

function clr = ternary_color(flag)
    % Returns green RGB for pass, red RGB for fail
    if flag
        clr = [0.20 0.70 0.30];
    else
        clr = [0.85 0.20 0.20];
    end
end