%% ========================================================================
%  2-D TRANSIENT CYLINDRICAL HEAT TRANSFER SOLVER — ROCKET NOZZLE WALL
%  -----------------------------------------------------------------------
%  Solves the 2-D unsteady conduction equation in cylindrical coordinates
%  over the rocket nozzle/combustion-chamber wall cross-section.
%
%  Governing equation (no circumferential variation assumed, i.e. ∂/∂θ = 0):
%
%    ρ·cp · ∂T/∂t = 1/r · ∂/∂r(r·k·∂T/∂r) + ∂/∂z(k·∂T/∂z)
%
%  Discretization : Finite-Volume Method (FVM), implicit backward-Euler
%                   → unconditionally stable for any dt
%
%  Linear solver  : Red-Black SOR (Successive Over-Relaxation)
%    • Grid is partitioned once into two interleaved checkerboard sets.
%      Red nodes only touch black neighbours and vice-versa, so each
%      half-sweep is fully data-independent and executes as a single
%      vectorised array operation — no inner node loop.
%    • Warm-started from the previous time step.
%    • Convergence: max|R| < solver_tol · max|B|  (global vectorised
%      residual, evaluated every res_check_freq full sweeps).
%
%  Boundary conditions:
%    • Inner wall (gas side)  : q = hG·(Taw − T)           [firing]
%                               q = h_nat_inner·(Tamb − T)  [cooling]
%    • Outer wall (cold side) : q = h_nat·(Tamb − T)
%    • Axial ends (j=1, j=m) : adiabatic (zero flux)
%
%  Performance vs. simple SOR:
%    • Coefficient assembly: fully vectorised — no inner i/j loops
%    • RB half-sweeps: two array operations per iteration, O(1) in MATLAB
%    • Static precomputation: insert mask, RB masks, face-area ratios
%      all computed once before the time loop
%    • Residual check throttled to every res_check_freq sweeps
%
%  External inputs:
%    turbopump_properties.csv  — thermodynamic properties at key stations
%    contour.csv               — (x[mm], r[mm]) nozzle contour points
%
%  Author  : Rafael Macia Titos
%  Date    : 2026
%% ========================================================================
clear; clc; close all;
%% ===== 1. USER SETTINGS ==================================================
% ---- Time integration ----
dt          = 1e-4;    % [s]  time step (implicit → stable for large dt)
runtime     = 2.0;     % [s]  engine firing duration
tf          = 10;      % [s]  total simulation time (firing + cooling)
% ---- Spatial discretization ----
n           = 100;      % number of radial control volumes 
                       % m (axial count) is set automatically from the contour
% ---- Material selection ----
%   1 = Copper (CuCrZr-type, polynomial fit)
%   2 = AISI 8630 steel (constant props)
%   3 = Graphite (constant props)
mtype1      = 2;       % base wall material
mtype_ins   = 3;       % insert material
% ---- Insert geometry ---- Changed for Test1 16:02 28/02
ins_x_start = 0.28685;   % [m]  axial start position from injector
ins_length  = 0.1016;    % [m]  axial length of the insert
ins_thick   = 0.03992;   % [m]  radial thickness measured from r_throat
% ---- Boundary conditions ----
h_nat       = 5;       % [W/m²K]  outer-wall natural convection
h_nat_inner = 10;      % [W/m²K]  inner-wall natural convection (post-firing)
T_amb       = 300;     % [K]      ambient temperature
% ---- Mass flow ----
mdot        = 9;       % [kg/s]  total propellant mass flow rate
% ---- Combustion-chamber wall thickness offset ----
cc_wall_offset = 0.03; % [m]  extra radial thickness in CC section
% ---- Visualization / video ----
plot_interval   = 20;
record_video    = true;
video_filename  = 'nozzle_heat_transfer.mp4';
video_fps       = 15;
video_quality   = 100;
% ---- Red-Black SOR solver settings ----
%  solver_tol      Relative residual tolerance.  max|R| < tol·max|B|.
%                  Loosen to 1e-4 for speed when dt is small (warm start).
%  solver_max_iter Safety ceiling on full RB sweeps per time step.
%  sor_omega       Relaxation factor ω ∈ (1, 2).
%                    ω = 1.0  → plain Gauss-Seidel
%                    ω = 1.5  → good default for heat-equation problems
%                    ω → 2.0  → aggressive; reduce to 1.3–1.4 on
%                                non-uniform or high-contrast-k grids
%  res_check_freq  Evaluate global residual every N sweeps to amortise
%                  the cost of the residual pass.
solver_tol      = 1e-6;
solver_max_iter = 500;
sor_omega       = 1.5;
res_check_freq  = 5;
%% ===== 2. LOAD GEOMETRY ==================================================
contour_tbl     = readtable('contour.csv');
[x_raw, I]      = sort((contour_tbl.x - min(contour_tbl.x)) / 1000);
r_raw           = contour_tbl.y(I) / 1000;
[~, idx_unique] = uniquetol(x_raw, 1e-12);
x_contour       = x_raw(idx_unique);
r_contour       = r_raw(idx_unique);
m               = length(x_contour);
fprintf('Geometry loaded  : %d axial stations.\n', m);
fprintf('Duplicates removed: %d points.\n', length(x_raw) - m);
%% ===== 3. KEY GEOMETRIC FEATURES =========================================
[r_throat, throat_idx] = min(r_contour);
Throat_D = 2 * r_throat;
Throat_A = pi * r_throat^2;
slope  = diff(r_contour) ./ (diff(x_contour) + 1e-9);
cc_idx = find(slope < -1e-4, 1, 'first');
if isempty(cc_idx), cc_idx = 1; end
fprintf('CC end  : idx=%d  x=%.3f m\n', cc_idx,     x_contour(cc_idx));
fprintf('Throat  : idx=%d  x=%.3f m\n', throat_idx, x_contour(throat_idx));
%% ===== 4. LOAD THERMODYNAMIC PROPERTIES ==================================
param         = readtable('turbopump_properties.csv');
param.Station = categorical(param.Station);
pos = [0, x_contour(cc_idx), x_contour(throat_idx), x_contour(end)];
iprop = @(p) [interp1(pos(1:2), p(1:2), x_contour(1:cc_idx),      'linear'); ...
              interp1(pos(2:4), p(2:4), x_contour(cc_idx+1:end),   'spline')];
gamma = iprop(param.Gamma);
Pr    = iprop(param.Prandtl);
mu    = iprop(param.Viscosity_Pa_s);
cp_g  = iprop(param.Cp_J_kgK);
%% ===== 5. MACH-NUMBER DISTRIBUTION =======================================
M                 = ones(m, 1);
M(1:throat_idx)   = 0.3;
M(throat_idx:end) = 2.0;
for iter = 1:1000
    AoA    = (((gamma+1)/2).^(-(gamma+1)./(2*(gamma-1)))) ...
            .* (1 + (gamma-1)/2.*M.^2).^((gamma+1)./(2*(gamma-1))) ./ M;
    res_M  = pi*r_contour.^2 / Throat_A - AoA;
    dM_eps = 1e-6;
    AoA_dM = (((gamma+1)/2).^(-(gamma+1)./(2*(gamma-1)))) ...
            .* (1 + (gamma-1)/2.*(M+dM_eps).^2).^((gamma+1)./(2*(gamma-1))) ./ (M+dM_eps);
    M      = M + 0.5 * res_M ./ ((AoA_dM - AoA) / dM_eps);
    if all(abs(res_M) < 1e-6), break; end
end
fprintf('Mach converged: %d iters, max res = %.2e\n', iter, max(abs(res_M)));
%% ===== 6. BARTZ CORRELATION ==============================================
Pc      = param.Pressure_Pa(param.Station == 'Injector');
Tc      = param.Temperature_K(param.Station == 'Injector');
cstar   = Pc * Throat_A / mdot;
T_g     = Tc ./ (1 + (gamma-1)/2 .* M.^2);
hG_base = (0.026 / Throat_D^0.2) ...
        .* (mu.^0.2 .* cp_g ./ Pr.^0.6) ...
        .* (Pc / cstar)^0.8 ...
        .* (Throat_A ./ (pi*r_contour.^2)).^0.9;
Taw     = T_g .* (1 + (gamma-1)/2 .* M.^2 .* Pr.^(1/3));
%% ===== 7. FVM GRID =======================================================
thickness = cc_wall_offset + r_contour(cc_idx) - r_contour;   % (m×1)
dr        = thickness(:)' / n;                                 % (1×m)
r_f       = zeros(n+1, m);
for j = 1:m
    r_f(:,j) = linspace(r_contour(j), r_contour(j)+thickness(j), n+1)';
end
rw = r_f(1:n,   :);
re = r_f(2:n+1, :);
dx_vec      = diff(x_contour);
dx_vec      = [dx_vec; dx_vec(end)];   % (m×1)
dx_vec      = dx_vec(:)';              % (1×m)
Vp = pi*(re.^2 - rw.^2) .* dx_vec;
Aw = 2*pi*rw .* dx_vec;
Ae = 2*pi*re .* dx_vec;
Ar = pi*(re.^2 - rw.^2);
% Axial conductances use centre-to-centre distances (dx_S, dx_N)
% correct for non-uniform mesh
dx_S        = zeros(1, m);
dx_N        = zeros(1, m);
dx_S(2:m)   = (dx_vec(1:m-1) + dx_vec(2:m)) / 2;
dx_N(1:m-1) = (dx_vec(1:m-1) + dx_vec(2:m)) / 2;
dx_S(1)     = dx_vec(1);     % unused — AS zeroed at axial boundary
dx_N(m)     = dx_vec(m);     % unused — AN zeroed at axial boundary
%% ===== 8. MATERIAL PROPERTY POLYNOMIALS ==================================
% Row 1: k   [W/mK],  Row 2: rho [kg/m³],
% Row 3: cp  [J/kgK], Row 4: T_max [K] (col 1 only)
materials = { ...
    [507.735, -0.540163,  7.06823e-4, -3.13884e-7; ...
       9075.0, -0.512415,  0,           0;          ...
      188.548,  1.00219,  -1.32148e-3,  5.75175e-7; ...
       1200,    0,          0,           0          ], ...
    [42.7,  0, 0, 0; 7850, 0, 0, 0;  475, 0, 0, 0; 5000, 0, 0, 0], ...
    [168,   0, 0, 0; 2500, 0, 0, 0;  717, 0, 0, 0; 5000, 0, 0, 0] ...
};
Mat1   = materials{mtype1};
MatIns = materials{mtype_ins};
% Inline polynomial evaluator (operates on full 2-D arrays)
poly4 = @(row, T) row(1) + row(2)*T + row(3)*T.^2 + row(4)*T.^3;
%% ===== 9. STATIC PRECOMPUTATIONS (computed once, reused every step) ======
% --- Insert geometry mask (logical, n×m) ---
[II, JJ]  = ndgrid(1:n, 1:m);
ins_mask  = (repmat(x_contour(:)', n, 1) >= ins_x_start)             & ...
            (repmat(x_contour(:)', n, 1) <= ins_x_start + ins_length) & ...
            (r_f(1:n,:) - r_throat       <= ins_thick)                & ...
            (r_f(1:n,:)                  >= repmat(r_contour(:)', n, 1));
base_mask = ~ins_mask;
% --- Red-Black partition masks ---
%     Red  nodes: mod(i+j,2)==0   Black nodes: mod(i+j,2)==1
%     Red nodes only touch black neighbours → half-sweep is data-independent
red_mask = mod(II + JJ, 2) == 0;
blk_mask = ~red_mask;
% --- Face-area / distance prefactors (geometry is fixed in time) ---
dr_2d    = repmat(dr,     n, 1);   % (n×m)
dx_S_2d  = repmat(dx_S,   n, 1);   % (n×m)  centre-to-centre south
dx_N_2d  = repmat(dx_N,   n, 1);   % (n×m)  centre-to-centre north
Aw_ov_dr = Aw ./ dr_2d;
Ae_ov_dr = Ae ./ dr_2d;
Ar_ov_dxS = Ar ./ dx_S_2d;        
Ar_ov_dxN = Ar ./ dx_N_2d;        
% --- Step count threshold for BC switch ---
% integer comparison avoids floating-point drift in time accumulation
firing_steps = round(runtime / dt);
fprintf('Static precomputation done.\n');
fprintf('Red nodes: %d   Black nodes: %d\n', nnz(red_mask), nnz(blk_mask));
%% ===== 10. INITIALISATION ================================================
T_2d       = 300 * ones(n, m);
time       = 0;
step_count = 0;
%% ===== 11. FIGURE & VIDEO SETUP ==========================================
fig = figure('Color', 'w', 'Position', [50, 50, 1080, 720]);
tiledlayout(2, 1, 'TileSpacing', 'compact');
ax1 = nexttile;
[Z_grid, R_norm] = meshgrid(x_contour, linspace(0, 1, n));
R_grid = zeros(size(R_norm));
for j = 1:m
    R_grid(:,j) = r_contour(j) + R_norm(:,j) * thickness(j);
end
h_plot = pcolor(ax1, Z_grid, R_grid, T_2d);
shading(ax1, 'interp');
cb = colorbar(ax1);  cb.Label.String = 'Temperature [K]';
colormap(ax1, hot);  caxis(ax1, [300, 3000]);
xlabel(ax1, 'Axial Position [m]');  ylabel(ax1, 'Radius [m]');
title(ax1, '2D Wall Temperature Distribution');
axis(ax1, 'equal', 'tight');
hold(ax1, 'on');
rectangle(ax1, 'Position', [ins_x_start, r_throat, ins_length, ins_thick], ...
    'EdgeColor', 'w', 'LineWidth', 0.5, 'LineStyle', ':');
ax2 = nexttile;
h_wall_line = plot(ax2, x_contour, T_2d(1,:), 'r-', 'LineWidth', 2);
hold(ax2, 'on');
xline(ax2, ins_x_start,            'b--', 'Start Insert', 'LabelVerticalAlignment', 'bottom');
xline(ax2, ins_x_start+ins_length, 'b--', 'End Insert',   'LabelVerticalAlignment', 'bottom');
h_taw_line = plot(ax2, x_contour, Taw, 'k:', 'LineWidth', 1.0);
grid(ax2, 'on');
xlabel(ax2, 'Axial Position [m]');  ylabel(ax2, 'Inner Wall Temperature [K]');
title(ax2, 'Wall Temperature Profile');
xlim(ax2, [min(x_contour), max(x_contour)]);  ylim(ax2, [300, 3500]);
legend(ax2, 'Inner Wall', 'Gas Recovery (Taw)', 'Location', 'northeast');
if record_video
    v_vid = VideoWriter(video_filename, 'MPEG-4');
    v_vid.FrameRate = video_fps;  v_vid.Quality = video_quality;
    open(v_vid);
end
%% ===== 12. MAIN TRANSIENT LOOP ===========================================
%
%  Each time step:
%
%  (a) VECTORISED conductivity field  k(T)
%      Material selection via precomputed ins_mask / base_mask.
%
%  (b) VECTORISED coefficient assembly
%      All five arrays (aW, aE, aS, aN, aP) and Bm assembled in one pass
%      using shifted-array indexing.  No inner i/j loops.
%      AS and AN use Ar_ov_dxS / Ar_ov_dxN (centre-to-centre).
%      BC switch uses step_count < firing_steps.
%
%  (c) RED-BLACK SOR
%      Two vectorised half-sweeps per iteration:
%        Red   half-sweep: update all red  nodes at once (black T untouched)
%        Black half-sweep: update all black nodes at once (red T already fresh)
%      Global residual evaluated with a clean vectorised pass every
%      res_check_freq full sweeps — no stale-neighbour contamination.
%
fprintf('\nStarting transient simulation (Red-Black SOR, ω = %.2f)...\n', sor_omega);
tic;
try
    while time < tf
        % ------------------------------------------------------------------
        % (a) Vectorised conductivity field  k(T)
        % ------------------------------------------------------------------
        T_c_base = min(T_2d, Mat1(4,1));
        T_c_ins  = min(T_2d, MatIns(4,1));
        k_field  = base_mask .* poly4(Mat1(1,:),   T_c_base) + ...
                   ins_mask  .* poly4(MatIns(1,:),  T_c_ins);
        % ------------------------------------------------------------------
        % (b) Vectorised coefficient assembly
        % ------------------------------------------------------------------
        % rho and cp fields
        rho_f = base_mask .* poly4(Mat1(2,:),   T_c_base) + ...
                ins_mask  .* poly4(MatIns(2,:),  T_c_ins);
        cp_f  = base_mask .* poly4(Mat1(3,:),   T_c_base) + ...
                ins_mask  .* poly4(MatIns(3,:),  T_c_ins);
        % Harmonic-mean conductivities at the four faces
        % Boundary padding replicates edge values; corresponding
        % conductances are zeroed below via BC enforcement.
        kW_nb = [k_field(1,:);     k_field(1:n-1,:)];
        kE_nb = [k_field(2:n,:);   k_field(n,:)    ];
        kS_nb = [k_field(:,1),     k_field(:,1:m-1)];
        kN_nb = [k_field(:,2:m),   k_field(:,m)    ];
        kw_h  = 2*k_field.*kW_nb ./ (k_field + kW_nb);
        ke_h  = 2*k_field.*kE_nb ./ (k_field + kE_nb);
        ks_h  = 2*k_field.*kS_nb ./ (k_field + kS_nb);
        kn_h  = 2*k_field.*kN_nb ./ (k_field + kN_nb);
        % Raw conductance arrays
        % axial conductances use centre-to-centre distances
        aW = kw_h .* Aw_ov_dr;
        aE = ke_h .* Ae_ov_dr;
        aS = ks_h .* Ar_ov_dxS;   
        aN = kn_h .* Ar_ov_dxN;   
        % Adiabatic axial boundaries and physical radial boundaries
        aW(1,:) = 0;   % inner wall: replaced by convective BC below
        aE(n,:) = 0;   % outer wall: replaced by convective BC below
        aS(:,1) = 0;   % axial left:  adiabatic
        aN(:,m) = 0;   % axial right: adiabatic
        % Gas-side BC: vectorised over all j simultaneously
        % step_count comparison instead of floating-point time
        if step_count < firing_steps
            term_Ma  = 1 + (gamma-1)/2 .* M.^2;                      % (m×1)
            sigma_v  = 1 ./ ((0.5*(T_2d(1,:)'./T_g).*term_Ma + 0.5).^0.68 ...
                             .* term_Ma.^0.12);
            hg_vec   = hG_base .* sigma_v;                            % (m×1)
            Tref_vec = Taw;                                            % (m×1)
        else
            hg_vec   = h_nat_inner * ones(m, 1);
            Tref_vec = T_amb       * ones(m, 1);
        end
        bin_row  = hg_vec(:)' .* Aw(1,:);   % (1×m)
        bout_row = h_nat      .* Ae(n,:);   % (1×m)
        % Source contributions from convective BCs
        src      = zeros(n, m);
        src(1,:) = bin_row  .* Tref_vec(:)';
        src(n,:) = src(n,:) + bout_row * T_amb;
        % Transient term
        AP0 = rho_f .* cp_f .* Vp / dt;
        % Central coefficient aP
        aP_bc      = zeros(n, m);
        aP_bc(1,:) = bin_row;
        aP_bc(n,:) = aP_bc(n,:) + bout_row;
        aP         = aW + aE + aS + aN + AP0 + aP_bc;
        % RHS
        Bm     = AP0 .* T_2d + src;
        B_norm = max(abs(Bm(:)));
        % ------------------------------------------------------------------
        % (c) Red-Black SOR
        %
        %  Each full iteration = one red half-sweep + one black half-sweep.
        %
        %  Half-sweep mechanics (shown for red; black is identical):
        %    1. Build the four neighbour arrays T_W, T_E, T_S, T_N
        %       using array slices of the current T_2d.
        %    2. Compute the Gauss-Seidel update for ALL nodes:
        %         T_gs = (Bm + aW·T_W + aE·T_E + aS·T_S + aN·T_N) / aP
        %    3. Apply SOR correction ONLY to the active colour:
        %         T_2d(red_mask) += ω·(T_gs − T_2d)(red_mask)
        %
        %  After the black half-sweep, red values are already fresh, so
        %  the black update sees the most up-to-date neighbourhood possible.
        %
        %  Global residual is computed in a separate vectorised pass
        %  every res_check_freq sweeps using current T_2d only — never mixed
        %  with mid-sweep stale values.
        % ------------------------------------------------------------------
        converged = false;
        for sor_iter = 1:solver_max_iter
            % ---- Red half-sweep ----
            T_W = [T_2d(1,:);    T_2d(1:n-1,:)];
            T_E = [T_2d(2:n,:);  T_2d(n,:)    ];
            T_S = [T_2d(:,1),    T_2d(:,1:m-1)];
            T_N = [T_2d(:,2:m),  T_2d(:,m)    ];
            T_gs           = (Bm + aW.*T_W + aE.*T_E + aS.*T_S + aN.*T_N) ./ aP;
            dT             = sor_omega * (T_gs - T_2d);
            T_2d(red_mask) = T_2d(red_mask) + dT(red_mask);
            % ---- Black half-sweep (uses freshly updated red values) ----
            T_W = [T_2d(1,:);    T_2d(1:n-1,:)];
            T_E = [T_2d(2:n,:);  T_2d(n,:)    ];
            T_S = [T_2d(:,1),    T_2d(:,1:m-1)];
            T_N = [T_2d(:,2:m),  T_2d(:,m)    ];
            T_gs           = (Bm + aW.*T_W + aE.*T_E + aS.*T_S + aN.*T_N) ./ aP;
            dT             = sor_omega * (T_gs - T_2d);
            T_2d(blk_mask) = T_2d(blk_mask) + dT(blk_mask);
            % Global residual — clean vectorised pass, throttled
            if mod(sor_iter, res_check_freq) == 0
                T_W_r = [T_2d(1,:);    T_2d(1:n-1,:)];
                T_E_r = [T_2d(2:n,:);  T_2d(n,:)    ];
                T_S_r = [T_2d(:,1),    T_2d(:,1:m-1)];
                T_N_r = [T_2d(:,2:m),  T_2d(:,m)    ];
                R_mat = aP.*T_2d - aW.*T_W_r - aE.*T_E_r ...
                                 - aS.*T_S_r - aN.*T_N_r - Bm;
                if max(abs(R_mat(:))) < solver_tol * B_norm
                    converged = true;
                    break;
                end
            end
        end
        if ~converged
            % Compute residual for the warning message if not already fresh
            if mod(sor_iter, res_check_freq) ~= 0
                T_W_r = [T_2d(1,:);    T_2d(1:n-1,:)];
                T_E_r = [T_2d(2:n,:);  T_2d(n,:)    ];
                T_S_r = [T_2d(:,1),    T_2d(:,1:m-1)];
                T_N_r = [T_2d(:,2:m),  T_2d(:,m)    ];
                R_mat = aP.*T_2d - aW.*T_W_r - aE.*T_E_r ...
                                 - aS.*T_S_r - aN.*T_N_r - Bm;
            end
            fprintf('WARNING: RB-SOR did not converge at t=%.4f s  (max|R|/|B| = %.2e, iters=%d)\n', ...
                    time, max(abs(R_mat(:))) / max(B_norm, 1e-30), sor_iter);
        end
        time       = time + dt;
        step_count = step_count + 1;
        if mod(step_count, plot_interval) == 0
            set(h_plot, 'CData', T_2d);
            state_str = 'FIRING';  if step_count >= firing_steps, state_str = 'COOLING'; end
            title(ax1, sprintf('%s | t = %.3f s', state_str, time));
            set(h_wall_line, 'YData', T_2d(1,:));
            if step_count >= firing_steps
                set(h_taw_line, 'YData', T_amb * ones(size(x_contour)));
            end
            drawnow;
            if record_video, writeVideo(v_vid, getframe(fig)); end
        end
    end
catch ME
    % Full stack trace so the error line is never lost
    fprintf('\nSimulation error at t = %.4f s:\n', time);
    fprintf('%s\n', getReport(ME, 'extended'));
    if record_video
        try, close(v_vid); catch, end
    end
    rethrow(ME);
end
elapsed = toc;
fprintf('\nSimulation complete.  Elapsed: %.1f s\n', elapsed);
if record_video, close(v_vid); end
%% ===== 13. FINAL TEMPERATURE FIELD ======================================
figure('Color', 'w', 'Position', [150, 150, 1100, 430]);
pcolor(Z_grid, R_grid, T_2d);
shading interp;  colorbar;  colormap(hot);
hold on;
plot(x_contour, r_contour, 'c--', 'LineWidth', 1.5);
rectangle('Position', [ins_x_start, r_throat, ins_length, ins_thick], ...
    'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '--');
hold off;
axis equal tight;
title(sprintf('Final Temperature Field  (t = %.2f s)', tf));
xlabel('Axial Position [m]');  ylabel('Radius [m]');
%% ===== 14. INNER / OUTER WALL TEMPERATURE PROFILE =======================
figure('Color', 'w');
plot(x_contour, T_2d(1,:), 'r-',  'LineWidth', 2, 'DisplayName', 'Inner wall (gas side)');
hold on;
plot(x_contour, T_2d(n,:), 'b--', 'LineWidth', 2, 'DisplayName', 'Outer wall (ambient side)');
xlabel('Axial Position [m]', 'FontSize', 12);
ylabel('Temperature [K]',    'FontSize', 12);
title(sprintf('Wall Temperature Profile at t = %.2f s', tf), 'FontSize', 13);
legend('Location', 'best');
grid on;  hold off;