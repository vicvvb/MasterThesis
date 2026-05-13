%% =========================================================
%  E-COMPRESSOR DYNAMIC TORQUE & POWER ANALYSIS
%  ---------------------------------------------------------
%  Section 4.5 — Power and Torque requirements
%
%  Spool-up from 0 → N_max = 150 000 rpm in t_spool = 0.3 s.
%  Speed profile: s-curve (sinusoidal / raised-cosine ramp).
%  Aerodynamic model: Euler turbomachinery equation (Eq. 4.6.7–4.6.9).
%  Two inertia cases (Section 4.5.2.5):
%    Case 1 — Impeller only      (I_total = I_impeller)
%    Case 2 — Impeller + 1×extra (I_total = 2 × I_impeller)
% =========================================================
 
clear; clc; close all;
 
%% ── 1. SYSTEM PARAMETERS ─────────────────────────────────
%  All values taken from Tables 1 and 2 of the thesis.
 
% Impeller geometry
r2         = 31.58e-3;           % [m]      impeller outlet radius  (Table 2)
sigma_s    = 0.8;                % [-]      slip factor             (Table 2)
 
% Inertia  (Section 4.5.2.2, Eq. 4.6.4)
I_impeller = 3.5275e-5;          % [kg·m²]  impeller inertia from CAD model
 
% Design-speed targets
N_max      = 150000;             % [rpm]    peak rotational speed   (Table 2)
omega_max  = N_max * (2*pi/60);  % [rad/s]  = 15 708 rad/s
 
% Spool-up specification
t_spool    = 0.3;                % [s]      0 → N_max  (chosen value)
 
% Steady-state aerodynamic shaft power (corrected for bearing losses)
%   The shaft power the PMSM must supply includes bearing friction.
%   Bearing efficiency η_bearing = 0.99  (Section 4.5.2.3, Eq. 4.6.10).
%   P_aero,ss below is the PMSM shaft power = P_aero / η_bearing.
eta_bearing = 0.99;              % [-]      bearing efficiency
P_aero_fluid = 42084;           % [W]      aerodynamic power delivered to fluid
P_aero_ss    = P_aero_fluid / eta_bearing;  % [W]  total shaft power incl. bearing losses
% Note: P_aero_ss ≈ 42 509 W.
 
% Known steady-state mass flow rate  (Section 4.5.2.3, Eq. 4.6.8)
m_dot_ss   = 0.2535;             % [kg/s]
 
% Time vector
dt   = 1e-5;                     % [s]   time step
t    = 0 : dt : t_spool + 0.05;  % [s]   run a little beyond spool-up
 
%% ── 2. SPEED PROFILE — S-CURVE  (Section 4.5.2.1, Eq. 4.6.5) ────────────
%
%  ω(t) = (ω_max / 2) · [1 − cos(π·t / t_spool)]   for 0 ≤ t ≤ t_spool
%  ω(t) = ω_max                                      for t  > t_spool
%
%  This is the standard sinusoidal / raised-cosine ramp used in industrial
%  servo drives; it guarantees zero acceleration at both endpoints,
%  minimising jerk and avoiding current spikes in the drive electronics.
 
omega          = zeros(size(t));
idx_ramp       = t <= t_spool;
omega(idx_ramp)  = (omega_max/2) .* (1 - cos(pi .* t(idx_ramp) ./ t_spool));
omega(~idx_ramp) = omega_max;
 
% Angular acceleration  α(t) = dω/dt  (Eq. 4.6.5)
%   Analytical expression:  α(t) = (π·ω_max)/(2·t_spool) · sin(π·t/t_spool)
%   Computed analytically inside the ramp window; forced to zero afterwards
%   to avoid numerical artefacts from the gradient function at the kink.
alpha              = zeros(size(t));
alpha(idx_ramp)    = (pi * omega_max) / (2 * t_spool) .* sin(pi .* t(idx_ramp) ./ t_spool);
% α = 0 for t > t_spool  (already initialised to zero)
 
N_rpm = omega * (60 / (2*pi));   % [rpm]  for plotting
 
%% ── 3. AERODYNAMIC TORQUE & POWER  (Section 4.5.2.3, Eq. 4.6.7–4.6.11) ──
%
%  Starting point: Euler turbomachinery equation.
%  For a radial impeller with no inlet pre-swirl (C_u1 = 0), the Euler
%  specific work is (Eq. 4.6.7):
%
%    w_Euler(ω) = σ_s · U₂(ω)² = σ_s · (ω · r₂)²          [J/kg]
%
%  Mass flow scales linearly with tip speed — constant flow coefficient
%  Φ = C_m/U₂  (Eq. 4.6.8):
%
%    ṁ(ω) = ṁ_ss · (ω / ω_max)                             [kg/s]
%
%  Aerodynamic shaft power  P = ṁ · w_Euler  (Eq. 4.6.9):
%
%    P_aero(ω) = ṁ_ss · σ_s · r₂² · ω³ / ω_max             [W]
%
%  Corrected for bearing losses  (Eq. 4.6.10):
%
%    P_motor,aero(ω) = P_aero(ω) / η_bearing
%
%  Aerodynamic torque  T = P / ω  (Eq. 4.6.11):
%
%    T_aero(ω) = ṁ_ss · σ_s · r₂² · ω² / ω_max   ∝  ω²     [N·m]
%
%  Note: the ω³ power / ω² torque results look like the fan-affinity laws
%  but are derived from Euler mechanics — valid for a compressible fluid
%  (PR = 3.2) where the affinity laws would not apply.
 
% Euler specific work cross-check
w_euler_ss_kinematic = sigma_s * (omega_max * r2)^2;  % from velocity triangles
w_euler_ss_energy    = P_aero_fluid / m_dot_ss;       % from energy balance P/m_dot
 
fprintf('=== Euler model cross-check ===\n');
fprintf('  m_dot_ss                             : %.4f  kg/s\n',  m_dot_ss);
fprintf('  w_Euler,ss  (velocity triangles)     : %.1f  J/kg\n', w_euler_ss_kinematic);
fprintf('  w_Euler,ss  (energy balance P/m_dot) : %.1f  J/kg\n', w_euler_ss_energy);
fprintf('  Discrepancy                          : %.1f %%\n',    ...
        100*abs(w_euler_ss_kinematic - w_euler_ss_energy)/w_euler_ss_energy);
fprintf('  P_aero (fluid)                       : %.1f  W\n',    P_aero_fluid);
fprintf('  P_aero,ss  (shaft incl. bearings)    : %.1f  W\n\n',  P_aero_ss);
 
% Speed-dependent aerodynamic quantities
m_dot   = m_dot_ss .* (omega ./ omega_max);                     % [kg/s]  Eq. 4.6.8
w_euler = sigma_s .* (omega .* r2).^2;                          % [J/kg]  Eq. 4.6.7
P_aero  = (m_dot .* w_euler) ./ eta_bearing;                    % [W]     Eq. 4.6.9+4.6.10
 
T_aero  = zeros(size(omega));                                    % [N·m]   Eq. 4.6.11
T_aero(omega > 0) = P_aero(omega > 0) ./ omega(omega > 0);
 
% Verify boundary condition
fprintf('=== Verification at ω_max ===\n');
fprintf('  P_aero,ss  (input)      : %.1f  W\n', P_aero_ss);
fprintf('  P_aero(ω_max) (calc)    : %.1f  W\n', P_aero(end));
fprintf('  T_aero,ss  (calc)       : %.2f  N·m\n\n', T_aero(end));
 
%% ── 4. INERTIA TORQUE — TWO CASES  (Section 4.5.2.5, Table 1) ───────────
%
%  Case 1: Impeller only          I_total = I_impeller
%  Case 2: Impeller + 1× extra    I_total = 2 × I_impeller
%           (rotor inertia ≈ I_impeller — see sensitivity study)
%
%  T_inertia(t) = I_total · α(t)   (Eq. 4.6.6)
 
I_cases = [1, 2] * I_impeller;           % [kg·m²]  total inertia per case
case_labels = {"Impeller only", ...
               "Impeller + 1\times extra rotor"};
colors  = ["#0072BD", "#D95319"];        % blue, orange — matches thesis figures
 
T_inertia = cell(1, 2);
T_total   = cell(1, 2);
P_total   = cell(1, 2);
 
fprintf('=== Peak values per inertia case ===\n');
for k = 1:2
    T_inertia{k} = I_cases(k) .* alpha;                % [N·m]  Eq. 4.6.6
    T_total{k}   = T_inertia{k} + T_aero;              % [N·m]  Eq. 4.6.12
    P_total{k}   = T_total{k}  .* omega;               % [W]    Eq. 4.6.13
 
    [pk_Ti, ~]  = max(T_inertia{k});
    [pk_Tt, ~]  = max(T_total{k});
    [pk_P,  ~]  = max(P_total{k});
 
    fprintf('  %s  (I = %.4e kg·m²)\n', case_labels{k}, I_cases(k));
    fprintf('    Peak inertia torque  : %6.2f  N·m\n', pk_Ti);
    fprintf('    Peak total torque    : %6.2f  N·m\n', pk_Tt);
    fprintf('    Peak shaft power     : %6.1f  W  (%5.1f kW)\n\n', pk_P, pk_P/1e3);
end
 
%% ── 5. PLOTS ─────────────────────────────────────────────
%  Six separate figures matching Section 4.5 of the thesis:
%    Fig 16 — Rotor speed
%    Fig 17 — Angular acceleration
%    Fig 18 — Inertia torque (two cases)
%    Fig 19 — (System sketch — not reproduced here)
%    Fig 20 — Aerodynamic (flow) torque
%    Fig 21 — Total required torque (two cases)
%    Fig 22 — Total required shaft power (two cases)
 
t_ms   = t * 1e3;                      % [ms]  time axis for all plots
xlims  = [0, (t_spool + 0.03)*1e3];   % shared x-range
FW = 560;  FH = 420;                   % figure size [px]
fsAxis  = 12;  fsTitle = 13;           % font sizes (fsTitle kept for compatibility)
fsAnnot = 12;                          % annotation font size (increased)
LW      = 2.0;                         % line width
 
% ── Fig 16 · Rotor speed ──────────────────────────────────
figure('Name','Fig 16 – Rotor speed', 'Position',[50 620 FW FH], 'Color','w');
plot(t_ms, N_rpm/1e3, 'k', 'LineWidth', LW);
hold on;
xline(t_spool*1e3, '--r', 't_{spool}', ...
      'LabelVerticalAlignment','bottom', 'FontSize', fsAxis-1);
% Annotation at t_spool
[~, idx_sp] = min(abs(t - t_spool));
plot(t_ms(idx_sp), N_rpm(idx_sp)/1e3, 'ok', 'MarkerFaceColor','k', 'MarkerSize',7);
text(t_ms(idx_sp), N_rpm(idx_sp)/1e3, ...
     sprintf('%.0f krpm  ', N_rpm(idx_sp)/1e3), ...
     'FontSize', fsAnnot, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
xlabel('Time  (ms)',    'FontSize', fsAxis);
ylabel('Speed  (krpm)', 'FontSize', fsAxis);
grid on;  box on;  set(gca,'FontSize', fsAxis);  xlim(xlims);
ylim([0, N_max/1e3 * 1.20]);   % 20% headroom above max speed
 
% ── Fig 17 · Angular acceleration ─────────────────────────
figure('Name','Fig 17 – Angular acceleration', 'Position',[630 620 FW FH], 'Color','w');
plot(t_ms, alpha/1e3, 'k', 'LineWidth', LW);
hold on;
xline(t_spool*1e3, '--r', 't_{spool}', ...
      'LabelVerticalAlignment','bottom', 'FontSize', fsAxis-1);
% Peak annotation
[pk_alpha, idx_alpha] = max(alpha);
plot(t_ms(idx_alpha), pk_alpha/1e3, 'ok', 'MarkerFaceColor','k', 'MarkerSize',7);
text(t_ms(idx_alpha), pk_alpha/1e3, ...
     sprintf('  %.1f krad/s²', pk_alpha/1e3), ...
     'FontSize', fsAnnot, 'VerticalAlignment','bottom');
xlabel('Time  (ms)',         'FontSize', fsAxis);
ylabel('\alpha  (krad/s²)',  'FontSize', fsAxis);
grid on;  box on;  set(gca,'FontSize', fsAxis);  xlim(xlims);
 
% ── Fig 18 · Inertia torque ───────────────────────────────
figure('Name','Fig 18 – Inertia torque', 'Position',[50 150 FW FH], 'Color','w');
hold on;
pk_Ti_all = zeros(1,2);
for k = 1:2
    plot(t_ms, T_inertia{k}, 'Color', colors(k), ...
         'LineWidth', LW, 'DisplayName', case_labels{k});
    [pk_Ti, idx_Ti] = max(T_inertia{k});
    pk_Ti_all(k) = pk_Ti;
    plot(t_ms(idx_Ti), pk_Ti, 'o', 'MarkerSize',7, ...
         'MarkerFaceColor', colors(k), 'MarkerEdgeColor', colors(k), ...
         'HandleVisibility','off');
    % Place label above the peak marker, inside the expanded y-axis
    text(t_ms(idx_Ti), pk_Ti, sprintf('  %.2f N·m', pk_Ti), ...
         'FontSize', fsAnnot, 'Color', colors(k), 'VerticalAlignment','bottom');
end
xline(t_spool*1e3, '--r', 't_{spool}', ...
      'LabelVerticalAlignment','bottom', 'FontSize', fsAxis-1);
xlabel('Time  (ms)',         'FontSize', fsAxis);
ylabel('T_{inertia}  (N·m)', 'FontSize', fsAxis);
grid on;  box on;  set(gca,'FontSize', fsAxis);  xlim(xlims);
ylim([0, max(pk_Ti_all) * 1.20]);   % 20% headroom so labels are never clipped
lg = legend('Location','northwest','FontSize',9,'Box','on','NumColumns',1,'Interpreter','tex');
lg.ItemTokenSize = [15, 8];
 
% ── Fig 20 · Aerodynamic torque ───────────────────────────
figure('Name','Fig 20 – Aerodynamic torque', 'Position',[630 150 FW FH], 'Color','w');
plot(t_ms, T_aero, 'Color','#7E2F8E', 'LineWidth', LW);
hold on;
xline(t_spool*1e3, '--r', 't_{spool}', ...
      'LabelVerticalAlignment','bottom', 'FontSize', fsAxis-1);
% Annotate steady-state value — placed to the LEFT of the marker to avoid right-edge clipping
plot(t_ms(idx_sp), T_aero(idx_sp), 'o', 'MarkerSize',7, ...
     'MarkerFaceColor','#7E2F8E', 'MarkerEdgeColor','#7E2F8E');
text(t_ms(idx_sp), T_aero(idx_sp), ...
     sprintf('T_{aero,ss} = %.2f N·m  ', T_aero(idx_sp)), ...
     'FontSize', fsAnnot, 'Color','#7E2F8E', ...
     'VerticalAlignment','bottom', 'HorizontalAlignment','right');
xlabel('Time  (ms)',      'FontSize', fsAxis);
ylabel('T_{aero}  (N·m)', 'FontSize', fsAxis);
grid on;  box on;  set(gca,'FontSize', fsAxis);  xlim(xlims);
ylim([0, T_aero(end) * 1.20]);   % 20% headroom
 
% ── Fig 21 · Total required torque ────────────────────────
figure('Name','Fig 21 – Total required torque', 'Position',[1210 620 FW FH], 'Color','w');
hold on;
pk_Tt_all = zeros(1,2);
for k = 1:2
    plot(t_ms, T_total{k}, 'Color', colors(k), ...
         'LineWidth', LW, 'DisplayName', case_labels{k});
    [pk_Tt, idx_Tt] = max(T_total{k});
    pk_Tt_all(k) = pk_Tt;
    plot(t_ms(idx_Tt), pk_Tt, 'o', 'MarkerSize',7, ...
         'MarkerFaceColor', colors(k), 'MarkerEdgeColor', colors(k), ...
         'HandleVisibility','off');
    % Place label above the peak marker, inside the expanded y-axis
    text(t_ms(idx_Tt), pk_Tt, sprintf('  %.2f N·m', pk_Tt), ...
         'FontSize', fsAnnot, 'Color', colors(k), 'VerticalAlignment','bottom');
end
xline(t_spool*1e3, '--r', 't_{spool}', ...
      'LabelVerticalAlignment','bottom', 'FontSize', fsAxis-1);
xlabel('Time  (ms)',       'FontSize', fsAxis);
ylabel('T_{total}  (N·m)', 'FontSize', fsAxis);
grid on;  box on;  set(gca,'FontSize', fsAxis);  xlim(xlims);
ylim([0, max(pk_Tt_all) * 1.20]);   % 20% headroom
lg = legend('Location','northwest','FontSize',9,'Box','on','NumColumns',1,'Interpreter','tex');
lg.ItemTokenSize = [15, 8];
 
% ── Fig 22 · Total required shaft power ───────────────────
figure('Name','Fig 22 – Total required shaft power', 'Position',[1210 150 FW FH], 'Color','w');
hold on;
pk_P_all = zeros(1,2);
for k = 1:2
    plot(t_ms, P_total{k}/1e3, 'Color', colors(k), ...
         'LineWidth', LW, 'DisplayName', case_labels{k});
    [pk_P, idx_P] = max(P_total{k});
    pk_P_all(k) = pk_P/1e3;
    plot(t_ms(idx_P), pk_P/1e3, 'o', 'MarkerSize',7, ...
         'MarkerFaceColor', colors(k), 'MarkerEdgeColor', colors(k), ...
         'HandleVisibility','off');
    % Place label above the peak marker, inside the expanded y-axis
    text(t_ms(idx_P), pk_P/1e3, sprintf('  %.1f kW', pk_P/1e3), ...
         'FontSize', fsAnnot, 'Color', colors(k), 'VerticalAlignment','bottom');
end
yline(P_aero_fluid/1e3, ':k', sprintf('P_{ss} (%.1f kW)', P_aero_fluid/1e3), ...
      'LabelHorizontalAlignment','left', 'FontSize', fsAxis-1);
xline(t_spool*1e3, '--r', 't_{spool}', ...
      'LabelVerticalAlignment','bottom', 'FontSize', fsAxis-1);
xlabel('Time  (ms)', 'FontSize', fsAxis);
ylabel('Power  (kW)', 'FontSize', fsAxis);
grid on;  box on;  set(gca,'FontSize', fsAxis);  xlim(xlims);
ylim([0, max(pk_P_all) * 1.20]);   % 20% headroom
lg = legend('Location','northwest','FontSize',9,'Box','on','NumColumns',1,'Interpreter','tex');
lg.ItemTokenSize = [15, 8];
 
fprintf('Done. Figures 16–18, 20–22 generated.\n');