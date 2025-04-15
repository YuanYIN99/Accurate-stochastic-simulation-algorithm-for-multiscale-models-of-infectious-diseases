% This script produces Figure 5: The resolution of harvesting the 
% within-host viral load using a constant time step ∆t = δ and a constant 
% probability step δΨ.

% BH parameters: -----------------------------------------------
BH_parms = struct;
BH_parms.mu = 10 ^ (-5);           % natural death rate
BH_parms.lambda = 10 ^ (-5);       % S introduction rate
% Coefficient linking the WH viral load to the individual disease
% transmission rate:
BH_parms.l = 10^(-9.8);

% WH IC: -----------------------------------------------
T0_infect = 10 ^ (8);          % healthy target cells
Tstar0_infect = 10 ^ (0);      % infectious cells
V0_infect = 10 ^ (6);          % viral particles
TIV0_infect = [T0_infect; Tstar0_infect; V0_infect]; % this are initial
                               % conditions for the INFECTIOUS agents

% WH parameters: -----------------------------------------------
WH_parms = struct;
WH_parms.p = 10 ^ (1);           % viral production rate (p per day)
WH_parms.c = 10 ^ (0);           % viral clearance rate (c per day)
WH_parms.k = 10 ^ (-7);          % infection coefficient of target cells 
                                 % (k per virus per day)
WH_parms.mu_c = 10 ^ (-1);       % mortality rate of cells (mu_c per day)
WH_parms.delta_c = 10 ^ (1);     % extra mortality rate of infectious cells 
                                 % (delta_c per day)
WH_parms.R_0_WH = 1.1;           % within-host basic reproduction number
WH_parms.Lambda_c = WH_parms.R_0_WH * ...
    (WH_parms.mu_c * (WH_parms.mu_c + WH_parms.delta_c) * WH_parms.c)...
    / (WH_parms.k * WH_parms.p);

% Simulation setup: -----------------------------------------------
t0 = 0; 
t_endWH = 8;       
delta_ts = [0.0000937, 0.28, 1.0]; % different resolutions

colors = {'magenta', 'black', 'red'};
linestyles = {'-', '-.', ':'};
markers = {'None', 'd', 'x'};
marker_sizes = [15, 15, 25];
f1 = figure(); 
set(gcf, 'Position', get(0, 'Screensize')); % full screen the figure
set(groot,'defaultAxesFontName','Verdana');
set(groot,'defaultAxesFontSize',34);
set(0, 'DefaultLineLineWidth', 4);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
tlt1 = tiledlayout(1, 2);
title(tlt1, 'Coarse-graining of the within-host viral load', 'interpreter', 'latex', 'FontSize', 42)

for j = 1 : length(delta_ts)
    delta_t = delta_ts(j);
    color = colors{j};
    linestyle = linestyles{j};
    marker_size = marker_sizes(j);
    marker = markers{j};
    [t_WH, V_infect, Psi_table, const_Psi, t_table] = ...
        WH_exact(TIV0_infect, WH_parms, t_endWH, delta_t, BH_parms, color, linestyle, marker, marker_size);
end


function dydt = WH_infect_odesystem_Log(t, y, WH_parms)
% The ODE equations describing the dynamics of WH subsystem among I populations. 
% Note that:
% y(1) := T
% y(2) := Tstar
% y(3) := V.
dydt = [
    WH_parms.Lambda_c*exp(-y(1)) - WH_parms.k*exp(y(3)) - WH_parms.mu_c;
    WH_parms.k*exp(y(1)+y(3)-y(2)) - WH_parms.mu_c - WH_parms.delta_c;
    -WH_parms.c + WH_parms.p*exp(y(2)-y(3))
];
end


function [t_WH, V_infect, Psi_table, const_Psi, t_table] = ...
    WH_exact(TIV0_infect, WH_parms, t_endWH, delta_t, BH_parms, color, linestyle, marker, marker_size)
% This function solves the TIV model for the infectious population and
% outputs WH infomation depending on the age of infection from 0 to 't_endWH'.

Initial_Conds_infect = log(TIV0_infect); % log trandsform the ICs to solve ODE
% time points based on uniform t discretisation for extracting WH info:
t_WH = (delta_t/2 : delta_t : ceil(2*t_endWH/delta_t)*delta_t/2);
t_WH = [0, t_WH]; % include ICs to solve the WH ODE.

opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, y1] = ode15s(@(t,y) WH_infect_odesystem_Log(t, y, WH_parms), t_WH, ...
    Initial_Conds_infect, opts);
t_WH = t_WH(2 : end);
V_infect_Log = y1(: , 3); 
V_infect = exp(V_infect_Log);
V_infect = V_infect(2 : end); % WH viral load at t_WH

% Record the 'Psi' table based on uniform t discretisation:
Psi_table = cumtrapz(t_WH, V_infect) * BH_parms.l;

% Record 'Psi' with constant probability discretisation but variable t
% values:
Psi_table_min = min(Psi_table);
Psi_table_max = max(Psi_table);
const_Psi = linspace(Psi_table_min, Psi_table_max, length(Psi_table));
delta_Psi = const_Psi(2) - const_Psi(1);
% obtatain variable time points at which WH info is extracted:
t_table = interp1(Psi_table, t_WH, const_Psi, 'linear', 'extrap');
t_table = [0, t_table];
[~, y2] = ode15s(@(t,y) WH_infect_odesystem_Log(t, y, WH_parms), t_table, ...
    Initial_Conds_infect, opts);
t_table = t_table(2 : end);
V_infect_Log_ = y2(: , 3); 
V_infect_ = exp(V_infect_Log_);
V_infect_ = V_infect_(2 : end); % WH viral load at t_table

nexttile(1)
% plot V versus constant time discretisation:  
plot(t_WH, V_infect, 'LineStyle', linestyle, 'Marker', marker, 'Color', ...
    color, 'DisplayName', ['$\Delta t = \delta = $ ', num2str(delta_t), ' day'], ...
    'MarkerSize', marker_size)
ylabel('number of viral particles $V$')
xlabel('time since infection $\delta t$ (days)')
legend('Location', 'northeast');
title('Constant time step')
hold on
nexttile(2)
% plot V versus varying time discretisation (based on constant Psi):
plot(t_table, V_infect_, 'LineStyle', linestyle, 'Marker', marker, 'Color',...
    color, 'DisplayName', ['$\delta \Psi = $ ', num2str(delta_Psi)], ...
    'MarkerSize', marker_size)
xlabel('time since infection $\delta t$ (days)')
legend('Location', 'northeast');
title('Constant probability step')
hold on

end

