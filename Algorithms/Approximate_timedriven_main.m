function Approximate_timedriven_main(k)
% This function implements the approximate time-driven algorithm based on
% the kth initial condition.

% To be submitted to the Spartan HPC in the University of Melbourne:
%pc = parcluster('local');
%pc.JobStorageLocation = getenv('SCRATCH');
%pc.NumWorkers = 1; % 30 parallel workers

% To run on a local computer with a 14-core CPU: 
pool = gcp('nocreate'); % 'nocreate' prevents creating a new pool if one doesn't exist
if ~isempty(pool)
    delete(pool); % Delete the existing pool
end
parpool(14);

% BH IC: -----------------------------------------------
S0_values = [100, 200, 400, 800];
I0_values = [10, 20, 40, 80];
S0 = S0_values(k); I0 = I0_values(k);
SI0 = [S0; I0];  

% BH parameters: -----------------------------------------------
BH_parms = struct;
BH_parms.mu = 10 ^ (-5);           % natural death rate
BH_parms.lambda = 10 ^ (-5);       % S introduction rate
% Coefficient linking the WH viral load to the individual disease
% transmission rate:
l_values = [10^(-9.8), 10^(-9.8), 10^(-9.8), 10^(-9.8)]; 
BH_parms.l = l_values(k);

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
t_end_values = [20, 20, 20, 20];
t_endBH = t_end_values(k);
t_endWH = t_endBH; 

% Resolution to extract the poplation-level SI dynamics over 800
% repetitions:
rep = 800; 
delta_t_uniform = 0.1;
t_uniform = (t0 : delta_t_uniform : t_endBH);
   
% The fixed time steps Delta t: --------------------------------
pt_N = 130;
delta_t_values = [0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N)];

% Resolution to extract WH info: ------------------------------------------
delta_t = delta_t_values(k, :);     

% Store the setup parameters and the IC: --------------------------------
setup_filename = ['timedriven_multi_main_Setup_S', num2str(S0), 'I', ...
    num2str(I0), '_tend', num2str(t_endBH), '_rep', num2str(rep), '_l', ...
    num2str(BH_parms.l), '_k', num2str(k), '.mat'];
save(setup_filename, "BH_parms", "WH_parms", "SI0", "TIV0_infect", "t_endBH", ...
    "t0", "t_uniform", "t_endWH", "delta_t", "rep", 'k');


% RUN the approximate time-driven algorithm: ------------------------------
Approx_Tdriven_Runs_deltat(BH_parms, SI0, WH_parms, TIV0_infect, ...
    delta_t, t0, t_endBH, t_uniform, rep, k)

end



function Approx_Tdriven_Runs_deltat(BH_parms, SI0, WH_parms, TIV0_infect, ...
    delta_t, t0, t_endBH, t_uniform, rep, k)
% This function runs the approximate time-driven algorithm for 'rep' times for 
% different delta t grid sizes.

% Initialisation: to store the total run times for all repetitions given
% different delta t's.
times_tabWH_approx = zeros(length(delta_t), 1); % run time to tabulate the WH info
runs_elapsed_time_approx = zeros(length(delta_t), 1); % run time to run the algorithm

% To store the overall info extrcated at 't_uniform' timepoints:
stocha_info_S_approx_deltat = zeros(rep, length(t_uniform), length(delta_t)); 
stocha_info_I_approx_deltat = zeros(rep, length(t_uniform), length(delta_t));

% RUN:
for i = 1 : length(delta_t)
    delta_t_i = delta_t(i);
    % for the approximate time-driven algorithm -- tabulate viral load info 
    % at [delta_t/2, 3delta_t/2, ...]:
    time_WH = tic; 
    [t_approxTdriven_WH, V_approxTdriven] = WH_approxTdriven(delta_t_i, ...
        t_endBH, TIV0_infect, WH_parms);
    times_tabWH_approx(i) = toc(time_WH);

    % run rep times:
    runs_time_approx = tic; 
    [stocha_info_S_approxTdriven, stocha_info_I_approxTdriven] = ...
        Approx_Tdriven_Runs(delta_t_i, rep, BH_parms, t0, t_endBH, SI0, ...
        t_approxTdriven_WH, V_approxTdriven, t_uniform);
    runs_elapsed_time_approx(i) = toc(runs_time_approx);
    stocha_info_S_approx_deltat(:, :, i) = stocha_info_S_approxTdriven;
    stocha_info_I_approx_deltat(:, :, i) = stocha_info_I_approxTdriven;
end

% Store info into relevant .mat files:
ApprxTimeDriven_filename = ['Approx_RunsDeltats_S', num2str(SI0(1)), 'I', ...
    num2str(SI0(2)), '_tendBH', num2str(t_endBH), '_rep', num2str(rep), '_l', ...
    num2str(BH_parms.l), '_k', num2str(k), '.mat'];
save(ApprxTimeDriven_filename, 't_uniform', 'stocha_info_S_approx_deltat', ...
    'stocha_info_I_approx_deltat', 'times_tabWH_approx', ...
    'runs_elapsed_time_approx', '-v7.3');


% One can also choose a specific delta_t resolution 'idx' to visualise the
% SI dynamics:
%figure;
%hold on;
%for r = 1:rep
%    plot(t_uniform, squeeze(stocha_info_S_approxTdriven(r, :, idx)), 'b');
%    plot(t_uniform, squeeze(stocha_info_I_approxTdriven(r, :, idx)), 'r');
%end

end


function [stocha_info_S_approxTdriven, stocha_info_I_approxTdriven] = ...
    Approx_Tdriven_Runs(delta_t_i, rep, BH_parms, t0, t_endBH, SI0, ...
    t_approxTdriven_WH, V_approxTdriven, t_uniform) 
% This function runs the approximate time-driven for 'rep' times in parallel:

% Initialisation: to store SI info at 't_uniform':
stocha_info_S_approxTdriven = zeros(rep, length(t_uniform)); 
stocha_info_I_approxTdriven = zeros(rep, length(t_uniform));

parfor run = 1 : rep
    [S_approxTdriven_uniform, I_approxTdriven_uniform] = ...
        Approx_Tdriven_OneRun(BH_parms, t0, delta_t_i, t_endBH, SI0, t_approxTdriven_WH, ...
        V_approxTdriven, t_uniform);
    % store relevant info:
    stocha_info_S_approxTdriven(run, :) = S_approxTdriven_uniform;
    stocha_info_I_approxTdriven(run, :) = I_approxTdriven_uniform;
end

end


function [S_approxTdriven_uniform, I_approxTdriven_uniform] = ...
    Approx_Tdriven_OneRun(BH_parms, t0, delta_t, t_endBH, SI0, t_approxTdriven_WH, ...
    V_approxTdriven, t_uniform)
% This function runs the approximate time-driven algorithm once.

% Initilise the system:
[t_approxTdriven, S_approxTdriven, I_approxTdriven, t_index, k] = ...
    init_system_approxTdriven(t0, delta_t, t_endBH, SI0, t_approxTdriven_WH);

while t0+(t_index-1)*delta_t <= t_endBH
    % Calculate propensities:
    alpha_b = BH_parms.lambda * delta_t;
    alpha_ds = BH_parms.mu * delta_t * S_approxTdriven(t_index-1);
    alpha_di = BH_parms.mu * delta_t * I_approxTdriven(t_index-1);
    alpha_i = delta_t * S_approxTdriven(t_index-1) * BH_parms.l * ...
        sum(V_approxTdriven(k(1:I_approxTdriven(t_index-1))));
    alpha_total = alpha_b + alpha_ds + alpha_di + alpha_i;

    % If there is any event to happen:
    if rand <= alpha_total
        % Decide which event to happen:
        u = rand * alpha_total;

        if u <= alpha_b % introduction of S into the system
            S_approxTdriven(t_index) = S_approxTdriven(t_index-1) + 1;
            I_approxTdriven(t_index) = I_approxTdriven(t_index-1);
            k = [k; 0];
            k(1:I_approxTdriven(t_index)) = k(1:I_approxTdriven(t_index)) + 1;

        elseif (u > alpha_b) && (u <= alpha_b+alpha_ds) % S dies
            S_approxTdriven(t_index) = S_approxTdriven(t_index-1) - 1;
            I_approxTdriven(t_index) = I_approxTdriven(t_index-1);
            k = Delete_ele(length(k), k);
            k(1:I_approxTdriven(t_index)) = k(1:I_approxTdriven(t_index)) + 1;

        elseif (u > alpha_b+alpha_ds) && (u <= alpha_b+alpha_ds+alpha_di) % Idies
            S_approxTdriven(t_index) = S_approxTdriven(t_index-1);
            I_approxTdriven(t_index) = I_approxTdriven(t_index-1) - 1;
            % randomly select a infectious person to remove:
            I_remove_index = randi(I_approxTdriven(t_index-1));
            k = Delete_ele(I_remove_index, k);
            k(1:I_approxTdriven(t_index)) = k(1:I_approxTdriven(t_index)) + 1;

        else % infection happens
            S_approxTdriven(t_index) = S_approxTdriven(t_index-1) - 1;
            I_approxTdriven(t_index) = I_approxTdriven(t_index-1) + 1;
            k(1:I_approxTdriven(t_index-1)) = k(1:I_approxTdriven(t_index-1)) + 1;
            k(I_approxTdriven(t_index)) = 1;
        end

    else % no event happens
        S_approxTdriven(t_index) = S_approxTdriven(t_index-1);
        I_approxTdriven(t_index) = I_approxTdriven(t_index-1);
        k(1:I_approxTdriven(t_index)) = k(1:I_approxTdriven(t_index)) + 1;
    end
    
    % advance in time
    t_index = t_index + 1;
end

% Harvest the SI info at 't_uniform' post-simulation:
S_approxTdriven_uniform = zeros(size(t_uniform))';
I_approxTdriven_uniform = zeros(size(t_uniform))';
for m = 1 : length(t_uniform)
    t_m = t_uniform(m);
    t_stocha_index = max(find(t_approxTdriven <= t_m)); 
    % the SI dyanmics at t = t_i are the same as the ones at t_stocha_index;
    S_approxTdriven_uniform(m) = S_approxTdriven(t_stocha_index);
    I_approxTdriven_uniform(m) = I_approxTdriven(t_stocha_index);
end

end


function [t_approxTdriven, S_approxTdriven, I_approxTdriven, t_index, k] = ...
    init_system_approxTdriven(t0, delta_t, t_endBH, SI0, t_approxTdriven_WH)
% This function initialises the system for approximate time-driven
% algorithm.

t_approxTdriven = (t0 : delta_t : t_endBH);
t_index = 1;
S_approxTdriven = zeros(size(t_approxTdriven)); S_approxTdriven(t_index) = SI0(1);
I_approxTdriven = zeros(size(t_approxTdriven)); I_approxTdriven(t_index) = SI0(2);

% If the age of infection for all infectious ppl at t = 0 is 0:
% index array for looking up V at those discrete time steps.
k = zeros(sum(SI0), 1);
for i = 1 : SI0(2)
    age_infect_i = t_approxTdriven(t_index) - 0; 
    [~, index] = min(abs(age_infect_i-t_approxTdriven_WH));
    k(i) = index;
end

t_index = t_index + 1;

end


function [t_approxTdriven_WH, V_approxTdriven] = WH_approxTdriven(delta_t, ...
    t_endBH, TIV0_infect, WH_parms)
% This function solved the WH viral load numerically, and stores 'V' into a
% table for later look up. 

% time points for extracting WH info:
t_approxTdriven_WH = (delta_t/2 : delta_t : ceil(2*t_endBH/delta_t)*delta_t/2);
t_approxTdriven_WH = [0, t_approxTdriven_WH]; % include ICs to solve the WH ODE.

opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
Initial_Conds_infect = log(TIV0_infect); % log trandsform to solve ODE
[~, y1] = ode15s(@(t,y) WH_infect_odesystem_Log(t, y, WH_parms), ...
    t_approxTdriven_WH, Initial_Conds_infect, opts);
V_infect_Log = y1(: , 3); 
V_approxTdriven = exp(V_infect_Log);

V_approxTdriven = V_approxTdriven(2 : end);
t_approxTdriven_WH = t_approxTdriven_WH(2 : end);

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


function matrix = Delete_ele(index, matrix)
% Remove the 'index_th' row in the matrix (note that if "matrix" is an array, 
% we need have a column vector!!).
if size(matrix, 1) == 1
    matrix = [];
else
    if index > 1 && index < size(matrix, 1)
        matrix = [matrix(1:index-1, :); matrix(index+1:size(matrix, 1), :)];
    else
        if index == 1
            matrix = matrix(index+1:size(matrix, 1), :);
        elseif index == size(matrix, 1)
            matrix = matrix(1:size(matrix, 1)-1, :);
        end
    end
end
end
