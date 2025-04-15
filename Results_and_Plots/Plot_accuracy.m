function result_main()
% This function extracts the accuracy information of the novel exact
% algorithm and the existing approximate time-driven algorithm, with respect
% to a 'golden standard' algorithm and plots out the results. 

% Extract summary statistics info and stored that into 'info_filename.mat':
info_filename = extract_summary_stats();
% Load the information:
load(info_filename, "delta_t", "approx_error_mean", "exact_error_mean");

% Define a critical threshold 'err_rel_threshold' for the relative error. 
% For Delta_t's larger than the Delta_t_{critical} corresponding to 'err_rel_threshold', 
% the algorithms are said to be inaccurate. 
err_rel_threshold = 0.05;
[exact_crit_stepsize, exact_crit_err, apprx_crit_stepsize, apprx_crit_err] ...
    = find_crit_info(approx_error_mean, exact_error_mean, err_rel_threshold, delta_t);

% Accuracy plot:
accuracy_plot(delta_t, approx_error_mean, exact_error_mean, ...
    exact_crit_stepsize, exact_crit_err, apprx_crit_stepsize, apprx_crit_err)

end



function info_filename = extract_summary_stats()
% This function extracts the accuracy info of two algorithms:
% 1). novel exact algorithm;
% 2). existing approximate time-driven algorithm;
% with respect to a 'golden standard' algorithm.
% It stores it into 'info_filename.mat' files for future plotting.


% Set up files which stores the parameters & ICs used in the simulations:
setup = {
    'timedriven_multi_main_Setup_S100I10_tend20_rep800_l1.5849e-10_k1.mat';
    'timedriven_multi_main_Setup_S200I20_tend20_rep800_l1.5849e-10_k2.mat';
    'timedriven_multi_main_Setup_S400I40_tend20_rep800_l1.5849e-10_k3.mat';
    'timedriven_multi_main_Setup_S800I80_tend20_rep800_l1.5849e-10_k4.mat'};

% The different Delta t values used in all simulations:
delta_t = 0.0000937 * 1.084.^(0:130); 
rep = 800; 

% Initialisation: to store the accuracy summary statistics:
% For the novel exact algorithm:
Exact_SIapprx_T = zeros(length(setup), length(delta_t), rep); 
exact_error_mean = zeros(length(setup), length(delta_t));
% For the existing approximate time-driven algorithm:
ApprxTimeDriven_SIapprx_T = zeros(length(setup), length(delta_t), rep); 
approx_error_mean = zeros(length(setup), length(delta_t));
% For the 'Golden Standard' algorithm:
GS_T = zeros(length(setup), 1); 

% Initialisation: to store the average run time for two algorithms of
% interest:
%exact_iterTime_mean = zeros(length(setup), length(delta_t));
%approx_iterTime_mean = zeros(length(setup), length(delta_t));

for k = 1 : length(setup) % loop through different simulation set ups
    setup_filename = setup{k};
    load(setup_filename);
    
    % Load simulation outputs for the exact algorithm and the approxiamte 
    % time-driven algorithm:
    ExactMultiscale_filename = ['ExactMulti_RunsDeltats_S', num2str(SI0(1)), 'I', ...
        num2str(SI0(2)), '_tendBH', num2str(t_endBH), '_rep', num2str(rep), '_l', ...
        num2str(BH_parms.l), '_k', num2str(k), '.mat'];
    load(ExactMultiscale_filename, 'stocha_info_S_effiMulti_deltat', ...
        'stocha_info_I_effiMulti_deltat', 'runs_elapsed_time_effiMulti');
    ApprxTimeDriven_filename = ['Approx_RunsDeltats_S', num2str(SI0(1)), 'I', ...
        num2str(SI0(2)), '_tendBH', num2str(t_endBH), '_rep', num2str(rep), '_l', ...
        num2str(BH_parms.l), '_k', num2str(k), '.mat'];
    load(ApprxTimeDriven_filename, 'stocha_info_S_approx_deltat', ...
        'stocha_info_I_approx_deltat', 'runs_elapsed_time_approx');

    % ACCURACY: -----------------------------------------------------------
    % For exact algorithm, find the first time point where S~I:
    for j = 1 : length(delta_t)
        for i = 1 : rep
            S_exact = stocha_info_S_effiMulti_deltat(i, 1:end, j);
            I_exact = stocha_info_I_effiMulti_deltat(i, 1:end, j);
            % Find the index of the first entry where S~I:
            [~, index_exact] = min(abs(S_exact - I_exact));
            Exact_SIapprx_T(k, j, i) = t_uniform(index_exact);    
        end
    end
    % For approximate time-driven algorithm, find the first time point 
    % where S～I:
    for p = 1 : length(delta_t)
        for q = 1 : rep
            S_apprx = stocha_info_S_approx_deltat(q, 1:end, p); 
            I_apprx = stocha_info_I_approx_deltat(q, 1:end, p);
            % Find the index of the first entry where S~I:
            [~, index_apprx] = min(abs(S_apprx - I_apprx));
            ApprxTimeDriven_SIapprx_T(k, p, q) = t_uniform(index_apprx);
        end
    end

    % Obtain the timepoint at which S~I for the 'golden-standard' algorithm:
    GS_T(k) = mean(ApprxTimeDriven_SIapprx_T(k, 1, :));

    % Obtain the mean relative error:
    for m = 1 : length(delta_t)
        approx_error_mean(k, m) = mean((ApprxTimeDriven_SIapprx_T(k, m, :) - GS_T(k)) / GS_T(k));
        exact_error_mean(k, m) = mean((Exact_SIapprx_T(k, m, :) - GS_T(k)) / GS_T(k));
    end

end

info_filename = 'bothAlgorithms_criticalResolution_accuracy.mat';
save(info_filename, "delta_t", "rep", "approx_error_mean", "exact_error_mean");

end



function [exact_crit_stepsize, exact_crit_err, apprx_crit_stepsize, apprx_crit_err] ...
    = find_crit_info(approx_error_mean, exact_error_mean, err_rel_threshold, delta_t)
% This function finds the information at the critical Delta t value where two
% algorithms become inaccurate. 

experiment_num = 4; % simulations are performed over 4 setups

% Initialisation: to store the critical Delta_t's that gives
% 'err_rel_threshold' together with the error when compared to the 'golden-
% standard' algorithm at those critical Delta _t's.
exact_crit_stepsize = zeros(1, experiment_num);  
exact_crit_err = zeros(1, experiment_num); 
apprx_crit_stepsize = zeros(1, experiment_num); 
apprx_crit_err = zeros(1, experiment_num); 
for i = 1 : experiment_num
    apprx_index = find(abs(approx_error_mean(i, :))>err_rel_threshold, 1, 'first');
    exact_index = find(abs(exact_error_mean(i, :))>err_rel_threshold, 1, 'first');

    apprx_crit_stepsize(i) = delta_t(apprx_index);
    apprx_crit_err(i) = approx_error_mean(i, apprx_index);
    exact_crit_stepsize(i) = delta_t(exact_index);
    exact_crit_err(i) = exact_error_mean(i, exact_index);
end
end



function accuracy_plot(delta_t, approx_error_mean, exact_error_mean, ...
    exact_crit_stepsize, exact_crit_err, apprx_crit_stepsize, apprx_crit_err)
% This function plots 'Figure 7' in the manuscript. It visualises accuracy 
% results for the approximate time-driven algorithm and the novel exact 
% algorithm compared with the 'golden-standard' one, across four simulation
% setups with different initial populations. 


% The four simulation setups have S initial conditions as follows:
S_population = [100, 200, 400, 800];
% Format different line colours for the four setups: 
colors4 = {'#db6104', '#008744', '#0057e7', '#cf4bcd'};

% Initialisation: to format legends
exact_pl_legends = zeros(length(S_population), length(delta_t));
approx_pl_legends = zeros(length(S_population), length(delta_t));

% Figure format:
f1 = figure(); 
set(gcf, 'Position', get(0, 'Screensize')); % full screen the figure
set(groot,'defaultAxesFontName','Verdana');
set(groot,'defaultAxesFontSize',40);
set(0, 'DefaultLineLineWidth', 6);
set(groot,'defaultLineMarkerSize', 8);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
tlt1 = tiledlayout(1, 1);
nexttile 
% For formatting legends:
for k = 1 : length(S_population)
    exact_pl_legends(k, :) = semilogx(delta_t, exact_error_mean(k, :), ...
        '-', 'Color', colors4{k}, 'HandleVisibility', 'off'); hold on 
    approx_pl_legends(k, :) = semilogx(delta_t, approx_error_mean(k, :), ...
        ':', 'Color', colors4{k}, 'HandleVisibility', 'off'); hold on 
    xlim([min(delta_t), max(delta_t)])
    xlabel(['$\Delta t= \delta $', ' (days)'], 'interpreter', 'latex'); 
    ylabel('$(t_{S \approx I}^{\mathrm{Algo. 1, 2}} - t_{S \approx I}^{\mathrm{GS}})/t_{S \approx I}^{\mathrm{GS}}$', ...
        'interpreter', 'latex')
end

semilogx(ones(1, 2), ones(1, 2)*NaN, 'w', 'DisplayName', 'Algo. 1'); hold on 
for k = 1 : length(S_population)
    semilogx(ones(1, 2), ones(1, 2)*NaN, ':', 'Color', colors4{k}, 'DisplayName', ['S = ', num2str(S_population(k)), ' ppl.']); hold on 
end
semilogx(ones(1, 2), ones(1, 2)*NaN, 'w', 'DisplayName', 'Algo. 2'); hold on 
for k = 1 : length(S_population)
    semilogx(ones(1, 2), ones(1, 2)*NaN, '-', 'Color', colors4{k}, 'DisplayName', ['S = ', num2str(S_population(k)), ' ppl.']); hold on 
end
lgd1 = legend('Location', 'northwest');
lgd1.NumColumns = 2;

% Plot results: 
for k = 1 : length(S_population)
    plot(exact_crit_stepsize(k), exact_crit_err(k), '*', 'Color', colors4{k}, 'MarkerSize', 40, 'HandleVisibility', 'off'); hold on
    plot(apprx_crit_stepsize(k), apprx_crit_err(k), '.', 'Color', colors4{k}, 'MarkerSize', 100, 'HandleVisibility', 'off'); hold on
end
title('Accuracy comparison given different population sizes');

% Save output as a .jpg file:
%saveas(f1, 'Accuracy_mean_err_k1234.jpg');

end

