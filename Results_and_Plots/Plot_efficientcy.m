% This script compares the runtime of two algorithms, and produces Figure
% 8.

% Different ICs:
S_population = [100, 200, 400, 800];
I_population = [10, 20, 40, 80];
colors4 = {'#db6104', '#008744', '#0057e7', '#cf4bcd'};

% find the critical resolutions where the relative errors are bounded by 5%:
load("bothAlgorithms_criticalResolution_accuracy.mat")
err_rel_threshold = 0.05;
[exact_crit_stepsize, exact_crit_err, apprx_crit_stepsize, apprx_crit_err] = ...
    find_crit_info(approx_error_mean, exact_error_mean, err_rel_threshold);

delta_ts = 0.0000937 * 1.084.^(0:130); 
rep = 800; % we have a total of 800 repititions for each run

f1 = figure(); 
set(gcf, 'Position', get(0, 'Screensize')); % full screen the figure
set(groot,'defaultAxesFontName','Verdana');
set(groot,'defaultAxesFontSize',40);
set(0, 'DefaultLineLineWidth', 6);
set(groot,'defaultLineMarkerSize', 8);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
tlt1 = tiledlayout(1, 2);
title(tlt1, 'Efficiency comparison (excluding within-host tabulation time)', 'interpreter', 'latex', 'FontSize', 40);

% load file of interest: approximate time-driven algorithm:
nexttile 
for k = 1 : length(S_population)
    approx_filename = ['Approx_RunsDeltats_S', num2str(S_population(k)), 'I', ...
            num2str(I_population(k)), '_tendBH20_rep800_l1.5849e-10_k', num2str(k), '.mat'];
    load(approx_filename, 'runs_elapsed_time_approx', 'times_tabWH_approx');
    loglog(delta_ts, runs_elapsed_time_approx/rep, ':', 'Color', colors4{k}, 'DisplayName', ['S = ', num2str(S_population(k)), ' ppl.']);
    hold on
    ylim([10^(-5), 10^0])
    xlim([10^(-5), 10^2])
    hold on
end
% mark the critical resolution which guarantees accuracy:
for k = 1 : length(S_population)
    approx_filename = ['Approx_RunsDeltats_S', num2str(S_population(k)), 'I', ...
            num2str(I_population(k)), '_tendBH20_rep800_l1.5849e-10_k', num2str(k), '.mat'];
    load(approx_filename, 'runs_elapsed_time_approx', 'times_tabWH_approx');
    loglog(delta_ts(apprx_crit_stepsize(k)), runs_elapsed_time_approx(apprx_crit_stepsize(k))/800, '.', 'Color', colors4{k}, 'MarkerSize', 100, 'HandleVisibility', 'off'); hold on
    ylim([10^(-5), 10^0])
    xlim([10^(-5), 10^2])
end
legend('Location', 'northeast');
xlabel('$\Delta t$ (days)'); 
ylabel('Ave. runtime per repetition (seconds)'); 
title('Algorithm 1')

% load file of interest: exact algorithm:
nexttile 
for k = 1 : length(S_population)
    ExactMultiscale_filename = ['ExactMulti_RunsDeltats_S', num2str(S_population(k)), 'I', ...
            num2str(I_population(k)), '_tendBH20_rep800_l1.5849e-10_k', num2str(k), '.mat'];
    load(ExactMultiscale_filename, 'runs_elapsed_time_effiMulti', 'times_tabWH_exact');
    loglog(delta_ts, runs_elapsed_time_effiMulti/rep, 'Color', colors4{k}, 'DisplayName', ['S = ', num2str(S_population(k)), ' ppl.']);
    hold on
    ylim([10^(-5), 10^0])
    xlim([10^(-5), 10^2])
    hold on
end
% mark the critical resolution which guarantees accuracy:
for k = 1 : length(S_population)
    ExactMultiscale_filename = ['ExactMulti_RunsDeltats_S', num2str(S_population(k)), 'I', ...
            num2str(I_population(k)), '_tendBH20_rep800_l1.5849e-10_k', num2str(k), '.mat'];
    load(ExactMultiscale_filename, 'runs_elapsed_time_effiMulti', 'times_tabWH_exact');
    loglog(delta_ts(exact_crit_stepsize(k)), runs_elapsed_time_effiMulti(exact_crit_stepsize(k))/800, '*', 'Color', colors4{k}, 'MarkerSize', 40, 'HandleVisibility', 'off'); 
    hold on
    ylim([10^(-5), 10^0])
    xlim([10^(-5), 10^2])
end
legend('Location', 'northeast');
xlabel('$\delta$ (days)'); 
yticks([]);
title('Algorithm 2')




function [exact_crit_stepsize, exact_crit_err, apprx_crit_stepsize, apprx_crit_err] ...
    = find_crit_info(approx_error_mean, exact_error_mean, err_rel_threshold)
% This function finds the information at the critical Delta t value where two
% algorithms become inaccurate. 


experiment_num = 4; % simulations are performed over 4 different ICs

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

    apprx_crit_stepsize(i) = apprx_index;
    apprx_crit_err(i) = approx_error_mean(i, apprx_index);
    exact_crit_stepsize(i) = exact_index;
    exact_crit_err(i) = exact_error_mean(i, exact_index);
end
end

