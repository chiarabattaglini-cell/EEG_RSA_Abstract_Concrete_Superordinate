% =========================================================================
% TIME-RESOLVED REPRESENTATIONAL SIMILARITY ANALYSIS (RSA)
% Simple and Partial Spearman Correlation with Permutation Testing
% =========================================================================
%
% Description:
%   This script performs time-resolved RSA by correlating a time-series
%   of EEG pairwise Euclidean distance matrices (RDMs) with model RDMs representing
%   different theoretical accounts of conceptual representation:
%     - Distributional semantics  (word2vec Euclidean pairwise distances)
%     - Sensorimotor features     (sensory-motor Euclidean pairwise distances)
%     - Abstraction               (abstraction ratings Euclidean pairwise distances)
%     - Abstractness              (abstractness ratings Euclidean pairwise distances)
%
%   Two analysis modes are available (set USE_PARTIAL_CORR below):
%     MODE 1 — Simple RSA:   Spearman correlation between EEG RDM and each
%                             model RDM independently at each timepoint.
%     MODE 2 — Partial RSA:  Spearman partial correlation controlling for
%                             all other model RDMs simultaneously.
%
% Statistical inference:
%   - Permutation test (n = 10,000): row/column permutation of model RDM
%   - Family-wise correction across timepoints and models: FDR-BH (q < 0.01)
%   - One-tailed p-value (positive correlation)
%
% Input files (all .mat, located in Data/RSA):
%   Dist_eu_EEG.mat          -> Deu_EEG   [timepoints x pairs]
%   Disteu_abstraction.mat   -> Abtion_pdist  [1 x pairs]
%   Abtness_pdist.mat        -> Abtness_pdist [1 x pairs]
%   word2vec_eu.mat          -> word2vec_pdist_eu [1 x pairs]
%   Dist_eu_sensory.mat      -> C          [1 x pairs]
%   timecode.mat             -> timecode   [1 x timepoints]
%
% Output:
%   RSA_results_[simple|partial].mat  saved in DATAPATH
%
% Dependencies:
%   fdr_bh.m  (Benjamini & Hochberg FDR correction)
%   Available at:
%   https://www.mathworks.com/matlabcentral/fileexchange/27418 and located
%   in CODE/functions
%
% Reference:
%   Kriegeskorte et al. (2008). Representational similarity analysis –
%   connecting the branches of systems neuroscience. Frontiers in Systems
%   Neuroscience, 2, 4, https://doi.org/10.3389/neuro.06.004.2008.
%
%   Benjamini & Hochberg (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
%   Journal of the Royal Statistical Society B, 57(1), 289-300, https://doi.org/10.1111/j.2517-6161.1995.tb02031.x.
%
% Authors:  [Chiara Battaglini Giacomo Handjaras]
% Date:    [22/05/26]
% Version: 1.0
% =========================================================================

%% 0. INITIALISATION
% -------------------------------------------------------------------------
clear; close all; clc


%% 1. USER PARAMETERS  ← edit only this section
% -------------------------------------------------------------------------

% Path to data files and output destination
DATAPATH = pwd;   % <-- change to your data folder

% Analysis mode:
%   false  →  Simple RSA    (independent Spearman correlations)
%   true   →  Partial RSA   (partial Spearman, all others as covariates)
USE_PARTIAL_CORR = false;

% Permutation test settings
N_PERM  = 10000;   % number of permutations
RNG_SEED = 42;     % random seed for reproducibility (set [] to disable)

% FDR correction settings
FDR_Q      = 0.01;         % FDR threshold
FDR_METHOD = 'pdep';       % 'pdep' (positive dependence) or 'dep' (general)

% Figure export
SAVE_FIGURE  = true;
FIG_FORMAT   = '-dpdf';    % '-dpng', '-depsc', '-dpdf'
FIG_RESOLUTION = '-r300';  % DPI (ignored for vector formats)


%% 2. LOAD DATA
% -------------------------------------------------------------------------
fprintf('[INFO] Loading data from: %s\n', DATAPATH);

tmp = load(fullfile(DATAPATH, 'Dist_eu_EEG.mat'));
Deu_EEG = tmp.Deu_EEG;
timepoints = size(Deu_EEG, 1);
n_pairs    = size(Deu_EEG, 2);

tmp = load(fullfile(DATAPATH, 'Dist_eu_abstraction.mat'));
abstraction = tmp.Abtion_pdist(:)';         % enforce row vector

tmp = load(fullfile(DATAPATH, 'Dist_eu_abstractness.mat'));
abstractness = tmp.Abtness_pdist(:)';

tmp = load(fullfile(DATAPATH, ['Dist_eu_word2vec' ...
    '.mat']));
word2vec = tmp.word2vec_pdist_eu(:)';

tmp = load(fullfile(DATAPATH, 'Dist_eu_sensory.mat'));
sensory_motor = tmp.C(:)';

tmp = load(fullfile(DATAPATH, 'timecode.mat'));
timecode = tmp.timecode(:)';



%% 3. FIXED PERMUTATION MATRIX
% -------------------------------------------------------------------------
% Pre-generate all permutations once so results are exactly reproducible
% when the same seed is used across runs.

if ~isempty(RNG_SEED)
    rng(RNG_SEED, 'twister');
    fprintf('[INFO] Random seed set to %d\n', RNG_SEED);
end

fprintf('[INFO] Generating permutation matrix (%d x %d)...\n', n_pairs, N_PERM);
perm_matrix = zeros(n_pairs, N_PERM, 'uint32');
for p = 1:N_PERM
    perm_matrix(:, p) = randperm(n_pairs);
end


%% 4. MODEL DEFINITIONS
% -------------------------------------------------------------------------
% Collect models in a struct array for compact, loop-friendly processing.
% To add/remove models, edit only this section.

models(1).name   = 'word2vec';
models(1).label  = 'Distributional';
models(1).rdm    = word2vec;
models(1).color  = [0.1797  0.5430  0.3398];

models(2).name   = 'sensory_motor';
models(2).label  = 'Sensorimotor';
models(2).rdm    = sensory_motor;
models(2).color  = [0.8516  0.6445  0.1250];

models(3).name   = 'abstraction';
models(3).label  = 'Abstraction';
models(3).rdm    = abstraction;
models(3).color  = [0.9961  0.4102  0.7031];

models(4).name   = 'abstractness';
models(4).label  = 'Abstractness';
models(4).rdm    = abstractness;
models(4).color  = [0.6953  0.1328  0.1328];

n_models = numel(models);


%% 5. TIME-RESOLVED CORRELATION WITH PERMUTATION NULL DISTRIBUTION
% -------------------------------------------------------------------------
% This computation takes some time
fprintf('[INFO] Computing correlations (%d timepoints x %d permutations)...\n', ...
        timepoints, N_PERM);

% Pre-allocate
for m = 1:n_models
    models(m).rho      = nan(timepoints, 1);
    models(m).rho_perm = nan(timepoints, N_PERM);
end

t_start = tic;

for t = 1:timepoints
t % this gives info on the timepoint that is being processed
    eeg_rdm = Deu_EEG(t, :)';    % [pairs x 1]

    % --- Observed correlations ---
    for m = 1:n_models
        models(m).rho(t) = compute_corr(eeg_rdm, models(m).rdm', ...
                                         models, m, USE_PARTIAL_CORR);
    end

    % --- Permuted null distributions ---
    for p = 1:N_PERM
        idx = perm_matrix(:, p);
        for m = 1:n_models
            models(m).rho_perm(t, p) = compute_corr( ...
                eeg_rdm, models(m).rdm(idx)', models, m, USE_PARTIAL_CORR, idx);
        end
    end

    % Progress report every 50 timepoints
    if mod(t, 50) == 0
        elapsed = toc(t_start);
        eta     = elapsed / t * (timepoints - t);
        fprintf('  t = %4d / %d   (elapsed: %.0fs, ETA: %.0fs)\n', ...
                t, timepoints, elapsed, eta);
    end
end

fprintf('[INFO] Correlation step completed in %.1f s\n', toc(t_start));


%% 6. PERMUTATION P-VALUES (one-tailed, positive)
% -------------------------------------------------------------------------
fprintf('[INFO] Computing permutation p-values...\n');

for m = 1:n_models
    models(m).pval = nan(timepoints, 1);
    for t = 1:timepoints
        null_dist = [models(m).rho_perm(t, :), models(m).rho(t)];
        null_sort = sort(null_dist, 'descend');
        rank_pos  = find(null_sort == models(m).rho(t), 1, 'last');
        models(m).pval(t) = rank_pos / (N_PERM + 1);
    end
end


%% 7. FDR CORRECTION (across all timepoints and models jointly)
% -------------------------------------------------------------------------
fprintf('[INFO] Applying FDR correction (q = %.2f, method = %s)...\n', ...
        FDR_Q, FDR_METHOD);

all_pvals = cell2mat(arrayfun(@(m) m.pval, models, 'UniformOutput', false));
% all_pvals: [timepoints x n_models]

[~, crit_p, adj_p_all] = fdr_bh(all_pvals(:), FDR_Q, FDR_METHOD, 'yes');

fprintf('  Critical p-value after FDR: %.6f\n', crit_p);

% Distribute corrected p-values back to each model
adj_p_mat = reshape(adj_p_all, timepoints, n_models);
for m = 1:n_models
    models(m).pval_adj = adj_p_mat(:, m);
    % Significant timepoints (using the common critical threshold)
    models(m).sig_idx  = find(models(m).pval <= crit_p);
    % Corrected rho (NaN at non-significant timepoints)
    models(m).rho_sig  = models(m).rho;
    models(m).rho_sig(models(m).pval > crit_p) = NaN;
end


%% 8. SAVE RESULTS
% -------------------------------------------------------------------------
if USE_PARTIAL_CORR
    out_suffix = 'partial';
else
    out_suffix = 'simple';
end

outfile = fullfile(DATAPATH, sprintf('RSA_results_%s.mat', out_suffix));
save(outfile, 'models', 'timecode', 'crit_p', 'adj_p_all', ...
     'N_PERM', 'RNG_SEED', 'FDR_Q', 'USE_PARTIAL_CORR', '-v7.3');
fprintf('[INFO] Results saved to: %s\n', outfile);


%% 9. PRINT SUMMARY TABLE
% -------------------------------------------------------------------------
fprintf('\n=== RSA RESULTS SUMMARY (%s) ===\n', upper(out_suffix));
fprintf('%-20s  %-12s  %-12s  %-20s\n', ...
        'Model', 'Peak rho', 'Peak time(s)', 'Significant windows');
fprintf('%s\n', repmat('-', 1, 70));

for m = 1:n_models
    [pk, pk_idx] = max(models(m).rho);
    sig_wins = get_sig_windows(models(m).sig_idx, timecode);
    fprintf('%-20s  %+.4f       %.3f s        %s\n', ...
            models(m).label, pk, timecode(pk_idx), sig_wins);
end
fprintf('\n');


%% 10. FIGURE
% -------------------------------------------------------------------------
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 12]);
hold on;

% --- Faded full time-course ---
for m = 1:n_models
    plot(timecode, models(m).rho, ...
         'Color', [models(m).color, 0.25], 'LineWidth', 2);
end

% --- Opaque significant portions ---
h_leg = gobjects(n_models, 1);
for m = 1:n_models
    h_leg(m) = plot(timecode, models(m).rho_sig, ...
                    'Color', models(m).color, 'LineWidth', 3);
end

% --- Vertical markers for significant clusters of Abstraction ---
% (shown for abstraction only; adjust model index as needed)
m_abs = find(strcmp({models.name}, 'abstraction'));
plot_cluster_markers(models(m_abs).sig_idx, timecode);

% --- Zero line ---
yline(0, '--k', 'LineWidth', 0.8, 'Alpha', 0.4);

% --- Axes formatting ---
ax = gca;
set(ax, 'Box', 'off', 'FontName', 'Palatino', 'FontSize', 12, ...
        'TickLabelInterpreter', 'none', ...
        'YTick', [-0.1, 0, 0.1, 0.2], ...
        'YTickLabel', {'-0.1', '0', '0.1', '0.2'});
set(gcf, 'Color', 'w');

xlabel('Time (s)',     'FontSize', 14, 'FontName', 'Palatino');
if USE_PARTIAL_CORR
    ylabel('Partial correlation (\rho)', 'FontSize', 14, 'FontName', 'Palatino');
    title_str = 'Time-resolved partial RSA';
else
    ylabel('Spearman correlation (\rho)', 'FontSize', 14, 'FontName', 'Palatino');
    title_str = 'Time-resolved RSA';
end
title(title_str, 'FontSize', 16, 'FontName', 'Palatino');

legend(h_leg, {models.label}, ...
       'FontSize', 14, 'FontName', 'Palatino', ...
       'Location', 'northeastoutside');
legend('boxoff');

pbaspect([1.5, 1, 1]);
hold off;

% --- Export figure ---
if SAVE_FIGURE
    fig_name = fullfile(DATAPATH, sprintf('RSA_figure_%s', out_suffix));
    print(fig, fig_name, FIG_FORMAT, FIG_RESOLUTION);
    fprintf('[INFO] Figure saved to: %s\n', fig_name);
end





