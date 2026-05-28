% =========================================================================
% CLUSTER-BASED PERMUTATION TEST FOR ERP DATA (FieldTrip)
% =========================================================================
%
% Description:
%   This script performs a cluster-based non-parametric permutation test
%   (Monte Carlo approximation) on ERP data across three
%   conditions: Abstract, Concrete, and Superordinate words.
%   A dependent-samples multivariate F-statistic is used to assess
%   differences in scalp topography over time.
%
% Statistical approach:
%   - Test:            Dependent-samples multivariate F (one-way repeated
%                      measures ANOVA across 3 conditions)
%   - Correction:      Cluster-based permutation (maxsum statistic)
%   - Permutations:    10000
%   - Alpha (cluster): 0.05 (one-tailed)
%   - Min. neighbours: 2 channels
%
% Input data:
%   - Grand-average ERP structures WITH individual subject data
%     (keepindividual = 'yes' in ft_timelockgrandaverage)
%   - Grand-average ERP structures WITHOUT individual subject data
%     (used for visualisation only)
%
% Output:
%   - stat: FieldTrip statistics structure saved as .mat file
%
% Requirements:
%   - FieldTrip toolbox (tested with fieldtrip-20200130)
%   - Electrode layout file: GSN_HydroCel_65_noF.sfp (Available in
%   Data/EEG/chanloc)
%
% Reference:
%   Maris & Oostenveld (2007). Nonparametric statistical testing of
%   EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190, doi: 10.1016/j.jneumeth.2007.03.024.
%
% Authors:  [Chiara Battaglini Davide Bottari Martina Berto]
% Date:    [22/05/26]
% Version: 1.0
% =========================================================================

%% 0. INITIALISATION
% -------------------------------------------------------------------------
clear; close all; clc

% Add FieldTrip to path and initialise defaults
% addpath('/path/to/fieldtrip')   % <-- set your FieldTrip path here
ft_defaults


%% 1. PATHS
% -------------------------------------------------------------------------
% All paths are defined relative to MAINPATH so that the script is
% portable across machines. Adjust MAINPATH if needed.

MAINPATH = fileparts(pwd);

% Input: grand-average structures WITH individual subject data
PATHIN_INDIV  = fullfile(MAINPATH, '03_TIME_DOMAIN_GA_keepSUB', filesep);

% Input: grand-average structures WITHOUT individual subject data
PATHIN_AVG    = fullfile(MAINPATH, '03_TIME_DOMAIN_GA',         filesep);

% Output: cluster permutation results
PATHOUT       = fullfile(MAINPATH, 'your output folder', filesep);

% Channel layout and electrode positions
CHANPATH      = fullfile(MAINPATH, 'chanloc', 'GSN_HydroCel_65_noF.sfp');

% Create output directory if absent
if ~exist(PATHOUT, 'dir')
    mkdir(PATHOUT);
    fprintf('[INFO] Created output directory: %s\n', PATHOUT);
end


%% 2. PARAMETERS
% -------------------------------------------------------------------------
% Conditions to compare (must match filenames)
conditions = {'abstract', 'concrete', 'superordinate'};
N_COND     = numel(conditions);

% Channels of interest (pre-selected a priori based on prior literature;
% see manuscript Methods for rationale)
CHANNELS_OI = {'E3',  'E4',  'E6',  'E7',  'E9',  'E12', ...
                'E14', 'E15', 'E16', 'E20', 'E21', 'E41', ...
                'E50', 'E51', 'E53', 'E54', 'E57', 'E60', 'E65'};

% Statistical parameters
STAT_PARAMS.method           = 'montecarlo';
STAT_PARAMS.statistic        = 'ft_statfun_depsamplesFmultivariate';
STAT_PARAMS.correctm         = 'cluster';
STAT_PARAMS.clusteralpha     = 0.05;
STAT_PARAMS.clusterstatistic = 'maxsum';
STAT_PARAMS.minnbchan        = 2;
STAT_PARAMS.tail             = 1;    % one-tailed (F-statistic is always >= 0)
STAT_PARAMS.clustertail      = 1;
STAT_PARAMS.alpha            = 0.05;
STAT_PARAMS.numrandomization = 1000;
STAT_PARAMS.latency          = [];   % full epoch; set to [tmin tmax] to restrict

% Neighbourhood distance (in units of the layout coordinates)
NEIGH_DIST = 0.11;


%% 3. LOAD GRAND-AVERAGE DATA
% -------------------------------------------------------------------------
fprintf('\n[INFO] Loading grand-average ERP data...\n');

grandavg_indiv = cell(1, N_COND);   % with individual trials (for stats)
grandavg_avg   = cell(1, N_COND);   % average only         (for plots)

for c = 1:N_COND
    cond = conditions{c};

    % -- With individual data --
    fname = fullfile(PATHIN_INDIV, sprintf('Gavg_ERP_%s.mat', cond));
    tmp   = load(fname, 'grandavg');
    grandavg_indiv{c} = tmp.grandavg;
    fprintf('  Loaded (indiv): %s  [%d subjects]\n', fname, ...
            size(grandavg_indiv{c}.individual, 1));

    % -- Average only --
    fname = fullfile(PATHIN_AVG, sprintf('Gavg_ERP_%s.mat', cond));
    tmp   = load(fname, 'grandavg');
    grandavg_avg{c} = tmp.grandavg;
    fprintf('  Loaded (avg):   %s\n', fname);
end

% Convenience aliases
grandavg_ABS  = grandavg_indiv{1};
grandavg_CNC  = grandavg_indiv{2};
grandavg_SUP  = grandavg_indiv{3};

n_subjects = size(grandavg_ABS.individual, 1);
fprintf('\n[INFO] Number of subjects: %d\n', n_subjects);


%% 4. ELECTRODE LAYOUT AND NEIGHBOURHOOD STRUCTURE
% -------------------------------------------------------------------------
fprintf('\n[INFO] Preparing electrode layout and neighbourhood structure...\n');

cfg_lay        = [];
cfg_lay.layout = CHANPATH;
lay            = ft_prepare_layout(cfg_lay);

% Visual check of layout (comment out for batch processing)
cfg_plot        = [];
cfg_plot.layout = lay;
ft_layoutplot(cfg_plot);

% Define neighbourhoods based on spatial distance
cfg_neigh              = [];
cfg_neigh.method       = 'distance';
cfg_neigh.neighbourdist= NEIGH_DIST;
cfg_neigh.elec         = ft_read_sens(CHANPATH);
cfg_neigh.feedback     = 'yes';   % set to 'no' to suppress figures
neighbours             = ft_prepare_neighbours(cfg_neigh, grandavg_ABS);


%% 5. DESIGN MATRIX
% -------------------------------------------------------------------------
% Rows:
%   Row 1 (uvar): subject index   [1 ... n_subjects repeated N_COND times]
%   Row 2 (ivar): condition index [1 ... N_COND, each block of n_subjects]
%
% Layout (N_COND = 3, n_subjects = S):
%   [ 1  2  ...  S   1  2  ...  S   1  2  ...  S  ]
%   [ 1  1  ...  1   2  2  ...  2   3  3  ...  3  ]

design = [repmat(1:n_subjects, [1, N_COND]); ...
          repelem(1:N_COND, n_subjects)];

fprintf('\n[INFO] Design matrix: %d conditions x %d subjects (%d observations total)\n', ...
        N_COND, n_subjects, size(design, 2));


%% 6. CLUSTER-BASED PERMUTATION TEST
% -------------------------------------------------------------------------
fprintf('\n[INFO] Running cluster-based permutation test...\n');
fprintf('       Statistic : %s\n',   STAT_PARAMS.statistic);
fprintf('       Correction : %s\n',  STAT_PARAMS.correctm);
fprintf('       N perm     : %d\n',  STAT_PARAMS.numrandomization);
fprintf('       Alpha      : %.3f\n',STAT_PARAMS.alpha);
fprintf('       Tail       : %d\n',  STAT_PARAMS.tail);

cfg_stat                  = STAT_PARAMS;    % copy all statistical params
cfg_stat.neighbours       = neighbours;
cfg_stat.channel          = CHANNELS_OI;
cfg_stat.design           = design;
cfg_stat.uvar             = 1;              % row of subject indices
cfg_stat.ivar             = 2;              % row of condition indices

stat = ft_timelockstatistics(cfg_stat, grandavg_ABS, grandavg_SUP, grandavg_CNC);


%% 7. SAVE RESULTS
% -------------------------------------------------------------------------
outname = fullfile(PATHOUT, 'ABS_SUP_CNC__cluster_permutation.mat');
save(outname, 'stat', 'cfg_stat', 'design', 'conditions', 'CHANNELS_OI', '-v7.3');
fprintf('\n[INFO] Results saved to: %s\n', outname);
