%% =========================================================================
%  CLUSTER-BASED PERMUTATION TEST & TOPOGRAPHIC VISUALIZATION
%  -------------------------------------------------------------------------
%  This script:
%    1. Loads grand-average ERP data for three semantic conditions
%       (abstract, concrete, superordinate)
%    2. Builds a channel neighbourhood structure
%    3. Runs a non-parametric cluster-based permutation test (Monte Carlo)
%       between two user-specified conditions over a chosen time window
%    4. Plots grand-average difference topographies with significant
%       cluster channels highlighted
%
%  Three latency windows were tested in the reference paper:
%    (1) Full epoch   : 0 – 1.058 s
%    (2) Early window : 0.250 – 0.500 s
%    (3) Late window  : 0.650 – 0.750 s
%
%  Change LATENCY_PRESET (USER PARAMETERS section) to switch windows
%  without modifying any other part of the script.
%
%  Dependencies : FieldTrip toolbox (ft_defaults must be on MATLAB path)
%  Reference    : Maris & Oostenveld (2007). Nonparametric statistical
%                 testing of EEG- and MEG-data. Journal of Neuroscience
%                 Methods, 164(1), 177–190.
%  Authors      : [Chiara Battaglini, Davide Bottari, Martina Berto]
%  Last updated : [22/05/26]
% =========================================================================

clear; close all; clc;

%% =========================================================================
%  USER PARAMETERS  –  edit this section only
% =========================================================================

% --- FieldTrip path -------------------------------------------------------
FIELDTRIP_PATH = '';    % leave empty if already on MATLAB path
                        % e.g. '/Users/username/fieldtrip2020'

% --- Paths ----------------------------------------------------------------
OUTPUT_FOLDER_STAT = 'stat_output';    % subfolder for permutation results
OUTPUT_FOLDER_PLOT = 'stat_plots';     % subfolder for figure output
OUTPUT_FILENAME    = 'stat_result.mat';% filename for permutation results

% --- Channel location file (relative to MAINPATH) ------------------------
CHANLOC_FILE = 'chanloc/GSN_HydroCel_65_noF.sfp';
LAYOUT_FILE  = 'chanloc/GSN-HydroCel-65_lay.mat';

% --- Contrast to test -----------------------------------------------------
% Choose which two conditions to compare.
% Options: 'ABS_vs_CNC' | 'ABS_vs_SUP' | 'CNC_vs_SUP'
CONTRAST = 'ABS_vs_CNC';

% --- Latency window (seconds) --------------------------------------------
% Choose ONE of the three windows used in the paper, or define a custom one.
%
%   'full'   : 0 – 1.058 s  (entire post-stimulus epoch)
%   'early'  : 0.250 – 0.500 s
%   'late'   : 0.650 – 0.750 s
%   'custom' : set CUSTOM_WIN below
%
LATENCY_PRESET = 'late';        % 'full' | 'early' | 'late' | 'custom'
CUSTOM_WIN     = [0.65 0.75];   % used only when LATENCY_PRESET = 'custom'

% --- Electrode cluster of interest ----------------------------------------
% Fronto-central electrodes (GSN-HydroCel-65 montage)
CHANNELS = {'E3','E4','E6','E7','E9','E12','E14','E15','E16', ...
            'E20','E21','E41','E50','E51','E53','E54','E57','E60','E65'};

% --- Statistical parameters -----------------------------------------------
NEIGHBOUR_DIST = 0.11;   % neighbourhood distance threshold
CLUSTER_ALPHA  = 0.05;   % cluster-forming alpha (sample level)
TEST_ALPHA     = 0.025;  % permutation test alpha (one-tailed)
TAIL           = 0;      % 0 = two-sided, 1 = right, -1 = left
N_PERMUTATIONS = 1000;   % number of Monte Carlo draws
MIN_NB_CHAN    = 2;      % minimum neighbours to form a cluster

% --- Topography plot parameters -------------------------------------------
ZLIM           = [-1 0.5];  % colour axis (µV); use 'maxmin' for auto
TIMESTEP       = 0.01;      % bin width for topography plots (s)
SAMPLING_RATE  = 250;       % recording sampling rate (Hz)
SUBPLOT_ROWS   = 5;
SUBPLOT_COLS   = 5;

% =========================================================================
%  END OF USER PARAMETERS
% =========================================================================


%% -------------------------------------------------------------------------
%  INITIALISE FIELDTRIP
% -------------------------------------------------------------------------

if ~isempty(FIELDTRIP_PATH)
    addpath(FIELDTRIP_PATH);
end
ft_defaults


%% -------------------------------------------------------------------------
%  PATHS
% -------------------------------------------------------------------------

MAINPATH = fileparts(pwd);
PATHIN   = fullfile(MAINPATH, '03_TIME_DOMAIN_GA_keepSUB', filesep); % with individual data
PATHIN1  = fullfile(MAINPATH, '03_TIME_DOMAIN_GA',         filesep); % condition averages only
PATHSTAT = fullfile(MAINPATH, OUTPUT_FOLDER_STAT, filesep);
PATHOUT  = fullfile(MAINPATH, OUTPUT_FOLDER_PLOT, filesep);

if ~exist(PATHSTAT, 'dir'); mkdir(PATHSTAT); end
if ~exist(PATHOUT,  'dir'); mkdir(PATHOUT);  end


%% -------------------------------------------------------------------------
%  RESOLVE LATENCY WINDOW
% -------------------------------------------------------------------------

switch lower(LATENCY_PRESET)
    case 'full'
        t_start = 0;
        t_end   = 1.058;
    case 'early'
        t_start = 0.250;
        t_end   = 0.500;
    case 'late'
        t_start = 0.650;
        t_end   = 0.750;
    case 'custom'
        t_start = CUSTOM_WIN(1);
        t_end   = CUSTOM_WIN(2);
    otherwise
        error('Unknown LATENCY_PRESET. Use ''full'', ''early'', ''late'', or ''custom''.');
end

LATENCY = [t_start t_end];
j = t_start : TIMESTEP : t_end;               % time bin edges (s)
fprintf('Latency window: %.3f – %.3f s\n', t_start, t_end);


%% =========================================================================
%  SECTION 1 – LOAD GRAND AVERAGES
% =========================================================================

conditions = {'abstract', 'concrete', 'superordinate'};

% Pre-allocate as struct arrays (one entry per condition)
grandavg_cond    = struct([]);   % includes individual subject data
grandavg_CONDavg = struct([]);   % condition-level averages only

for c = 1 : length(conditions)

    % With individual data (needed for the permutation test)
    fname = fullfile(PATHIN, ['Gavg_ERP_' conditions{c} '.mat']);
    tmp   = load(fname, 'grandavg');
    grandavg_cond(c) = tmp.grandavg;

    % Without individual data (needed for topography plots)
    fname = fullfile(PATHIN1, ['Gavg_ERP_' conditions{c} '.mat']);
    tmp   = load(fname, 'grandavg');
    grandavg_CONDavg(c) = tmp.grandavg;

end

% Named variables for readability
grandavg_ABS      = grandavg_cond(1);
grandavg_CNC      = grandavg_cond(2);
grandavg_SUP      = grandavg_cond(3);

grandavg_COND_ABS = grandavg_CONDavg(1);
grandavg_COND_CNC = grandavg_CONDavg(2);
grandavg_COND_SUP = grandavg_CONDavg(3);

% Number of subjects (from individual data array)
n_subj = size(grandavg_ABS.individual, 1);
fprintf('Number of subjects: %d\n', n_subj);


%% =========================================================================
%  SECTION 2 – CHANNEL LAYOUT AND NEIGHBOURHOOD STRUCTURE
% =========================================================================

% Prepare layout for topography plots
load(fullfile(MAINPATH, LAYOUT_FILE), 'lay');
cfg_lay        = [];
cfg_lay.layout = lay;
lay            = ft_prepare_layout(cfg_lay);

% Compute channel neighbourhood structure
cfg_nb               = [];
cfg_nb.method        = 'distance';
cfg_nb.neighbourdist = NEIGHBOUR_DIST;
cfg_nb.elec          = ft_read_sens(fullfile(MAINPATH, CHANLOC_FILE));
cfg_nb.feedback      = 'yes';
neighbours           = ft_prepare_neighbours(cfg_nb, grandavg_ABS);


%% =========================================================================
%  SECTION 3 – CLUSTER-BASED PERMUTATION TEST
% =========================================================================

% Select conditions for the chosen contrast
switch upper(CONTRAST)
    case 'ABS_VS_CNC'
        condA = grandavg_ABS;  condB = grandavg_CNC;
        contrast_label = 'ABSvsCNC';
    case 'ABS_VS_SUP'
        condA = grandavg_ABS;  condB = grandavg_SUP;
        contrast_label = 'ABSvsSUP';
    case 'CNC_VS_SUP'
        condA = grandavg_CNC;  condB = grandavg_SUP;
        contrast_label = 'CNCSvsSUP';
    otherwise
        error('Unknown CONTRAST. Use ''ABS_vs_CNC'', ''ABS_vs_SUP'', or ''CNC_vs_SUP''.');
end

fprintf('Running permutation test: %s\n', contrast_label);

% Within-subject design matrix:
%   Row 1 – subject indices (repeated across the two conditions)
%   Row 2 – condition labels (1 = condition A, 2 = condition B)
design = [repmat(1:n_subj, [1, 2]); ...
          ones(1, n_subj),           2*ones(1, n_subj)];

cfg_stat                  = [];
cfg_stat.neighbours       = neighbours;
cfg_stat.channel          = CHANNELS;
cfg_stat.latency          = LATENCY;
cfg_stat.method           = 'montecarlo';
cfg_stat.statistic        = 'depsamplesT';   % dependent-samples t-statistic
cfg_stat.correctm         = 'cluster';
cfg_stat.clusteralpha     = CLUSTER_ALPHA;
cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.minnbchan        = MIN_NB_CHAN;
cfg_stat.tail             = TAIL;
cfg_stat.clustertail      = TAIL;
cfg_stat.alpha            = TEST_ALPHA;
cfg_stat.numrandomization = N_PERMUTATIONS;
cfg_stat.design           = design;
cfg_stat.uvar             = 1;   % unit variable  (subjects)
cfg_stat.ivar             = 2;   % independent variable (conditions)

stat = ft_timelockstatistics(cfg_stat, condA, condB);

% Save permutation test results
stat_out = fullfile(PATHSTAT, OUTPUT_FILENAME);
save(stat_out, 'stat');
fprintf('Permutation results saved: %s\n', stat_out);


%% =========================================================================
%  SECTION 4 – TOPOGRAPHIC VISUALIZATION
% =========================================================================

% Sample indices of time-bin edges
sample_count = length(stat.time);
m = 1 : TIMESTEP * SAMPLING_RATE : sample_count;

% Select grand-average structures matching the chosen contrast
switch upper(CONTRAST)
    case 'ABS_VS_CNC'
        GA_A = grandavg_COND_ABS;  GA_B = grandavg_COND_CNC;
    case 'ABS_VS_SUP'
        GA_A = grandavg_COND_ABS;  GA_B = grandavg_COND_SUP;
    case 'CNC_VS_SUP'
        GA_A = grandavg_COND_CNC;  GA_B = grandavg_COND_SUP;
end

% Grand-average difference map (A minus B)
cfg_math           = [];
cfg_math.operation = 'subtract';
cfg_math.parameter = 'avg';
GA_diff = ft_math(cfg_math, GA_A, GA_B);

% Align channel order between difference map and statistical output
[i1, i2] = match_str(GA_diff.label, stat.label);


% -------------------------------------------------------------------------
%  FIGURE 1 – POSITIVE CLUSTERS  (condition A > condition B)
% -------------------------------------------------------------------------

pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust         = find(pos_cluster_pvals < TEST_ALPHA);
pos_mask          = ismember(stat.posclusterslabelmat, pos_clust);

fig1 = figure('Position', [71, 110, 1047, 687]);

for k = 1 : length(j) - 1
    subplot(SUBPLOT_ROWS, SUBPLOT_COLS, k);

    % Channels significant throughout the entire bin
    pos_int           = zeros(numel(GA_diff.label), 1);
    pos_int(i1)       = all(pos_mask(i2, m(k) : m(k+1)), 2);

    cfg_plot                  = [];
    cfg_plot.xlim             = [j(k)  j(k+1)];
    cfg_plot.zlim             = ZLIM;
    cfg_plot.showlabels       = 'yes';
    cfg_plot.marker           = 'off';
    cfg_plot.colorbar         = 'yes';
    cfg_plot.highlight        = 'on';
    cfg_plot.highlightchannel = find(pos_int);
    cfg_plot.highlightsize    = 8;
    cfg_plot.comment          = 'xlim';
    cfg_plot.commentpos       = 'middlebottom';
    cfg_plot.layout           = lay;

    ft_topoplotER(cfg_plot, GA_diff);
end

out_pos = fullfile(PATHOUT, ...
    sprintf('%s_%s_%dms-%dms_POS', ...
            contrast_label, LATENCY_PRESET, ...
            round(t_start*1000), round(t_end*1000)));
print(fig1, out_pos, '-dtiff', '-r300');
fprintf('Saved: %s.tif\n', out_pos);

