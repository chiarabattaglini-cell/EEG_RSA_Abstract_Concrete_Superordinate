% =========================================================================
% ERP Grand Average and Time-Domain Plot
%
% Description:
%   Loads epoched EEG files (.set) for three semantic conditions
%   (abstract, concrete, superordinate), converts them from EEGLAB to
%   FieldTrip format, computes grand averages, and plots ERPs with
%   standard error shading for a cluster of electrodes of interest.
%
% Dependencies:
%   - EEGLAB (tested with eeglab2020_0)
%   - FieldTrip (tested with fieldtrip2020)
%   - boundedline (kakearney-boundedline-pkg)
%
% Input:
%   Epoched and baseline-corrected .set files organised in:
%   <MAINPATH>/TIME_DOMAIN_01/
%   Named as: *abstract.set, *concrete.set, *superordinate.set
%
% Output:
%   .tiff figure saved to PATHOUT
%
% Authors: [Chiara Battaglini Davide Bottari Martina Berto]
% Last modified: [22/06/26]
% =========================================================================

clear; close all

% -------------------------------------------------------------------------
% USER SETTINGS — edit these paths before running
% -------------------------------------------------------------------------

EEGLAB_PATH     = 'eeglab2020_0';
FIELDTRIP_PATH  = 'fieldtrip2020';
BOUNDEDLINE_PATH= 'kakearney-boundedline-pkg-8179f9a';

% Output subfolder name (created automatically if it does not exist)
OUTPUT_FOLDER   = 'Your folder name';

% Channel location file (relative to MAINPATH)
CHANLOC_FILE    = 'chanloc/GSN_HydroCel_65_noF.sfp';

% -------------------------------------------------------------------------
% ELECTRODE CLUSTER OF INTEREST
% -------------------------------------------------------------------------
% Fronto central electrodes
electrodes = [3, 4, 6, 7, 9, 12, 14, 15, 16, 20, 21, 41, 50, 51, 53, 54, 57, 60, 65];

% -------------------------------------------------------------------------
% PLOT COLOURS  [R G B], values in 0–1
% -------------------------------------------------------------------------
% Active selection: All three conditions
color = [0.97, 0.45, 0.45;   % abstract
         0.35, 0.38, 0.99;   % concrete
         0.32, 0.66, 0.18];  % superordinate

% Other colour sets (uncomment to use instead):
% % SUP vs CNC
% color = [0.35, 0.38, 0.99;   % concrete
%         0.32, 0.66, 0.18];  % superordinate
%
% % ABS vs CNC
% color = [0.97, 0.45, 0.45;   % abstract
%          0.35, 0.38, 0.99];  % concrete
%
% % ABS vs SUP
% color = [0.97, 0.45, 0.45;   % abstract
%          0.32, 0.66, 0.18];  % superordinate

% =========================================================================
% INITIALISE TOOLBOXES
% =========================================================================

addpath(EEGLAB_PATH);
eeglab;
close;                   % close GUI; paths are now set

addpath(FIELDTRIP_PATH);
ft_defaults;

addpath(BOUNDEDLINE_PATH);

% =========================================================================
% PATHS
% =========================================================================

MAINPATH = fileparts(pwd);
PATHIN   = fullfile(MAINPATH, 'TIME_DOMAIN_01', filesep);
PATHOUT  = fullfile(MAINPATH, OUTPUT_FOLDER, filesep);

if ~exist(PATHOUT, 'dir')
    mkdir(PATHOUT);
end

% =========================================================================
% LOAD DATA, CONVERT TO FIELDTRIP, COMPUTE GRAND AVERAGE
% =========================================================================

cd(PATHIN)

conditions = ["*abstract.set", "*concrete.set", "*superordinate.set"];

for c = 1:length(conditions)

    list    = dir(char(conditions(c)));          % find matching .set files
    list    = list(~strncmp({list.name},'._',2)); % exclude macOS hidden files
    subject = {list.name};

    for s = 1:length(subject)

        loadname = subject{s};
        EEG      = pop_loadset('filename', loadname, 'filepath', PATHIN);

        % Convert EEGLAB structure to FieldTrip
        data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

        % Compute trial average
        cfg            = [];
        cfg.channel    = 'all';
        cfg.keeptrials = 'no';
        timelock       = ft_timelockanalysis(cfg, data);

        % Attach electrode positions
        timelock.elec  = ft_read_sens(fullfile(MAINPATH, CHANLOC_FILE));

        array{s} = timelock;
        clear timelock data EEG

    end

    % Grand average across subjects
    cfg                  = [];
    cfg.keepindividual   = 'yes';
    grandavg             = ft_timelockgrandaverage(cfg, array{:});
    grandavg.elec        = ft_read_sens(fullfile(MAINPATH, CHANLOC_FILE));

    % Store individual data per condition
    switch c
        case 1
            individual_data_abs = grandavg.individual;
        case 2
            individual_data_con = grandavg.individual;
        case 3
            individual_data_sup = grandavg.individual;
            time                = grandavg.time;
    end

    clear grandavg array

end

% =========================================================================
% COMPUTE MEAN AND STANDARD ERROR OVER ELECTRODE CLUSTER
% =========================================================================

for i = 1:length(time)

    ABS = individual_data_abs(:, :, i);
    CON = individual_data_con(:, :, i);
    SUP = individual_data_sup(:, :, i);

    M_ABS(:, i) = mean(ABS(:, electrodes), 2);
    M_CON(:, i) = mean(CON(:, electrodes), 2);
    M_SUP(:, i) = mean(SUP(:, electrodes), 2);

    n = size(M_ABS, 1);
    SE_ABS(i) = std(M_ABS(:, i)) / sqrt(n);
    SE_CON(i) = std(M_CON(:, i)) / sqrt(n);
    SE_SUP(i) = std(M_SUP(:, i)) / sqrt(n);

end

% =========================================================================
% PLOT
% =========================================================================

figure;

% Active plot: All three conditions
l = boundedline(time, mean(M_ABS), SE_ABS, '-', ...
                time, mean(M_CON), SE_CON, '-', ...
                time, mean(M_SUP), SE_SUP, '-', ...
                'alpha', 'transparency', 0.15, 'cmap', color);


% Other plot combinations (uncomment to use instead):
% %  SUP vs CNC
% l = boundedline(time, mean(M_CON), SE_CON, '-', ...
%                 time, mean(M_SUP), SE_SUP, '-', ...
%                 'alpha', 'transparency', 0.15, 'cmap', color);
%
% % ABS vs CNC
% l = boundedline(time, mean(M_ABS), SE_ABS, '-', ...
%                 time, mean(M_CON), SE_CON, '-', ...
%                 'alpha', 'transparency', 0.15, 'cmap', color);
%
% % ABS vs SUP
% l = boundedline(time, mean(M_ABS), SE_ABS, '-', ...
%                 time, mean(M_SUP), SE_SUP, '-', ...
%                 'alpha', 'transparency', 0.15, 'cmap', color);

% Axes formatting
xlim([0 0.80]);
set(gca, 'FontSize', 16, 'Box', 'off', 'LineWidth', 1, 'YDir', 'reverse');
set(l, 'LineWidth', 2.5);

% Active legend: all three conditions
legend('abstract', 'concrete', 'superordinate', 'Location', 'southeast');  

% Other legends (uncomment to match chosen plot):
% legend('concrete', 'superordinate', 'Location', 'southeast');              %SUP vs CNC
% legend('abstract', 'concrete', 'Location', 'southeast');                   % ABS vs CNC
% legend('abstract', 'superordinate', 'Location', 'southeast');              % ABS vs SUP

% =========================================================================
% SAVE FIGURE
% =========================================================================

name = fullfile(PATHOUT, 'filename');
print(gcf, name, '-dtiff', '-r300');
