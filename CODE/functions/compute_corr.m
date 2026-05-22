
function rho = compute_corr(eeg_col, model_col, models, focal_idx, use_partial, perm_idx)
% COMPUTE_CORR  Spearman (partial) correlation between EEG RDM and one model.
%
%   eeg_col    : [pairs x 1] EEG distances at one timepoint
%   model_col  : [pairs x 1] focal model distances (possibly permuted)
%   models     : struct array of all models
%   focal_idx  : index of the focal model
%   use_partial: logical — partial correlation if true
%   perm_idx   : (optional) permutation index applied to focal model only
%
% In partial mode, covariate models keep their ORIGINAL order (unpermuted),
% because only the focal model is under the null hypothesis.

    if nargin < 6
        perm_idx = [];
    end

    if ~use_partial
        rho = corr(eeg_col, model_col, 'type', 'spearman');
    else
        % Build covariate matrix from all models except the focal one
        n_models = numel(models);
        cov_idx  = setdiff(1:n_models, focal_idx);
        covariates = cell2mat(arrayfun(@(i) models(i).rdm(:), cov_idx, ...
                                       'UniformOutput', false));
        rho = partialcorr(eeg_col, model_col, covariates, 'type', 'spearman');
    end
end
