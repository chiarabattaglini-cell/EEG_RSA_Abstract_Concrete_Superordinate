
function plot_cluster_markers(sig_idx, timecode)
% PLOT_CLUSTER_MARKERS  Draw vertical grey lines at the edges of each
%                       contiguous significant cluster.

    if isempty(sig_idx)
        return
    end

    gaps          = diff(sig_idx);
    cluster_start = [sig_idx(1);              sig_idx(find(gaps > 1) + 1)];
    cluster_end   = [sig_idx(find(gaps > 1)); sig_idx(end)];

    for c = 1:numel(cluster_start)
        xline(timecode(cluster_start(c)), '-', ...
              'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.7);
        xline(timecode(cluster_end(c)),   '-', ...
              'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.7);
    end
end