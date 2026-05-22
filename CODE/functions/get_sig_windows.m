
function str = get_sig_windows(sig_idx, timecode)
% GET_SIG_WINDOWS  Return a string listing significant time windows.

    if isempty(sig_idx)
        str = 'none';
        return
    end

    gaps          = diff(sig_idx);
    cluster_start = [sig_idx(1);              sig_idx(find(gaps > 1) + 1)];
    cluster_end   = [sig_idx(find(gaps > 1)); sig_idx(end)];

    parts = cell(numel(cluster_start), 1);
    for c = 1:numel(cluster_start)
        parts{c} = sprintf('[%.3f–%.3f]', ...
                           timecode(cluster_start(c)), timecode(cluster_end(c)));
    end
    str = strjoin(parts, '  ');
end