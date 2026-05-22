
% example
% 
% y1 and y2 are your data vector with avg of the group at each duration (5 columns)
% 
% e1 and e2 are your stand err vector with the group stnd err at each duration (5 columns)


color=[0 1 0.2
    1 0 0.2];
figure
[l,p] = boundedline(t, y1, e1, '-b',  t, y2, e2, '-r', 'alpha', 'transparency', 0.12, 'cmap', color);
%ylim([0 1]);
title('Right Auditory ROI');

