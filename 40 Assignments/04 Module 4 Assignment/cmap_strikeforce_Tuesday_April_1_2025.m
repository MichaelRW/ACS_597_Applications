

close all;  clear;  clc;


rng( 0 );
    cvars = randi( [10 450], [128, 1] );  % Column vector with random integers.
        % Minimum:  12
        % Maximum:  449


n = 512;
    cmap = jet(n);  % Colormap with 512 colors.


% cvars has a different range from cmap range (1:512), so you have to map it. (similar to imagesc)
cvars_map = ( cvars - min(cvars) ) / ( max(cvars) - min(cvars) ) * (n-1) + 1;



% Now find the color of i-th cvars
i = 2;             
    cvars( i )  % ans = 136
    cvars_map(i)  % mapped cvars 147
% ans = 147

cmap( cvars(i), : )   % corresponding color