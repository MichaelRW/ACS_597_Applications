
function [ h_window ] = h_Hann_window( window_size )


% Modfied Hann Window

N = window_size;

h_Hann = hann( N, 'periodic' );

mean_square = sum( h_Hann.^2 ) / N;
average = sum( h_Hann ) / N;

h_Hann = h_Hann / sqrt( mean_square );

noise_bandwidth_factor = mean_square / average^2;


h_window.window = h_Hann;
h_window.enbw_factor = noise_bandwidth_factor;


