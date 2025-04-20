
function [ spectrogram ] = spectrogram_September_26_2023( signal, frame_length, frame_hop, sample_rate, APPLY_HANN_WINDOW )


frames = frame_set( signal, frame_length, frame_hop, sample_rate );

N = size( frames.frame_set, 2 );
    number_of_frames = frames.count;

nyquist_frequency = sample_rate / 2;

spectral_resolution_Hz = sample_rate / N;
    frequencies = 0:spectral_resolution_Hz:( sample_rate / 2 );
        frequencies = [ -fliplr( frequencies( 2:1:end - 1 ) )  frequencies ];

sample_period = 1 / sample_rate;


if ( APPLY_HANN_WINDOW == 0 )
    window = ones( 1, N );  % "Box" or "Rectangular" window.
else
    h_window = h_Hann_window( N );  % A periodic window.
        window = h_window.window';  % Change to a row vector.
        ebnw_factor = h_window.enbw_factor;    
end


% Pre-allocated variable for individiual frame linear spectrum calculations.
linear_spectrum.frames = nan( N, number_of_frames );

% Compute Gxx per frame.
for frame_index = 1:1:number_of_frames
    linear_spectrum.frames( :, frame_index ) = ( fft( frames.frame_set( frame_index, : ) .* window ) * sample_period ).';
end


spectrogram.frequencies = frequencies';


% Compute the double-sided and single-sided power spectral densities.
Sxx = nan( numel( frequencies ), number_of_frames );
Gxx = nan( N / 2 + 1, number_of_frames );

net_time = numel( signal ) * sample_period;
    time_indices = 0:( N / 2 * sample_period):net_time;


T = N * sample_period;  % Record duration.

for frame_index = 1:1:number_of_frames
    Sxx( :, frame_index ) = ( ( 1 / T ) .* conj( linear_spectrum.frames( :, frame_index ) ) .* linear_spectrum.frames( :, frame_index ) );
        Gxx( 2:1:( end - 1 ), frame_index ) = Sxx( 2:1:( N / 2 ), frame_index ) .* 2;
            Gxx( 1, frame_index ) = Sxx( 1, frame_index );
            Gxx( end, frame_index ) = Sxx( N / 2 + 1, frame_index );
end


% Set values in results structure.
spectrogram.sample_rate = sample_rate;
spectrogram.sample_period = sample_period;
spectrogram.nyquist_frequency = nyquist_frequency;
spectrogram.spectral_resolution = spectral_resolution_Hz;
spectrogram.Sxx.spectrum = real( circshift( Sxx, N / 2 - 1 ) );
spectrogram.Sxx.frequencies = frequencies;
spectrogram.Gxx.spectrum = real( Gxx );
spectrogram.Gxx.frequencies = frequencies( N/2:1:end );
spectrogram.time_indices = time_indices;



%% Reference(s)


