
function [ frames ] = frame_set( signal, frame_length, frame_hop, sample_rate )


%% Synopsis

% This function generates a set of frames for a time series.

% Edge Cases:
%
%   Checked - 1.) NF = N -> Size of the frame is equal to the length of the signal.
%
%   Checked - 2.) NF = NH = 1 -> Frame size of 1.  There should be N frames.
%
%   Checked - 3.) NF = 1 and NH = 2 -> Taking every other sample, decimating by 2.  N/2 frames if N
%   is even.  (N - 1)/2 if N is odd.



%% Processing

N = numel( signal );

nyquist_frequency = sample_rate / 2;

sample_period = 1 / sample_rate;


frame_count = 1 + floor( ( N - frame_length ) / frame_hop );

indices = zeros( frame_count, 2 );

frame_set = nan( frame_count, frame_length );

for index = 0:1:( frame_count - 1 )
    indices( index + 1, 1 ) = index * frame_hop + 1;
    indices( index + 1, 2 ) = index * frame_hop + frame_length;

    frame_set( index + 1, : ) = signal( indices( index + 1, 1 ):1:indices( index + 1, 2 ) );
end


frames.signal = signal;
frames.number_of_samples = N;
frames.sample_rate = sample_rate;
frames.sample_period = sample_period;
frames.count = frame_count;
frames.indices = indices;
frames.frame_set = frame_set;
frames.rate = sample_rate / frame_hop;
frames.period = 1 / frames.rate;



%% Reference(s)

% Framing:  https://brianmcfee.net/dstbook-site/content/ch09-stft/Framing.html

% See Coding Interlude 4.


