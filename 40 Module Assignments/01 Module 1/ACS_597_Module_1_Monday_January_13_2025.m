

%% Synopsis

% Homework Set 9, Fireworks

% Question 1 of 5, Familiarizing Ourselves with the Data



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( fullfile( '..', '..', '05 Support' ), '-begin' );

set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.0 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Fireworks Time Series - Information

audioinfo( 'Test_shoot_A.wav' );

% Filename:  '.\Test_shoot_A.wav'
% CompressionMethod:  'Uncompressed'
% NumChannels:  4
% SampleRate:  48,000
% TotalSamples:  2,064,000
% Duration:  43 seconds
% Title:  []
% Comment:  []
% Artist:  []
% BitsPerSample:  16



%% Fireworks Time Series

[ s, fs ] = audioread( 'Test_shoot_A.wav' );
%
assert( fs == 48e3, '*** Sample rate is not 48 kHz. ***' );

sample_period_seconds = 1 / fs;
    time_indices = ( 0:1:( size( s, 1 ) - 1 ) ) .* sample_period_seconds;

figure( ); ...
    plot( time_indices, s( :, 1 ) );  hold on;  % Channel 1
    plot( time_indices, s( :, 2 ) );  % Channel 2
    plot( time_indices, s( :, 3 ) );  % Channel 3
    plot( time_indices, s( :, 4 ) );  grid on;  % Channel 4
    legend( 'Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Location', 'NorthWest' );
    xlabel( 'Time [seconds]' );  ylabel( 'Amplitude [WU]' );  title( 'Fireworks Recording' );
    xlim( [ time_indices( [ 1 end ] ) ] );  ylim( [ -0.15 +0.34 ] );
    %
    if ( PRINT_FIGURES == 1 ), print( fullfile( '.', 'HS9_Q1_Time_Series.png' ), '-dpng' );  end
    


%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


