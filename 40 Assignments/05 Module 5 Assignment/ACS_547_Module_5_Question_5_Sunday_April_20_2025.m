


%% Synopsis

% Problem 5 - Order Tracking of Mobius Propeller Data

% See Lecture 26
%   ""D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\26 Monday, April 21, 2025\Lecture 26 - Shafting and bearings.pptx""



%% Note(s)



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

% set( groot, 'DefaultFigurePosition', [ 1.7e3  775    750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'norma' );
set( 0, 'DefaultLineLineWidth', 1.0 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Constants and Parameters


color_map = slanCM( 'ColorBlind', 4 );



%% Load Data

load( 'mobius_prop_data.mat' );  % Variable(s):  fs;  mic_angles_degrees;  p;  trigger



%% Visualize Data

% for index = 1:1:11
%     figure( );  plot( p( :, index ) );  grid on;
% end



%% Listen

% soundsc( p( :, 1 ), fs );  % 23 degrees

% i = 5;  figure( );  plot( p( :, i ) );  grid on;  soundsc( p( :, i ), fs );  % 0 degrees

% soundsc( p( :, 11 ), fs );  % -32 degrees



signal = p(:, 1);
% signal = p(:, 5);

frame_length = 2048;  frame_hop = 1024;

APPLY_HANN_WINDOW = 1;
    spectrogram_mrw = spectrogram_September_26_2023( signal, frame_length, frame_hop, fs, APPLY_HANN_WINDOW );

figure( ); ...
    h_1 = imagesc( spectrogram_mrw.time_indices, spectrogram_mrw.Sxx.frequencies, 10*log10( spectrogram_mrw.Sxx.spectrum ), [ -100 -10 ] );
        labelColorbar( 'PSD [dB re: 1 $\frac{Volts^2}{Hz}$]' );  grid on;
        colormap parula;  % Option:  Turbo
        set( gca, 'GridColor', 'w', 'GridAlpha', 0.3 );
    title( 'Order Tracking', 'Interpreter', 'Latex' );
    xlabel( 'Time [seconds]', 'Interpreter', 'Latex' );  ylabel( 'Frequency [Hz]', 'Interpreter', 'Latex' );
    set( gca, 'ydir', 'normal' );
    ylim( [ 0 fs/2] );
    % datatip( h_1, 2, +128 );
    % datatip( h_1, 2, -128, 'Location', 'southeast' );
    

    
%% Clean-up

if ( ~isempty( findobj( 'Type', 'figure' ) ) )
    monitors = get( 0, 'MonitorPositions' );
        if ( size( monitors, 1 ) == 1 )
            autoArrangeFigures( 3, 4, 1 );
        elseif ( 1 < size( monitors, 1 ) )
            autoArrangeFigures( 2, 2, 2 );
        end
end

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)

% p = 10 + 1i*5;  % Pascals;  sinusoid
% 
% p_mag = abs( p );  % 11.18 Pascals
% 
% p_rms = p_mag / sqrt(2);  % 7.91 Pascals RMS
% 
% p_dB_SPL = 20*log10( p_rms / 20e-6 );  % 111.94 dB SPL Z
% 
% 
% p_dB_SPL_verify = convert_complex_pressure_to_dB_SPL( p );




%% URLs


