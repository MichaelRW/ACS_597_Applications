


%% Synopsis

% Problem 5 - Order Tracking of Mobius Propeller Data

% See Lecture 26
%   ""D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\26 Monday, April 21, 2025\Lecture 26 - Shafting and bearings.pptx""



%% Note(s)

% Plots of all the data.

% Include labels on plots for each harmonic order (of blade rotation rate).



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

set( groot, 'DefaultFigurePosition', [ 230 730  750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.0 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

color_map = slanCM( 'ColorBlind', 4 );

pause( 1 );



%% Load Data

load( 'mobius_prop_data.mat' );  % Variable(s):  fs;  mic_angles_degrees;  p;  trigger

% fs - Sample rate of 80 kHz (scalar).
% mic_angles_degrees - Microphone angles relative to horizontal plane (11-by-1;  23, 18, 6, 0, -6, -12, -18, -23, -28, and -32).
% p - Recorded pressures for the corresponding microphone angles (3,600,000-by-11).
% trigger - Shaft rotation synchronization signal (3,600,000-by-11).



%% RMS dB SPL Per Pressure Recording

% p_rms = 20*log10( rms( p ) / 20e-6 );
% 
% figure( 'Name', 'RMS dB SPL Per Pressure Recording' ); ...
%     stem( mic_angles_degrees, fliplr( p_rms ) );  grid on;
%         xlabel( 'Microphone Angle [$^\circ$]' );
%         xticks( fliplr( mic_angles_degrees ) );
%         ylabel( 'Sound Pressure [dB SPL]' );  title( 'Sound Pressure Level Versus Microphone Angle' );



%% Visualize Experimental Data

% net_time = size( p, 1 ) / fs;  % 45 seconds


% time_indices = ( 0:1:( numel( trigger ) - 1 ) ) ./ fs;
% 
% figure( 'Name', 'Trigger signal and Data for Horizontal Plane' ); ...
% 
%     h1 = subplot( 2, 1, 1 ); ...
%         plot( time_indices, trigger, 'Color', color_map( 1, : ) );  grid on;
%             xlabel( 'Time [s]' );  ylabel( 'Amplitude [WU]' );  title( 'Shaft Trigger Signal' );
%             ylim( [ -0.5 4 ] );
%     h2 = subplot( 2, 1, 2 ); ...
%         plot( time_indices, p( :, 5 ), 'Color', color_map( 2, : ) );  grid on;
%             xlabel( 'Time [s]' );  ylabel( 'Pressure [Pa]' );  title( 'Recorded Pressure at 0$^\circ$ Elevation' );
%             ylim( [ -5 5.5 ] );
% 
%     linkaxes( [ h1 h2 ], 'x' );
%         xlim( [ -5 50 ] );



%% Spectrogram

% signal = p(:, 1);
% 
% frame_length = 8192;
%     frame_hop = floor( 0.20 * frame_length );
% 
% a_spectrogram = spectrogram_September_26_2023( signal, frame_length, frame_hop, fs, 0 );
% 
% figure( 'Name', 'Spectrogram of Pressure Signal' ); ...
% 
%     imagesc( a_spectrogram.time_indices, a_spectrogram.Sxx.frequencies, 20*log10( sqrt(a_spectrogram.Sxx.spectrum)./20e-6 ), [ 0 85 ] );
%         labelColorbar( 'Power Spectral Density [dB re: $\frac{20 \mu Pa^2}{Hz}$]' );  grid on;
%         colormap parula;  % Option:  Turbo
%         set( gca, 'GridColor', 'w', 'GridAlpha', 0.4 );
%     xlabel( 'Time [s]' );  ylabel( 'Frequency [Hz]' );  title( 'Order Tracking' );
%     set( gca, 'ydir', 'normal' );
%     ylim( [ 0 1e3] );



%% RESAMPLE - time-domain sampling to angular domain

threshold = 0.5;

start_angle = 0;  % radians
end_angle = 2*pi; % radians
    y_difference = end_angle - start_angle;


% Locate the leading edges of the trigger data.
smoothed_trigger_data = [ 0;  diff( smooth( trigger, 'moving', 5 ) ) ];
    edge_indices = find( threshold  <  smoothed_trigger_data );
        lead_edge = edge_indices( 1:2:end );


lead_edge( 1:10 );

temp = diff( lead_edge );
    temp( 1:10 );

[ indices_to_remove, indices ] = find( temp == 2 );    

temp2 = lead_edge;
    temp2( indices_to_remove ) = [ ];

temp2;    

% figure;  plot( temp2 );




theta = [ ]; 

start_index = [ ];  end_index = [ ];

% Linear interpolation for each interval between the leading edges.
for index = 1:1:numel( lead_edge )


    % fprintf( 1, '\n%d of %d', index, numel( lead_edge ) );


    if ( index == 1 )
        x_difference = lead_edge( index );
    else
        x_difference = lead_edge( index ) - lead_edge( index - 1 );
            start_index = [ start_index  lead_edge( index - 1 ) ];
            end_index = [ end_index  lead_edge( index )  ];
    end    
    
    slope = y_difference / x_difference;

    segment = (0:1:x_difference)*slope;
        theta = [ theta segment ];

    % keyboard
    
end


    

% figure( ); ...
%     h1 = subplot( 2, 1, 1 ); ...
%         plot( trigger( 1:1:end ) );  hold on;
%         %
%         for index = 1:1:numel(temp2)
%             line( [ temp2(index) temp2(index)], [ 0 4 ], 'Color', 'r' );  grid on;
%         end
% 
%     h2 = subplot( 2, 1, 2 ); ...
%         plot( theta );  grid on;
%     %
%     linkaxes( [ h1 h2 ], 'xy' );
%         axis( [ 1 5e5   0 6.5 ] );
%         % axis( [ 1 numel( trigger )   0 6.5 ] );

% return

%% Compute Fourier Coefficients Directly for Each Revolution

% block data for each revolution;  used block data directly



%% Listen

% soundsc( p( :, 1 ), fs );  % 23 degrees

% i = 5;  figure( );  plot( p( :, i ) );  grid on;  soundsc( p( :, i ), fs );  % 0 degrees

% soundsc( p( :, 11 ), fs );  % -32 degrees


temp2;  % 2,296-by-1 elements

% indices = reshape( temp2, 2, numel(temp2)/2 );
%     indices( :, 1:4 )

indices = temp2;
    indices( 1:6 )

start_indices = 1:1:( numel( indices ) - 1 );
end_indices = 2:1:numel( indices );

frame_set_indices = [ indices( start_indices(:) )  indices( end_indices(:) ) ];
    frame_set_indices( 1:4, : );


WORKING_INDEX = 1;
    p_temp = p( :, WORKING_INDEX );

for frame_indices = 1:1:size( frame_set_indices, 1 )
    segments{ frame_indices } = p_temp( frame_set_indices( frame_indices, 1 ):1:frame_set_indices( frame_indices, 2 ) );
    working_theta{ frame_indices } = (2*pi) ./ ( ( frame_set_indices( frame_indices, 2 ) - frame_set_indices( frame_indices, 1 ) + 1 ):-1:1 );
end

segments = segments(:);
working_theta = working_theta(:);

% return

%% Compute Fourier Coefficients

close all;

number_of_revolutions = size( segments, 1 );

An = nan( number_of_revolutions, 1 );  Bn = nan( number_of_revolutions, 1 );

TRACKING_ORDER = 2;


for revolution_index = 1:1:number_of_revolutions

    M = numel( segments{ revolution_index } );
        n = 0:1:( M - 1 );


    temp9 = segments{ revolution_index };
    temp10 = working_theta{ revolution_index };  temp10 = temp10(:);

    temp12 = cos( TRACKING_ORDER.*temp10 );
    temp13 = sin( TRACKING_ORDER.*temp10 );

    temp11 = ( 2/M .* temp9 ).' * temp12;
    temp14 = ( 2/M .* temp9 ).' * temp13;

        An( revolution_index ) = temp11;
        Bn( revolution_index ) = temp14;

    % keyboard

end

% Sample time is very high, at 80 kHz.

% return

c = sqrt( An.^2 + Bn.^2 );

figure; ...
    plot( An, 'Color', 'b' );  hold on;
    plot( Bn, 'Color', 'r' );  grid on;
        legend( 'Sine', 'Cosine' );

% figure;  plot( 10*log10(c.^2 ./ 20e-6 ) );  grid on;

% return
    
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

% https://www.mathworks.com/matlabcentral/answers/1439884-evaluate-rising-edge-sample-of-a-signal


