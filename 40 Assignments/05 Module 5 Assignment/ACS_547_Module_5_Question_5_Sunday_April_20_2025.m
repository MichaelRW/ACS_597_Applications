


%% Synopsis

% Problem 5 - Order Tracking of Mobius Propeller Data

% See Lecture 26
%   "D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\26 Monday, April 21, 2025\Lecture 26 - Shafting and bearings.pptx"



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

% set( groot, 'DefaultFigurePosition', [ 230 730  3*750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.6 );
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



%% Visualize Experimental Data

net_time = size( p, 1 ) / fs;  % 45 seconds

time_indices = ( 0:1:( numel( trigger ) - 1 ) ) ./ fs;
    time_indices = time_indices(:);


figure( 'Name', 'Trigger signal and Channel 1 Pressure Data' ); ...

    h1 = subplot( 2, 1, 1 ); ...
        plot( time_indices, trigger, 'Color', color_map( 1, : ) );  grid on;
            xlabel( 'Time [s]' );  ylabel( 'Voltage [V]' );  title( 'Shaft Trigger Signal' );
            ylim( [ -0.5 4 ] );
    h2 = subplot( 2, 1, 2 ); ...
        plot( time_indices, p( :, 1 ), 'Color', color_map( 2, : ) );  grid on;
            xlabel( 'Time [s]' );  ylabel( 'Pressure [Pa]' );  title( 'Recorded Pressure at 23$^\circ$ Elevation' );
            ylim( [ -5 5.5 ] );

    linkaxes( [ h1 h2 ], 'x' );
        xlim( [ -5 50 ] );



%% Spectrogram

signal = p(:, 1);

frame_length = 8192;
    frame_hop = floor( 0.20 * frame_length );

a_spectrogram = spectrogram_September_26_2023( signal, frame_length, frame_hop, fs, 0 );

figure( 'Name', 'Spectrogram of Channel 1 Pressure Signal' ); ...

    imagesc( a_spectrogram.time_indices, a_spectrogram.Sxx.frequencies, 20*log10( sqrt(a_spectrogram.Sxx.spectrum)./20e-6 ), [ 0 85 ] );
        labelColorbar( 'Power Spectral Density [dB re: $\frac{20 \mu Pa^2}{Hz}$]' );  grid on;
        colormap parula;  % Option:  Turbo
        set( gca, 'GridColor', 'w', 'GridAlpha', 0.4 );
    xlabel( 'Time [s]' );  ylabel( 'Frequency [Hz]' );
    set( gca, 'ydir', 'normal' );
    ylim( [ 0 1e3] );



%% Locate the Leading Edges of the Trigger Data

% Smooth data to reduce noise;  supports finding the leading eges.
smoothed_trigger_data = [ 0;  diff( smooth( trigger, 'moving', 5 ) ) ];

% Threshold smoothed data to located leading edges.
threshold = 0.5;
    edge_indices = find( threshold  <  smoothed_trigger_data );
        leading_edges = edge_indices( 1:2:end );

temp = diff( leading_edges );
    [ indices_to_remove, indices ] = find( temp == 2 );
        clear temp;

temp = leading_edges;
    temp( indices_to_remove ) = [ ];

leading_edges = temp;
    clear temp;

time_indices_leading_edge = leading_edges ./ fs;    


% figure( 'Name', 'Revolution Segmentation' ); ...
% 
%     plot( trigger( 1:1:end ), 'Color', 'k' );  hold on;
% 
%     for index = 1:1:numel( leading_edges )
%         line( [ leading_edges( index ) leading_edges( index ) ], [ 0 4 ], 'Color', 'b' );  grid on;
%     end
% 
%     axis( [ 1 numel( trigger )  -0.5  4.5 ] );



%% Segment Microphone Recording Channels

% Channel 1:  23 degrees
% Channel 2:  18 degrees
% Channel 3:  12 degrees
% Channel 4:  6 degrees
% Channel 5:  0 degrees
% Channel 6:  -6 degrees
% Channel 7:  -12 degrees
% Channel 8:  -18 degrees
% Channel 9:  -23 degrees
% Channel 10:  -28 degrees
% Channal 11:  -32 degrees


start_indices = 1:1:( numel( leading_edges ) - 1 );  end_indices = 2:1:numel( leading_edges );

frame_set_indices = [ leading_edges( start_indices(:) )  leading_edges( end_indices(:) ) ];
    number_of_revolutions = size( frame_set_indices, 1 );


for channel_index = 1:1:size( p, 2 )

    for frame_index = 1:1:size( frame_set_indices, 1 )
    
        channel_segments{ frame_index, channel_index } = p( frame_set_indices( frame_index, 1 ):1:frame_set_indices( frame_index, 2 ), channel_index );
    
        slope_theta = (2*pi) ./ size( channel_segments{ frame_index, channel_index }, 1 );
            theta{ frame_index, channel_index } = slope_theta .* ( 0:1:size( channel_segments{ frame_index, channel_index }, 1 ) - 1 );
    
    end

end



%% Compute Instantaneous RPM

for channel_index = 1:1:size( p, 2 )
    for frame_index = 1:1:size( frame_set_indices, 1 )    
        rpm_estimate{ frame_index, channel_index } = 60 * ( 1 / ( numel( channel_segments{ frame_index, channel_index } ) ./ fs ) );    
    end
end


figure( 'Name', 'Shaft Speed' ); ...

    yyaxis left; ...
        plot( [ rpm_estimate{ :, 1 } ] ./ 60 );  grid on;
        ylabel( 'Cycles Per Second [CPS]' );
        ylim( [ 0 65 ] );

    yyaxis right; ...
        plot( [ rpm_estimate{ :, 1 } ] );  grid on;
        ylabel( 'Rotations Per Minute [RPM]' );

    xlabel( 'Revolution [WU]' );

    axis( [ -50 number_of_revolutions+50  0 3.75e3 ] );
    shg;



%% Compute Fourier Coefficients

An = nan( size( frame_set_indices, 1 ), 11 );  Bn = An;

TRACKING_ORDER = 2;
% TRACKING_ORDER = 4;

for channel_index = 1:1:size( p, 2 )
    
    for frame_index = 1:1:size( frame_set_indices, 1 )        

        M = numel( channel_segments{ frame_index, channel_index } );
            
        An( frame_index, channel_index ) = ...
            2/M * channel_segments{ frame_index, channel_index }.' * cos( TRACKING_ORDER/2 .* theta{ frame_index, channel_index } ).';

        Bn( frame_index, channel_index ) = ...
            2/M * channel_segments{ frame_index, channel_index }.' * sin( TRACKING_ORDER/2 .* theta{ frame_index, channel_index } ).';

    end

end



%% Plot Fourier Coefficients An and Bn and the Overall Sound Pressure Level for Channel 1

if ( TRACKING_ORDER == 2 )

    figure( 'Name', 'Second-order Tracking and Associated Sound Pressure Level for Channel 1' ); ...
    
        subplot( 2, 1, 1 ); ...
            plot( An( :, 1 ), 'Color', 'b' );  hold on;
            plot( Bn( :, 1 ), 'Color', 'r' );  grid on;
                legend( '2x SR (cosine, $A_n$)', '2x SR (sine, $B_n$)', 'Location', 'South', 'Interpreter', 'Latex' );
            xlabel( 'Revolution [WU]' );  ylabel( 'Fourier Coefficient' );
            grid on;  axis( [ 1 number_of_revolutions  -0.3 0.4 ] );
    
        subplot( 2, 1, 2 ); ...
            plot( 10*log10( ( An( :, 1 ).^2 + Bn( :, 1 ).^2 ) ./ 2e-6^2 ), 'Color', 'k' );  grid on;
            xlabel( 'Revolution [WU]' );  ylabel( 'Sound Pressure dB re:20e-6 Pa' );
            axis( [ 1 number_of_revolutions  50 110 ] );
    
        shg;

else

    figure( 'Name', 'Fourth-order Tracking and Associated Sound Pressure Level for Channel 1' ); ...
    
        subplot( 2, 1, 1 ); ...
            plot( An( :, 1 ), 'Color', 'b' );  hold on;
            plot( Bn( :, 1 ), 'Color', 'r' );  grid on;
                legend( '4x SR (cosine, $A_n$)', '4x SR (sine, $B_n$)', 'Location', 'South', 'Interpreter', 'Latex' );
            xlabel( 'Revolution [WU]' );  ylabel( 'Fourier Coefficient' );
            grid on;  axis( [ 1 number_of_revolutions  -1 1 ] );
    
        subplot( 2, 1, 2 ); ...
            plot( 10*log10( ( An( :, 1 ).^2 + Bn( :, 1 ).^2 ) ./ 2e-6^2 ), 'Color', 'k' );  grid on;
            xlabel( 'Revolution [WU]' );  ylabel( 'Sound Pressure dB re:20e-6 Pa' );
            axis( [ 1 number_of_revolutions  50 120 ] );
    
        shg;


end



%% Plot Fourier Coefficients An and Bn and the Overall Sound Pressure Level for all Channels

ORDER = 3;  FRAMELEN = 25;

FOURIER_OFFSET = 0.5;
MAGNITUDE_OFFSET = 25;


if ( TRACKING_ORDER == 2 )

    figure( 'Name', 'Second-order Tracking and Associated Sound Pressure Level for all Channels' ); ...
    
    CHANNEL_INDEX = 11;
    
        h1 = subplot( 2, 2, 1 ); ...
    
            for channel_index = 1:1:CHANNEL_INDEX
    
                if ( channel_index < 5 )
                    h_above = plot( sgolayfilt( An( :, channel_index ), ORDER, FRAMELEN ) + FOURIER_OFFSET, 'Color', 'b', 'LineStyle', '--' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                elseif ( channel_index == 5 )
                    h_plane = plot( sgolayfilt( An( :, channel_index ), ORDER, FRAMELEN ), 'Color', 'b' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                else
                    h_below = plot( sgolayfilt( An( :, channel_index ), ORDER, FRAMELEN ) - FOURIER_OFFSET, 'Color', 'b', 'LineStyle', '-.' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                end
    
            end
            %
            legend( [ h_above h_plane h_below ], { '2x SR (cosine, $A_n$) - Positive Angles', '2x SR (cosine, $A_n$) - Plane', '2x SR (cosine, $A_n$) - Negative Angles' }, 'Location', 'South', 'Interpreter', 'Latex' );
            xlabel( 'Revolution [WU]' );  ylabel( 'Fourier Coefficient' );
            grid on;  axis( [ 1 number_of_revolutions  -2 1 ] );
    
        % fprintf( 1, '\n\n' );
    
        h2 = subplot( 2, 2, 2 ); ...
    
            for channel_index = 1:1:CHANNEL_INDEX
    
                if ( channel_index < 5 )
                    h_above = plot( sgolayfilt( Bn( :, channel_index ), ORDER, FRAMELEN ) + FOURIER_OFFSET, 'Color', 'r', 'LineStyle', '--' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                elseif ( channel_index == 5 )
                    h_plane = plot( sgolayfilt( Bn( :, channel_index ), ORDER, FRAMELEN ), 'Color', 'r' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                else
                    h_below = plot( sgolayfilt( Bn( :, channel_index ), ORDER, FRAMELEN ) - FOURIER_OFFSET, 'Color', 'r', 'LineStyle', '-.' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                end
    
            end
            %
            legend( [ h_above h_plane h_below ], { '2x SR (sine, $B_n$) - Positive Angles', '2x SR (sine, $B_n$) - Plane', '2x SR (sine, $B_n$) - Negative Angles' }, 'Location', 'South', 'Interpreter', 'Latex' );
            xlabel( 'Revolution [WU]' );  ylabel( 'Fourier Coefficient' );
            grid on;  axis( [ 1 number_of_revolutions  -2 1 ] );
    
        % fprintf( 1, '\n\n' );
    
        h3 = subplot( 2, 2, [ 3 4 ] ); ...
    
            for channel_index = 1:1:CHANNEL_INDEX
    
                if ( channel_index < 5 )
                    h_above = plot( sgolayfilt( 10*log10( ( An( :, channel_index ).^2 + Bn( :, channel_index ).^2 ) ./ 2e-6^2 ) + MAGNITUDE_OFFSET, ORDER, FRAMELEN ), 'Color', 'k', 'LineStyle', '--' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                elseif ( channel_index == 5 )
                    h_plane = plot( sgolayfilt( 10*log10( ( An( :, channel_index ).^2 + Bn( :, channel_index ).^2 ) ./ 2e-6^2 ), ORDER, FRAMELEN ), 'Color', 'k' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                else
                    h_below = plot( sgolayfilt( 10*log10( ( An( :, channel_index ).^2 + Bn( :, channel_index ).^2 ) ./ 2e-6^2 ) - MAGNITUDE_OFFSET, ORDER, FRAMELEN ), 'Color', 'k', 'LineStyle', '-.' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                end
    
            end
            %
            legend( [ h_above h_plane h_below ], { 'Magnitude Response - Positive Angles', 'Magnitude Response - Plane', 'Magnitude Response - Negative Angles' }, 'Location', 'South', 'Interpreter', 'Latex' );
            xlabel( 'Revolution [WU]' );  ylabel( 'Sound Pressure Level [dB re: 20~$\mu$ Pa]' );
            grid on;  axis( [ 1 number_of_revolutions  20 140 ] );
    
        shg;

else

    FOURIER_OFFSET = 1;

    figure( 'Name', 'Fourth-order Tracking and Associated Sound Pressure Level for all Channels' ); ...
    
    CHANNEL_INDEX = 11;
    
        h1 = subplot( 2, 2, 1 ); ...
    
            for channel_index = 1:1:CHANNEL_INDEX
    
                if ( channel_index < 5 )
                    h_above = plot( sgolayfilt( An( :, channel_index ), ORDER, FRAMELEN ) + FOURIER_OFFSET, 'Color', 'b', 'LineStyle', '--' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                elseif ( channel_index == 5 )
                    h_plane = plot( sgolayfilt( An( :, channel_index ), ORDER, FRAMELEN ), 'Color', 'b' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                else
                    h_below = plot( sgolayfilt( An( :, channel_index ), ORDER, FRAMELEN ) - FOURIER_OFFSET, 'Color', 'b', 'LineStyle', '-.' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                end
    
            end
            %
            legend( [ h_above h_plane h_below ], { '4x SR (cosine, $A_n$) - Positive Angles', '4x SR (cosine, $A_n$) - Plane', '4x SR (cosine, $A_n$) - Negative Angles' }, 'Location', 'South', 'Interpreter', 'Latex' );
            xlabel( 'Revolution [WU]' );  ylabel( 'Fourier Coefficient' );
            grid on;  axis( [ 1 number_of_revolutions  -4 2 ] );
    
        % fprintf( 1, '\n\n' );
    
        h2 = subplot( 2, 2, 2 ); ...
    
            for channel_index = 1:1:CHANNEL_INDEX
    
                if ( channel_index < 5 )
                    h_above = plot( sgolayfilt( Bn( :, channel_index ), ORDER, FRAMELEN ) + FOURIER_OFFSET, 'Color', 'r', 'LineStyle', '--' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                elseif ( channel_index == 5 )
                    h_plane = plot( sgolayfilt( Bn( :, channel_index ), ORDER, FRAMELEN ), 'Color', 'r' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                else
                    h_below = plot( sgolayfilt( Bn( :, channel_index ), ORDER, FRAMELEN ) - FOURIER_OFFSET, 'Color', 'r', 'LineStyle', '-.' );  hold on;
                        % fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                end
    
            end
            %
            legend( [ h_above h_plane h_below ], { '4x SR (sine, $B_n$) - Positive Angles', '4x SR (sine, $B_n$) - Plane', '4x SR (sine, $B_n$) - Negative Angles' }, 'Location', 'South', 'Interpreter', 'Latex' );
            xlabel( 'Revolution [WU]' );  ylabel( 'Fourier Coefficient' );
            grid on;  axis( [ 1 number_of_revolutions  -4 2 ] );
    
        fprintf( 1, '\n\n' );
    
        h3 = subplot( 2, 2, [ 3 4 ] ); ...
    
            for channel_index = 1:1:CHANNEL_INDEX
    
                if ( channel_index < 5 )
                    h_above = plot( sgolayfilt( 10*log10( ( An( :, channel_index ).^2 + Bn( :, channel_index ).^2 ) ./ 2e-6^2 ) + MAGNITUDE_OFFSET, ORDER, FRAMELEN ), 'Color', 'k', 'LineStyle', '--' );  hold on;
                        fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                elseif ( channel_index == 5 )
                    h_plane = plot( sgolayfilt( 10*log10( ( An( :, channel_index ).^2 + Bn( :, channel_index ).^2 ) ./ 2e-6^2 ), ORDER, FRAMELEN ), 'Color', 'k' );  hold on;
                        fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                else
                    h_below = plot( sgolayfilt( 10*log10( ( An( :, channel_index ).^2 + Bn( :, channel_index ).^2 ) ./ 2e-6^2 ) - MAGNITUDE_OFFSET, ORDER, FRAMELEN ), 'Color', 'k', 'LineStyle', '-.' );  hold on;
                        fprintf( 1, '\n%d - %3.1f', channel_index, mic_angles_degrees( channel_index ) );
                end
    
            end
            %
            legend( [ h_above h_plane h_below ], { 'Magnitude Response - Positive Angles', 'Magnitude Response - Plane', 'Magnitude Response - Negative Angles' }, 'Location', 'South', 'Interpreter', 'Latex' );
            xlabel( 'Revolution [WU]' );  ylabel( 'Sound Pressure Level [dB re: 20~$\mu$ Pa]' );
            grid on;  axis( [ 1 number_of_revolutions  20 150 ] );
    
        shg;

end



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

% https://www.mathworks.com/help/signal/ref/ordertrack.html#bvaw5ra-1-x

% https://www.mathworks.com/matlabcentral/answers/1439884-evaluate-rising-edge-sample-of-a-signal


