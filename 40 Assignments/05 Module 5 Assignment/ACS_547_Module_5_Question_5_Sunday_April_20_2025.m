


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

pause( 1 );



%% Constants and Parameters


% color_map = slanCM( 'ColorBlind', 4 );



%% Load Data

load( 'mobius_prop_data.mat' );  % Variable(s):  fs;  mic_angles_degrees;  p;  trigger



%% Visualize Data

% close all;  clc;  pause( 1 );
% 
% fs;  % sample rate is 80 kHz
% figure;  plot( mic_angles_degrees );  % 23, 18, 12, 6, 0, -6, -12, -18, -23, -28, -32
% figure;  plot( trigger );

% return

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



for frame_indices = 1:1:size( frame_set_indices, 1 )
    segments{ frame_indices } = p( frame_set_indices( frame_indices, 1 ):1:frame_set_indices( frame_indices, 2 ) );
    working_theta{ frame_indices } = (2*pi) ./ ( ( frame_set_indices( frame_indices, 2 ) - frame_set_indices( frame_indices, 1 ) + 1 ):-1:1 );
end

segments = segments(:);
working_theta = working_theta(:);



%% Compute Fourier Coefficients

close all;

number_of_revolutions = size( segments, 1 );

An = nan( number_of_revolutions, 1 );  Bn = nan( number_of_revolutions, 1 );

TRACKING_ORDER = 100;


for revolution_index = 1:1:number_of_revolutions

    M = numel( segments{ revolution_index } );
        n = 0:1:( M - 1 );


    temp9 = segments{ revolution_index };
    temp10 = working_theta{ revolution_index };  temp10 = temp10(:);

    temp12 = cos( TRACKING_ORDER.*temp10 );
    temp13 = sin( TRACKING_ORDER.*temp10 );

    temp11 = ( 2/M .* temp9 ) * temp12;
    temp14 = ( 2/M .* temp9 ) * temp13;

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

%% Hold

% signal = p(:, 1);
% % signal = p(:, 5);
% 
% frame_length = 2048;  frame_hop = 1024;
% 
% APPLY_HANN_WINDOW = 1;
%     spectrogram_mrw = spectrogram_September_26_2023( signal, frame_length, frame_hop, fs, APPLY_HANN_WINDOW );
% 
% figure( ); ...
%     h_1 = imagesc( spectrogram_mrw.time_indices, spectrogram_mrw.Sxx.frequencies, 10*log10( spectrogram_mrw.Sxx.spectrum ), [ -100 -10 ] );
%         labelColorbar( 'PSD [dB re: 1 $\frac{Volts^2}{Hz}$]' );  grid on;
%         colormap parula;  % Option:  Turbo
%         set( gca, 'GridColor', 'w', 'GridAlpha', 0.3 );
%     title( 'Order Tracking', 'Interpreter', 'Latex' );
%     xlabel( 'Time [seconds]', 'Interpreter', 'Latex' );  ylabel( 'Frequency [Hz]', 'Interpreter', 'Latex' );
%     set( gca, 'ydir', 'normal' );
%     ylim( [ 0 fs/2] );
%     % datatip( h_1, 2, +128 );
%     % datatip( h_1, 2, -128, 'Location', 'southeast' );


    
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


