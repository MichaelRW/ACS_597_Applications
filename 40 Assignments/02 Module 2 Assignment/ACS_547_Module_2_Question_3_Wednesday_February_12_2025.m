


%% Synopsis

% Question 3 - Transmission Loss Measurement


% Reference(s):
%
% Slide 8 - Noise Reduction and Transmission Loss


% Assumptions:
%
%   1.)  Lp1 depends on the transmission loss.
%
%       If the transmission loss is low then more energy goes to room 2 (i.e., the receiver room).
%
%       The noise reduction from the source room to the receiver.
%
%       Adding the barrier will change the level in the source room.  Typically making the sound level higher in the source room.



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.0 );
set( 0, 'DefaultLineMarker', 'x' );
set( 0, 'DefaultLineMarkerSize', 15 );
% set( 0, 'DefaultAxesLineStyleOrder', { '-' '--o' '+' } );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Dimensions of Rooms and Panel

room.length = 4;  % meters
room.width = 4;  % meters
room.height = 4;  % meters
    room.volume = room.length * room.width * room.height;  % 64 m^2
    room.area = 2*(room.length * room.width) + 2*(room.length * room.height) + 2*(room.width * room.height);  % 96 m^2

panel.width = 0.8;  % meters
panel.height = 0.8;  % meters
    panel.area = panel.width * panel.height;  % 0.64 m^2



%% Data

octave_band_frequencies = [ 125  250  500  1000  2000  4000 ].';  % Hz
T60 = [ 2.0  2.1  1.8  1.5  1.2  0.9 ].';  % seconds
spl.source_room = [ 90  95  103  105  100  93 ].';  % dB re: 20e-6 Pascals
spl.receiver_room = [ 50  50  46  50  50  38 ].';  % dB re: 20e-6 Pascals

c = 343;  % meters per second



%% Pressure Difference

spl.delta = spl.source_room - spl.receiver_room;

% figure( ); ...
%     stem( octave_band_frequencies, spl.delta, 'Marker', '.', 'MarkerSize', 12, 'Color', 'r' );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Transmission Loss [dB]' );
%     title( 'Pressure Difference Versus Octave Band Center Frequency' );
%     %
%     xticks( octave_band_frequencies );  xticklabels( num2cell( octave_band_frequencies ) );
%     set( gca, 'XScale', 'log' );
%     xlim( [ 80 6e3 ] );  ylim( [ 0 65 ] );

% https://www.mathworks.com/matlabcentral/answers/413686-how-to-set-log-scale-range



%% Determine Average Absorption in the Receiver Room using Reverberation Time Measurements

average_absorption = @( volume, area, c, T60 )  ( 55.25 .* volume ) ./ ( area .* c .* T60 );

receiver_room.average_absorption = average_absorption( room.volume, room.area, c, T60 );

% Assumption:  Calibration panel has very high transmission loss.

% figure( ); ...
%     stem( octave_band_frequencies, receiver_room.average_absorption, 'Marker', '.', 'MarkerSize', 12, 'Color', 'r' );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Average Absorption [Sabine]' );
%     title( 'Average Absorption in Receiver Room Versus Octave Band Center Frequency' );
%     %
%     xticks( octave_band_frequencies );  xticklabels( num2cell( octave_band_frequencies ) );
%     set( gca, 'XScale', 'log' );
%     xlim( [ 80 6e3 ] );  ylim( [ 0 0.14 ] );



%% Determine Receiver Room Constant

room_constant = @( average_absorption, area )  ( average_absorption * area ) ./ ( 1 - average_absorption );  % Unitless

receiver_room.room_constant = room_constant( receiver_room.average_absorption, room.area );



%% Determine the Transmission Loss in Each Octave Band

% The calibration plate isolate the receiver room and the area of the
% receiver room does not consider the calibration plate.

transmission_coefficient = @( receiver_room_pressure, source_room_pressure, panel_area, receiver_room_constant )  ( ( receiver_room_pressure ./ source_room_pressure ) .* receiver_room_constant ) ./ panel_area;

tau = transmission_coefficient( 10.^(spl.receiver_room./10)*20e-6, 10.^(spl.source_room/10)*20e-6, panel.area, receiver_room.room_constant );

TL = -10*log10( tau );

figure( ); ...
    stem( octave_band_frequencies, TL, '.', 'MarkerSize', 15, 'Color', 'r' );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Transmission Loss [dB]' );
    title( 'Transmission Loss Per Octave Band' );
    %
    xticks( octave_band_frequencies );  xticklabels( num2cell( octave_band_frequencies ) );
    set( gca, 'XScale', 'log' );
    xlim( [ 80 6e3 ] );  ylim( [ 0 55 ] );



%% Validation

TL_verify = spl.source_room - spl.receiver_room + 10*log10( panel.area ./ receiver_room.room_constant )

return

%% Clean-up

if ( ~isempty( findobj( 'Type', 'figure' ) ) )
    monitors = get( 0, 'MonitorPositions' );
        if ( size( monitors, 1 ) == 1 )
            autoArrangeFigures( 2, 2, 1 );
        elseif ( 1 < size( monitors, 1 ) )
            autoArrangeFigures( 2, 2, 1 );
        end
end


fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


