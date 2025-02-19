


%% Synopsis

% Question 3 - Transmission Loss Measurement


% Reference(s):
%
% Slide 8 - Noise Reduction and Transmission Loss


% Assumptions:
%
%   1.)  



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../40 Assignments/00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
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

receiver_room.room_constant = room_constant( receiver_room.average_absorption, room.area )

% return

%% Determine the Transmission Loss in Each Octave Band

transmission_coefficient = @( receiver_room_pressure, source_room_pressure, panel_area, receiver_room_constant )  ( ( receiver_room_pressure ./ source_room_pressure ) .* receiver_room_constant ) ./ panel_area;

tau = transmission_coefficient( 10.^(spl.receiver_room./20)*20e-6, 10.^(spl.source_room/20)*20e-6, panel.area, receiver_room.room_constant );

TL = -10*log10( tau )

figure( ); ...
    stem( octave_band_frequencies, TL, 'Marker', '.', 'MarkerSize', 12, 'Color', 'r' );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Transmission Loss [dB]' );
    title( 'Octave Band Transmission Loss' );
    %
    xticks( octave_band_frequencies );  xticklabels( num2cell( octave_band_frequencies ) );
    set( gca, 'XScale', 'log' );
    xlim( [ 80 6e3 ] );  ylim( [ -10 0 ] );

% return

%% Validation

TL_verify = spl.source_room - spl.receiver_room + 10*log10( panel.area ./ receiver_room.room_constant )

return

%% Source

D = 1;  % Unitless

T60 = 2;  % seconds

SPL_reverberant_field = 80;  % dB



%% Average Room Absorption and Room Constant

alpha_average = (55.25 * room.volume) / (room.area * 343 * T60 );  % 0.1 Sabines

R = (room.area * alpha_average) / (1 - alpha_average);  % 71.6 m^2



%% Solve for the Sound Power Level of the Source

% Lp = Lw + 10*log10( 4 / R );

Lw = 80 - 10*log10( 4 / R );  % 92.5 dB re: 1e-12 Watts
%
% Note(s):
%
%   1.)  Using the reverberant field level of 80 dB (r is large, so its associate term is zero).
%   2.)  The correction term define on slide 3 is not included here.



%% Partition

% Placed at the 20 meter mark of on the room length.

TL = 20;  % dB

transmission_coefficient = 10^( TL / -10 );  % 0.01 unitless



%% Sound Pressures in Each Half of the Room

% Room 1 (with the source)
alpha_new = ...
    ( ( room.area / 2) * 0.1  +  (room.width * room.height * transmission_coefficient ) ) / ...
    (room.area/2 + room.width * room.height);  % 0.09 Sabines

R_new = (room.area/2 * alpha_new) / ( 1 - alpha_new );  % 31.6 m^2

Lp1 = 92.5 + 10*log10( 4 / R_new );  % 83.5 dB


% Room 2
Lp2 = Lp1 - TL + 10*log10( room.width*room.height / 31.6 );  % 64.5 dB


delta_L = Lp1 - Lp2;  % 19 dB



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


