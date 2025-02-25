


%% Synopsis

% Question 3 - Transmission Loss Measurement


% Note(s):
%
%   1.)  Lp1 depends on the transmission loss.
%
%       If the transmission loss is low, then more energy goes to room 2 (i.e., the receiver room).
%
%       The noise reduction from the source room to the receiver room.
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
set( 0, 'DefaultLineLineWidth', 1.2 );
set( 0, 'DefaultLineMarker', 'x' );
set( 0, 'DefaultLineMarkerSize', 15 );
% set( 0, 'DefaultAxesLineStyleOrder', { '-' '--o' '+' } );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Dimensions of Rooms and Panel

room.length = 4;  % m
room.width = 4;  % m
room.height = 4;  % m
    room.volume = room.length * room.width * room.height;  % 64 m^3
    room.area = 2*(room.length * room.width) + 2*(room.length * room.height) + 2*(room.width * room.height);  % 96 m^2

panel.width = 0.8;  % m
panel.height = 0.8;  % m
    panel.area = panel.width * panel.height;  % 0.64 m^2

c = 343;  % m/s



%% Measurement Data

octave_band_frequencies = [ 125  250  500  1000  2000  4000 ].';  % Hz
T60 = [ 2.0  2.1  1.8  1.5  1.2  0.9 ].';  % seconds
spl.source_room = [ 90  95  103  105  100  93 ].';  % dB re: 20e-6 Pa
spl.receiver_room = [ 50  50  46  50  50  38 ].';  % dB re: 20e-6 Pa
    % [ octave_band_frequencies  T60  spl.source_room  spl.receiver_room  (spl.source_room - spl.receiver_room) ]



%% Calculate Average Absorption in the Receiver Room using Reverberation Time Measurements

average_absorption = @( volume, area, c, T60 )  ( 55.25 .* volume ) ./ ( area .* c .* T60 );

receiver_room.average_absorption = average_absorption( room.volume, room.area, c, T60 );

% Assumption:  Calibration panel has very high transmission loss.

figure( ); ...
    stem( octave_band_frequencies, receiver_room.average_absorption, 'Marker', '.', 'MarkerSize', 12, 'Color', 'b' );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Average Absorption [$m^2$]' );
    title( 'Average Absorption' );
    %
    xticks( octave_band_frequencies );  xticklabels( num2cell( octave_band_frequencies ) );
    set( gca, 'XScale', 'log' );
    xlim( [ 80 6e3 ] );  ylim( [ 0 0.14 ] );
    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.3 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Q3 Average Absorption', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex



%% Calculate the Receiver Room Constant

% The calibration plate isolates the receiver room.

% The receiver room does not consider the calibration plate.

room_constant = @( average_absorption, area )  ( average_absorption * area ) ./ ( 1 - average_absorption );  % Unitless

receiver_room.room_constant = room_constant( receiver_room.average_absorption, room.area );



%% Calculate the Transmission Loss in Each Octave Band

transmission_coefficient = @( receiver_room_pressure, source_room_pressure, panel_area, receiver_room_constant )  ( ( receiver_room_pressure ./ source_room_pressure ) .* receiver_room_constant ) ./ panel_area;

tau = transmission_coefficient( 10.^(spl.receiver_room./10)*20e-6, 10.^(spl.source_room/10)*20e-6, panel.area, receiver_room.room_constant );

TL = -10*log10( tau );

figure( ); ...
    stem( octave_band_frequencies, TL, '.', 'MarkerSize', 15, 'Color', 'b' );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Transmission Loss [dB]' );
    title( 'Transmission Loss' );
    %
    xticks( octave_band_frequencies );  xticklabels( num2cell( octave_band_frequencies ) );
    set( gca, 'XScale', 'log' );
    xlim( [ 80 6e3 ] );  ylim( [ 0 55 ] );
    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.3 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Q3 Transmission Loss', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex



%% Clean-up


fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


