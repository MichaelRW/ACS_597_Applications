


%% Synopsis

% Problem 2 - Sabine Room



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.5 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Define Room

room.width = 8;  % m
room.length = 6;  % m
room.height = 3;  % m
    room.volume = room.width * room.length * room.height;  % 144 m^3
    room.area = 2*(room.width*room.height)  +  2*(room.length*room.height)  + 2*(room.width*room.length);  % 180 m^2

alpha_average_walls_and_floor = 0.05;  % For the walls and the floor.
alpha_average_ceiling = 0.15;  % For the ceiling.

% For the 125 Hz octave band.



%% Part a - Estimate the Reverberant Sound Pressure Level

Lw = 10*log10( 25e-3 / 1e-12 );  % 103.98 dB

average_absorption_coefficient = ( (room.width*room.length)*alpha_average_ceiling  +  (room.width*room.length  +  2*(room.width*room.height)  + 2*(room.length*room.height) )*alpha_average_walls_and_floor ) / room.area;  % 0.076667 unitless

room_constant = room.area * average_absorption_coefficient / ( 1 - average_absorption_coefficient );  % 14.9 m^2 or Sabin


sound_pressure_level = Lw + 10*log10( 4 / room_constant );  % 98.3 dB



%% Part b - Calculate the Critical Distance

D0 = 1;

r = 0:0.05:12;  % m

h_Lp_direct = @( Lw, D0, r )  Lw + 10*log10( D0./(4.*pi.*r.^2) );
h_Lp_reverberant = @( Lw, room_constant )  Lw + 10*log10( 4 ./ room_constant );
%
h_Lp_net = @( Lw, D0, r, room_constant )  Lw + 10*log10( D0./(4.*pi.*r.^2) + 4/room_constant ) + 10*log10( 343*1.2/400 );


figure( ); ...
    h1 = plot( r, h_Lp_direct( Lw, D0, r ) );  hold on;
    h2 = plot( r, ones( size(r) ).*h_Lp_reverberant( Lw, room_constant) );
    h3 = plot( r, h_Lp_net( Lw, D0, r, room.volume ));  grid on;
        legend( [ h1, h2 h3 ], { 'Direct $L_p$', 'Reverberant $L_p$', 'Total $L_p$' }, 'Interpreter', 'Latex', 'Location', 'SouthWest' );        
    %
    text( 0.56, 105, 'Critical Distance $\approx$ 0.55 meters.', 'Interpreter', 'Latex' );
        line( [ 0.545 0.545 ], [ 90 110 ], 'Color', 'k', 'LineWidth', 0.6 );
    %
    xlabel( 'Distance from Source [meters]' );  ylabel( 'Sound Pressure Level [dB re:20e-6 Pa]' );
    title( 'Sound Pressure Levels of 125 Hz Point Source' );
    %
    set( gca, 'XScale', 'log' );
    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.3 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Sabine Room', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex


% Estimate the critical distance (see page 84 of "06-Indoors.pdf" notes for ACS 537).
rc = 0.141 * sqrt( D0 * room_constant );  % 0.5451 meters



%% Clean-up


fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


