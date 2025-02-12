


%% Synopsis

% Slide 8 - Noise Reduction and Transmission Loss



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



%% Define Room

room.length = 20;  % meters
room.width = 10;  % meters
room.height = 4;  % meters
    room.volume = room.length * room.width * room.height;  % m^2
    room.area = 2*(room.length * room.width) + 2*(room.length * room.height) + 2*(room.width * room.height);  % m^2



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


