


%% Synopsis

% ACS 547, Lecture 7 - Example 1, Sabine Rooms

% See slide 10 on "Lecture 07 - Sabine rooms - Filled.pptx".



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

PRINT_FIGURES = 0;



%% Information

room.width = 10;  room.length = 10;  room.height = 5;  % meters

frequency = 250;  % Hz

alpha.pre_treatment = 0.08;  % m^2 or Sabines
alpha.treatment = 0.5;  % m^2 or Sabines

directivity_factor = 2;  % Machine located in the center of the floor.



%% Calculate the Old Room Constant

S = 2*( room.width * room.length ) + 2*(room.height * room.length )  +  2*(room.height * room.width);  % 400 m^2

room_constant_old = S * alpha.pre_treatment / ( 1 - alpha.pre_treatment );  % 34.8 m^2 or Sabines



%% Calculate the New Room Constant

alpha_new = ( ( (room.length * room.width) + (4 * room.length * room.height ) )*alpha.treatment  +  (room.length * room.width)*alpha.pre_treatment ) / S;  % 0.395
%
% Room area does not change.

room_constant_new = S * alpha_new / ( 1 - alpha_new );  % 261.2 m^2 or Sabines


room_constant_new / room_constant_old;  % 7.5 times more absorption.



%% Pressure Difference

% In reverberant field (no term with depedence on distance from source).

pressure_difference = 10*log10(  room_constant_old / room_constant_new );  % -8.8 dB



%% Plot of Level Difference Versus Distance from Source

r = 0:0.1:10;  % meters

h_delta_Lp = @( D0, r, Rnew, Rold )  10*log10( D0./(4.*pi.*r.^2) + 4 / Rnew )  -  10*log10( D0./(4.*pi.*r.^2) + 4 / Rold );


figure( ); ...
    plot( r, h_delta_Lp( 2, r, room_constant_new, room_constant_old ) );  grid on;
    xlabel( 'Distance from Source [meters]' );  ylabel( 'Level Difference [dB]' );
    title( 'Pressure Difference for a Point Source Versus Distance from the Source' );
%
% For small distances, there is no benefit of adding new absorption material.
%
% For distance approaching infinity, the level difference approaches -8.8 dB.



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


