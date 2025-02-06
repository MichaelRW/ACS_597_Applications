


%% Synopsis

% ACS 547, Modal Behaviour of a Cylindrical Room Milestone

% See slide 32 on "Lecture 06 - Room Modes - Filled.pptx".



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

room.radius = 3;  room.length = 10;  % meters

return

%% ph

room.width = 8;  room.length = 6;  room.height = 3;  % meters
    room_volume = room.width * room.length * room.height;  % 144 m^3
    room_area = 2*(room.width*room.height)  +  2*(room.length*room.height)  + 2*(room.width*room.length);  % 180 m^2

alpha_average_walls_and_floor = 0.05;  % For the walls and the floor.
alpha_average_ceiling = 0.15;  % For the ceiling.
%
% For the 125 Hz octave band.



%% Part a - Estimate the reverberant sound pressure level.

Lw = 10*log10( 25e-3 / 1e-12 );  % 103.98 dB

average_absorption_coefficient = ( (room.width*room.length)*alpha_average_ceiling  +  (room.width*room.length  +  2*(room.width*room.height)  + 2*(room.length*room.height) )*alpha_average_walls_and_floor ) / room_area;  % 0.076667 unitless

room_constant = room_area * average_absorption_coefficient / ( 1 - average_absorption_coefficient );  % 14.9 m^2 or Sabines


sound_pressure_level = Lw + 10*log10( 4 / room_constant );  % 98.3 dB



%% Part b

D0 = 1;

r = 0:0.05:12;  % meters

h_Lp_direct = @( Lw, D0, r )  Lw + 10*log10( D0./(4.*pi.*r.^2) );
h_Lp_reverberant = @( Lw, room_constant )  Lw + 10*log10( 4 ./ room_constant );
%
h_Lp_net = @( Lw, D0, r, room_constant )  Lw + 10*log10( D0./(4.*pi.*r.^2) + 4/room_constant ) + 10*log10( 343*1.2/400 );


figure( ); ...
    plot( r, h_Lp_direct( Lw, D0, r ) );  hold on;
    plot( r, ones( size(r) ).*h_Lp_reverberant( Lw, room_constant) );
    plot( r, h_Lp_net( Lw, D0, r, room_volume ));  grid on;
        legend( 'Direct $L_p$', 'Reverberant $L_p$', 'Total $L_p$', 'Interpreter', 'Latex' );
    %
    text( 0.545, 100, 'Critical Distance $\approx$ 0.55 meters.', 'Interpreter', 'Latex' );
    %
    xlabel( 'Distance from Source [meters]' );  ylabel( 'Sound Pressure Level [dB re:20e-6 Pascals]' );
    title( 'Sound Pressure Components from Direct and Reverberant Fields from 125 Hz Point Source' );
    %
    set( gca, 'XScale', 'log' );


% Estimate the critical distance (see page 84 of "06-Indoors.pdf" notes for ACS 537).
rc = 0.141 * sqrt( D0 * room_constant );  % 0.5451 meters



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

%% Overall, A-weighted Level

% Note(s):
%
%   The above analysis was done using the unweighted sound level for the 500 Hz octave band.
%
%   If an overall level is to be calculated (i.e., across a set of octave bands), then this analysis
%   must be done for all octave band center frequencies.  Once the unweighted sound pressure
%   levels at the location are determned, the respective octave band A-weighting offsets are
%   applied.
%
%   The overall A-weighted sound pressure levels is then calculated logarithmically add using the
%   expression,
%
%   10*log10( sum( 10^(Lp_a/10) )


