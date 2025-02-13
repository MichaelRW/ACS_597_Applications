


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



%% Test Circular Mode Function

% psi = circular_mode_shape( 3, 1, 2, false );  % 3.7261 - CHECKED FROM CLASS (PLOT NOT CREATED)
% psi = circular_mode_shape( 3, 1, 2, true );  % 3.7261 - CHECKED FROM CLASS (PLOT CREATED)



%% Natural Frequencies Function

h_natural_frequencies = @( c, nx, ntheta, nr, Lx, cylinder_radius, plot_flag )  (c/2) .* sqrt( (nx/Lx).^2  + (circular_mode_shape(nr, ntheta, cylinder_radius, plot_flag )/cylinder_radius).^2 );


NX_SIZE = 20;
NTHETA_SIZE = 7;
NR_SIZE = 4;
    natural_frequencies = nan( NX_SIZE, NTHETA_SIZE, NR_SIZE );

for nx = 0:1:NX_SIZE
    for ntheta = 0:1:NTHETA_SIZE
        for nr = 0:1:NR_SIZE
            natural_frequencies( nx+1, ntheta+1, nr+1 ) = h_natural_frequencies( 343, nx, ntheta, nr, 10, 3, false );
        end
    end
end



%% Part a - Find 10 Lowest Resonance Frequencies

% https://www.mathworks.com/matlabcentral/answers/1883747-how-to-find-the-5-minimum-values-in-a-multidimensional-matrix-and-the-indices-to-which-these-entries

% nr MIGHT ALWAYS be zero.

NUMBER_OF_LOWEST_FREQUENCIES = 11;

[ sortedValues, sortedIndices ] = sort( natural_frequencies(:) );  % 21-by-8-by-5 -> 840 elements

smallestValues = sortedValues( 1:NUMBER_OF_LOWEST_FREQUENCIES );
%
% 0
% 17.15
% 33.505
% 34.3
% 37.64
% 47.949
% 51.45
% 55.577
% 58.163
% 61.398
% 65.31

smallestIndices = sortedIndices(1:NUMBER_OF_LOWEST_FREQUENCIES);


[ x, y, z ] = ind2sub( size(natural_frequencies), smallestIndices );
    [ x y z ] - 1;


% Verify the calculated mode indices.
h_natural_frequencies( 343, 0, 0, 0, 10, 3, false );  % 0 Hz
h_natural_frequencies( 343, 1, 0, 0, 10, 3, false );  % 17.2 Hz
h_natural_frequencies( 343, 0, 1, 0, 10, 3, false );  % 33.5 Hz
h_natural_frequencies( 343, 2, 0, 0 , 10, 3, false );  % 34.3 Hz
h_natural_frequencies( 343, 1, 1, 0 , 10, 3, false );  % 37.6 Hz
h_natural_frequencies( 343, 2, 1, 0 , 10, 3, false );  % 48.0 Hz
h_natural_frequencies( 343, 3, 0, 0 , 10, 3, false );  % 51.5 Hz
h_natural_frequencies( 343, 0, 2, 0 , 10, 3, false );  % 55.6 Hz
h_natural_frequencies( 343, 1, 2, 0 , 10, 3, false );  % 58.2 Hz
h_natural_frequencies( 343, 3, 1, 0 , 10, 3, false );  % 61.4 Hz
h_natural_frequencies( 343, 2, 2, 0 , 10, 3, false );  % 65.3 Hz



%% Part b - Two 

[ (1:11).'  abs( smallestValues - 53 ) ]

temp = [ x y z ] - 1;
    temp( 7:8, :, : )

% Modes:
%   ( 3, 0, 0 ) and ( 0, 2, 0 )

h_natural_frequencies( 343, 3, 0, 0, 10, 3, true );  % 51.5 Hz

h_natural_frequencies( 343, 0, 2, 0, 10, 3, true );  % 55.6 Hz



%% Part c

% For mode ( 3, 0, 0 ), place the source in the center of the cylinder at 5 meters.

% For mode ( 0, 2, 0 ), place the source in the center of the cylinder.



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


