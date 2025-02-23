


%% Synopsis

% Problem 1 - Modal Behaviour of a Cylindrical Room



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( './00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.5 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Define Room

room.radius = 3;  % m
room.length = 10;  % m



%% Test Circular Mode Function

% psi = circular_mode_shape( 3, 1, 2, false );  % 3.7261 - CHECKED FROM CLASS (PLOT NOT CREATED)
% psi = circular_mode_shape( 3, 1, 2, true );  % 3.7261 - CHECKED FROM CLASS (PLOT CREATED)



%% Define Anonymous Function for the Natural Frequencies

h_natural_frequencies = @( c, nx, ntheta, nr, Lx, cylinder_radius, plot_flag )  (c/2) .* sqrt( (nx/Lx).^2  + (circular_mode_shape(nr, ntheta, cylinder_radius, plot_flag )/cylinder_radius).^2 );



%% Calculate the Natural Frequencies

% The maximum number of radial modes is 5 (indexed from 0 to 4).
% The maximum number of angular modes is 8 (indexed from 0 to 7).

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

NUMBER_OF_LOWEST_FREQUENCIES = 11;
    mode_indices =  ( 1:1:NUMBER_OF_LOWEST_FREQUENCIES ).';

[ sortedValues, sortedIndices ] = sort( natural_frequencies(:) );  % 21-by-8-by-5 -> 840 elements

smallestValues = sortedValues( 1:NUMBER_OF_LOWEST_FREQUENCIES );
    % [ mode_indices   round( smallestValues, 1 ) ]
    %
    % 1            0
    % 2         17.2
    % 3         33.5
    % 4         34.3
    % 5         37.6
    % 6         47.9
    % 7         51.5
    % 8         55.6
    % 9         58.2
    % 10         61.4
    % 11         65.3

smallestIndices = sortedIndices( 1:NUMBER_OF_LOWEST_FREQUENCIES );

[ x, y, z ] = ind2sub( size(natural_frequencies), smallestIndices );
    % ( [ x y z ] - 1 )

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

% [ (1:11).'  abs( smallestValues - 53 ) ]

temp = [ x y z ] - 1;  temp( 7:8, :, : )

h_natural_frequencies( 343, 3, 0, 0, 10, 3, false );  % 51.5 Hz, (3, 0, 0)
h_natural_frequencies( 343, 0, 2, 0, 10, 3, false );  % 55.6 Hz, (0, 2, 0)



%% Part c

% See the report.



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

% https://www.mathworks.com/matlabcentral/answers/1883747-how-to-find-the-5-minimum-values-in-a-multidimensional-matrix-and-the-indices-to-which-these-entries


