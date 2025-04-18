


%% Synopsis

% Problem 4 - Traffic Noise Modeling

% See Lecture 22 on Monday, April 7, 2025
%   ""D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\22 Monday, April 7, 2025\Lecture 22 - Coherence effects - filled.pptx""



%% Note(s)

% Incoherent line source.

% Distance (source and receivers) are in units of meters.

% Source can be real-valued or complex-valued.  Complex-valued sources have magnitude and phase information.

% Sources have units of volume velocity, m^3/s.

% Complex pressures are in units of Pascals.



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

% set( groot, 'DefaultFigurePosition', [ 100  750    750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Define Anonymous Function(s)

h_distance_to_house = @( horizontal_distance, distance_to_house_from_line_source )  sqrt( horizontal_distance^2 + distance_to_house_from_line_source^2 );

h_source_pressure = @( Lw, distance )  Lw - 10*log10( distance );

h_sound_pressure_net = @( levels )  10*log10( sum( 10.^(levels./10) ) );



%% Contants and Parameters

source_separation_distance = 8;  % m

average_A_weighted_source_pressure = 86;  % dB SPL A at 1 m

distance_to_house_from_line_source = 200;  % m

D0 = 2;

Lw = 86 - 10*log10( D0 / (4*pi*1^2) );  % 94 dB



%% Part 4a

% SOURCE_OFFSET = 64;  % 75.74 dB

% SOURCE_OFFSET = 1e0;  % 75.74 dB
% SOURCE_OFFSET = 1e1;  % 84.07 dB
% SOURCE_OFFSET = 1e2;  % 91.18 dB
% SOURCE_OFFSET = 1e3;  % 94.38 dB
SOURCE_OFFSET = 1e4;  % 96.21 dB
% SOURCE_OFFSET = 1e5;  % 97.50 dB
% SOURCE_OFFSET = 1e6;  % 98.49 dB
% SOURCE_OFFSET = 1e7;  % 99.29 dB
% SOURCE_OFFSET = 1e8;  % 99.97 dB
% SOURCE_OFFSET = 1e9;  % 100.56 dB
    x_axis_multiplier = -SOURCE_OFFSET:1:SOURCE_OFFSET;


level_set = nan( numel( x_axis_multiplier ), 1 );

fprintf( 1, '\n' );

for index = 1:1:numel( x_axis_multiplier )

    horizontal_distance = abs( x_axis_multiplier(index)*source_separation_distance );

    working_distance = h_distance_to_house( horizontal_distance, distance_to_house_from_line_source );

    level_set( index ) = h_source_pressure( Lw, working_distance );

    % fprintf( 1, '%d\t%d m\t\tDistance:  %3.1f\t%3.1f\n', x_axis_multiplier(index), horizontal_distance, round( working_distance, 1 ), level_set(index) );

end

fprintf( 1, '\n' );

% h_sound_pressure_net( level_set )


number_of_sources = [ 1e2  1e3  1e4  1e5  1e6  1e7  1e8  1e9 ];
estimated_levels = [ 91.18  94.38  96.21  97.50  98.49  99.29  99.97  100.56 ];

figure( ); ...
    semilogx( number_of_sources, estimated_levels, 'Marker', '*' );  grid on;
    xlabel( 'Number of Sources' );  ylabel( 'Estimate Sound Pressure Level  [dB]' );
    xlim( [ 5e1 2e9 ] );



%% Part 4b

line_source_circumference = 2*pi*distance_to_house_from_line_source;

number_of_sources_around_the_house = line_source_circumference / 8;  % 157

single_source_pressure = h_source_pressure( Lw, 200 );  % 70.97 dB

sound_pressure_net_circle  = 10*log10( number_of_sources_around_the_house * 10.^(single_source_pressure./10) );  % 93 dB



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


