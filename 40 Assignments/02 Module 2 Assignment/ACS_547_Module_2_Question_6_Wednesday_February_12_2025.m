


%% Synopsis

% Lecture 11, Wednesday, February 19, 2025

% The compressor elevated above the ground.



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



%% Define Compressor

compressor.width = 1;  % m
compressor.depth = 1;  % m
compressor.height = 2;  % m
    compressor.area = 2*(compressor.width * compressor.depth) + 2*(compressor.width * compressor.height) + 2*(compressor.depth * compressor.height);  % 3 m^2
    compressor.volume = compressor.width * compressor.depth * compressor.height;  % m^3

compressor.power_level = 105;  % dB re: 1e-12 Watts
compressor.frequency = 50;  % Hz

c = 343;  % m/s

rho0 = 1.2;  % kg/m^3  CHECK



%% Sound Level Target

sound_level_target = 82;  % dB re: 20e-6 Pascals



%% Define Workshop

R = 40;  % m^2 or Sabins



%% Define Close-fitting Enclosure

helmholtz_factor = (2 * pi * compressor.frequency) / c;  % 0.92 m
%
%  For a small enclosure, k*d << 1.  Therefore d << 1.1.
d = 0.75;
    % helmholtz_factor * d;  % 0.69

enclosure.width = compressor.width + d;  % m
enclosure.depth = compressor.depth + d;  % m
enclosure.height = 3;  % m
    enclosure.area = 2*(enclosure.width * enclosure.depth) + 2*(enclosure.width * enclosure.height) + 2*(enclosure.depth * enclosure.height);
    enclosure.volume = enclosure.width * enclosure.depth * enclosure.height;

enclosure.E = 3.6e12;  % Pascals
enclosure.thickness = 3.81e-2;  % m
enclosure.density = 800;  % kg/m^3
enclosure.poisson_ratio = 0.25;  % Unitless

% Clamped boundary conditions.



%% Define Air Intage

air_intake_radius = 10e-2;  % m



%% Calculate Diffuse Sound Pressure Level

% Assume distance is beyond the critical distance, so the distance value is
% large and its associated term is not relevant.

sound_pressure_level = 105  +  10*log10( 4/R );  % 95 dB SPL



%% Calculate the Required Insertion Loss

target_insertion_loss = sound_pressure_level - 82;  % 13 dB



%% Insertion Loss

% For the insertion loss to be high, we need:
%
%   1.)  Compliance of the air to be high;  volume of enclosure must be large.
%   2.)  Compliance of each enclosure wall to be low;  low area, high stiffness, edges clamped).
%           AREA IS THE DOMINATE FACTOR OVER VOLUME.


% The correction factor for clamped walls.  See Figure 12.4 on slide 9 of the Lecture 11 notes.
aspect_ratio = enclosure.height / enclosure.width;  % 1.7
    correction_factor = 2;  % Approximate value read from the Figure 12.4.

bending_stiffness = ( enclosure.E * enclosure.thickness^3 ) / ( 12*( 1 - enclosure.poisson_ratio^2 ) );  % 1.78e7
    h_wall_compliance = @( wall_area, correction_factor )  ( 0.001 * wall_area^3 * correction_factor ) / bending_stiffness;

Ca = enclosure.volume / ( rho0 * c^2 );




% Top
top.area = enclosure.width * enclosure.depth;
top.aspect_ratio = max( enclosure.width, enclosure.depth ) / min( enclosure.width, enclosure.depth );
top.correction_factor = 3.8;
    top.compliance = h_wall_compliance( top.area, top.correction_factor );


% Side 1
side_1.area = enclosure.depth * enclosure.height;
side_1.aspect_ratio = max( enclosure.width, enclosure.height ) / min( enclosure.width, enclosure.height );
side_1.correction_factor = 2;
    side_1.compliance = h_wall_compliance( side_1.area, side_1.correction_factor );

% Side 2
side_2.compliance = side_1.compliance;


% Side 3
side_3.area = enclosure.width * enclosure.height;
side_3.aspect_ratio = max( enclosure.width, enclosure.height ) / min( enclosure.width, enclosure.height );
side_3.correction_factor = 2; 
    side_3.compliance = h_wall_compliance( side_3.area, side_3.correction_factor );

% Side 4
side_4.compliance = side_3.compliance;


estimated_insertion_loss = 20*log10( 1 + Ca / ( top.compliance + 2*side_1.compliance + 2* side_3.compliance ) );  % 59.2 dB



%% Compliance of the Air Intake





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


