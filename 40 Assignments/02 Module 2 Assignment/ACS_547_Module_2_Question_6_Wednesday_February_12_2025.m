


%% Synopsis

% Problem 6 - Close-fitting Enclosure Design



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



%% Define Compressor

compressor.width = 1;  % m
compressor.depth = 1;  % m
compressor.height = 2;  % m
    compressor.area = 2*(compressor.width * compressor.depth) + 2*(compressor.width * compressor.height) + 2*(compressor.depth * compressor.height);  % 3 m^2
    compressor.volume = compressor.width * compressor.depth * compressor.height;  % m^3

compressor.power_level = 105;  % dB re: 1e-12 Watts
compressor.frequency = 50;  % Hz

c = 343;  % m/s

rho0 = 1.21;  % kg/m^3



%% Sound Level Target

sound_level_target = 82;  % dB re: 20e-6 Pascals



%% Define Workshop

R = 40;  % m^2 or Sabins



%% Define Anonymous Functions

h_RA_term_1 = @( rho0, c , S, k, delta_mu, D, w )  ( rho0*c/S )  *  ( (k * sqrt( (2*3.178e-5) / (rho0*w) ) * D * 0.004 ) / (2*S) *1.4364 );
h_RA_term_2 = @( rho0, c , S, k, delta_mu, D, w, h )  ( rho0*c/S )  *  0.288*k*3.178e-5*log10((4*S)/(pi*h^2));
h_RA_term_3 = @( rho0, c , S, k, delta_mu, D, w, h )  ( rho0*c/S )  *  (0.5*S*k^2)/(2*pi);



%% Define Close-fitting Enclosure

helmholtz_factor = (2 * pi * compressor.frequency) / c;  % 0.92 m^-1
%
%  For a small enclosure, k*d << 1.  Therefore d << 1.1.

enclosure.width = compressor.width + 0.25;  % m
enclosure.depth = compressor.depth + 0.5;  % m
enclosure.height = 3;  % m;  compression height is 2 m
    enclosure.area = 2*(enclosure.width * enclosure.depth) + 2*(enclosure.width * enclosure.height) + 2*(enclosure.depth * enclosure.height);
    enclosure.volume = enclosure.width * enclosure.depth * enclosure.height;

enclosure.E = 3.6e9;  % Pascals
enclosure.thickness = 3.81e-2;  % m
enclosure.density = 800;  % kg/m^3
enclosure.poisson_ratio = 0.25;  % Unitless

% Clamped boundary conditions.



%% Calculate Diffuse Sound Pressure Level

% Assume distance is beyond the critical distance, so the distance value is
% large and its associated term is not relevant.

sound_pressure_level = 105  +  10*log10( 4/R );  % 95 dB SPL



%% Calculate the Required Insertion Loss

target_insertion_loss = sound_pressure_level - 82  % 13 dB



%% Insertion Loss

bending_stiffness = ( enclosure.E * enclosure.thickness^3 ) / ( 12*( 1 - enclosure.poisson_ratio^2 ) );
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
side_1.correction_factor = 0.6;
    side_1.compliance = h_wall_compliance( side_1.area, side_1.correction_factor );
%
% Side 2
side_2.compliance = side_1.compliance;


% Side 3
side_3.area = enclosure.width * enclosure.height;
side_3.aspect_ratio = max( enclosure.width, enclosure.height ) / min( enclosure.width, enclosure.height );
side_3.correction_factor = 0.5; 
    side_3.compliance = h_wall_compliance( side_3.area, side_3.correction_factor );
%
% Side 4
side_4.compliance = side_3.compliance;


estimated_insertion_loss = 20*log10( 1 + Ca / ( top.compliance + 2*side_1.compliance + 2* side_3.compliance ) );
    13 - estimated_insertion_loss



%% Compliance of the Air Intake

air_intake_radius = 10e-2;  % m
air_intake_thickness = enclosure.thickness;  % m
air_intake_frequency = 50;  % Hz
    air_intake_angular_frequency = 2*pi*air_intake_frequency;  % radians/s

viscosity = 1.5e-5;  % m^2/s

h = 0.3;

f = 50; 
    term_1 = h_RA_term_1( rho0, c, pi*(air_intake_radius*2)^2/4, 2*pi*f/c, sqrt( (2 * 3.178e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.1, 2*pi*f );
    term_2 = h_RA_term_2( rho0, c, pi*(air_intake_radius*2)^2/4, 2*pi*f/c, sqrt( (2 * 3.178e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.1, 2*pi*f, 0.3 );
    term_3 = h_RA_term_3( rho0, c, pi*(0.1)^2/4, 2*pi*f/c, sqrt( (2 * 3.178e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.1, 2*pi*f, 0.3 );
        impedance.real = term_1 + term_2 + term_3;


% Deng (1998)
epsilon = 1;
    L_o = air_intake_radius * ( 1.27 / (1 + 1.92*epsilon) - 0.086 );

L_e =  enclosure.thickness + 2*L_o;
    impedance.imaginary = 1j * rho0 * (2 * pi * f) * L_e / ( pi*0.1^2/4 );


impedance.net = impedance.real + impedance.imaginary;
    compliance_of_hole = 1 / impedance.net;
        Cl = abs( compliance_of_hole );

estimated_insertion_loss_with_hole = 20*log10( (Cl + Ca)  / ( Cl + ( top.compliance + 2*side_1.compliance + 2* side_3.compliance ) ) );
    13 - estimated_insertion_loss_with_hole


critical_frequency = c^2/(2*pi)*sqrt( enclosure.density * enclosure.thickness / bending_stiffness );  % 777.1 Hz



%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


