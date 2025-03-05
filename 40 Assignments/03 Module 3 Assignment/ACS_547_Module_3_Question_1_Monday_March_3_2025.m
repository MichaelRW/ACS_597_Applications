


%% Synopsis

% Problem 1 - Washing Machine Mount Design



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( './00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Placeholder

% logspace( log10(0.001), log10(1), 10 )

% linspace( 0.001, 1, 10 )

% [ FRF ] = nDOF_direct_solution( masses, stiffnesses, dampings, freq, FRF_type )



%% Define Machine

machine.mass = 100;  % kg
machine.rotation_speed = 31.4;  % radians\s
machine.load = 10;  % kg
machine.static_displacement = 1e-2;  % m



%% Define Force

dynamic_force.mass = 1;  % kg
dynamic_force.radius = 0.5;  % m


% Maximum force applied to foundation is 100 N.



%% Part 1a

h_force = @( force_mass, rotation_speed, load_distance  )  force_mass .* rotation_speed.^2 .* load_distance;


rpm = 0:1:350;  % rotations per minute
    rpm_conversion_to_radians_per_second = 0.10472;  % radians\s
        angular_velocity = rpm .* rpm_conversion_to_radians_per_second;


h_force( dynamic_force.mass, 300*0.10472, dynamic_force.radius );  % 494.4 N

force_frequency = 300 * 0.10472 / ( 2 * pi );  % Hz


% figure( ); ...
%     plot( rpm, h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ) );  hold on;
%     plot( 300, 494, 'Marker', '.', 'MarkerSize', 20 );
%     line( [ 300, 300 ], [ 400 600 ], 'Color', 'r' );
%     text( 305, 494, '494 N' );  grid on;
%     xlabel( 'Angular Speed [RPM]' );  ylabel( 'Dynamic Force Amplitude [N]' );
%     %
%     axis( [ -10 355  -5 705 ] );
    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Homework_3_Part_1', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex



%% Part 1b





%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


