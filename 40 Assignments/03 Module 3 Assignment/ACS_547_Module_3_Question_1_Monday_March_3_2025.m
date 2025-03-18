


%% Synopsis

% Problem 1 - Washing Machine Mount Design


% Part 1 - Single-stage Isolation
%
% The challenge here is that the mass does no go away.

% Part 2 - Two-stage Isolation
%
% A 2 degrees-of-freedom (DOF) system is better than a 1 DOF system.

% Part 3 - Dynamic Vibration Absorption (DVA) System
%
% At 300 RPM, this system will attenuate\reduce vibration amplitude more
% than the 2 DOF and 1 DOF systems.


% The DVA system acts like a Helmholtz resonator, working as a mechanical
% notch filter at a particular frequency.

% DVA systems work well in systems that have one dominate mode.  A good
% example of this is building sway and its compensation.  Based on the
% design and construction of a building, it will typically have a single,
% dominate mode.



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

PRINT_FIGURES = 0;



%% Define Machine

machine.mass = 100;  % kg
machine.rotation_speed = 31.4;  % radians\s
machine.load_mass = 10;  % kg
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


force_frequency = 300 * 0.10472 / ( 2 * pi );  % 5 Hz

h_force( dynamic_force.mass, 300*0.10472, dynamic_force.radius );  % 494.4 N (amplitude)


% figure( ); ...
%     plot( rpm, h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ) );  hold on;
%     plot( 300, 494, 'Marker', '.', 'MarkerSize', 20 );
%     line( [ 300, 300 ], [ 400 600 ], 'Color', 'r' );
%     text( 305, 494, '494.4 N' );  grid on;
%     xlabel( 'Angular Speed [RPM]' );  ylabel( 'Force Amplitude [N]' );
%     %
%     axis( [ -10 355  -5 705 ] );
% 
%     % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
%     % 
%     set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
%         pos = get( gcf, 'Position' );
%             set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
%                 if ( PRINT_FIGURES == 1 )
%                     print(gcf, 'Homework_3_Part_1a', '-dpdf', '-r0' );
%                 end
% %
% % https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex



%% Part 1b

total_mass = machine.mass + machine.load_mass;  % 110 kg

maximum_allowable_displacement = 1e-2;  % m

ks = (total_mass * 9.81) / maximum_allowable_displacement;  % 1.0791e5 N\m

wo = sqrt( ks / total_mass );  % 31.3 radians\s
    fo = wo / (2*pi);  % 4.98 Hz



%% Part 1c

frequency_set = 0.1:0.1:1e3;

r = frequency_set ./ fo;

damping_ratio = logspace( log10(0.001), log10(1), 6 );
    epsilon = damping_ratio ./ ( 2.*sqrt(ks.*total_mass) );

h_transmissibility = @( epsilon, r )  sqrt( ( 1 + (2.*epsilon.*r).^2 ) ./ ( ( 1 - r.^2 ).^2 + ( 2.*epsilon.*r ).^2 ) );


figure( ); ...
    plot( nan, nan );  hold on;
    %
    for index = 1:1:numel( epsilon )
        plot( r, h_transmissibility( epsilon( index ), r ) );
    end
    %
    grid on;
    set( gca, 'XScale', 'log', 'YScale', 'log' );
    xlabel( 'Frequency Ratio [$\frac{f}{fo}$;  unitless]' );  ylabel( 'Force Transmissibility [$T_F$;  unitless]' );
    legend( '0.001', '0.004', '0.015', '0.063', '0.25', '1', 'Location', 'NorthEast' );
    %
    axis( [ 0.1 1e2  1e-4 5e2 ] );

    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    % 
    set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
        pos = get( gcf, 'Position' );
            set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
                if ( PRINT_FIGURES == 1 )
                    print(gcf, 'Homework_3_Part_1c', '-dpdf', '-r0' );
                end
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex



return

%% Part 1d

return

%% Part 1e

return

%% Part 1f

return

%% Strikeforce

% [ FRF ] = nDOF_direct_solution( masses, stiffnesses, dampings, freq, FRF_type )

% The FRF output of this function is a 3-dimensional vector with,
%
%   First dimension:  The number of frequencies at which the calculation is made.
%   Second dimension:  The driving force on the masses.
%   Third dimension:  The response of the masses.

% This matrix will be symmetrical on the second and third dimensions.


%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)














% % Displacement Versus Frequency Ratio
% 
% Fo = 1;
% wf = 4;  % Force frequency, radians\s
% 
% lambda = 1;
%     k = (2*pi)/lambda;
% 
% phase_offset = 0;
% 
% w = 5;
%     r = wf / w;
% 
% time_indices = 0:1e-2:5;
% 
% h_x = @( Fo, k, wf, time_indices, phase_offset, r, epsilon ) ( Fo./k.*sin(4*time_indices - phase_offset) ) ./ ( sqrt( (1 - r^2)^2 + (2*r*epsilon)^2 ) );
% 
% figure( ); ...
%     epsilon = 0;
%         % plot( time_indices, h_x( Fo, k, wf, time_indices, phase_offset, r, epsilon ) );  hold on;
%     %
%     for wf = 10:-1:1
%         wf
%         epsilon = 0.25;
%         r = wf / w;
%             plot( h_x( Fo, k, wf, time_indices, phase_offset, r, epsilon ) );  hold on;
%         keyboard
%     end
%     %
%     legend( '', '' );
%     xlabel( 'Frequency Ratio [WU]' );  ylabel( 'Normalized Maximum Displacment [WU]' );




