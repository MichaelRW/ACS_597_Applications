


%% Synopsis

% Problem 1 - Washing Machine Mount Design - 1 Degree of Freedom (DOF)


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



%% To Do

% 1c - Transmission loss is in dB.  It is the 10*log10 of transmissibility.
%       Check plots and make the context clear in the write-up.

%



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format LongG;

pause( 1 );

PRINT_FIGURES = 0;



%% Washing Machine Information

machine.mass = 100;  % kg
machine.rotation_speed = 31.4;  % radians\s or 300 rotations per minute (RPM)
machine.load_mass = 10;  % kg
machine.static_displacement = 1e-2;  % m



%% Load Force Information

dynamic_force.mass = 1;  % kg
dynamic_force.radius = 0.5;  % m

% Maximum force applied to foundation is 100 N.



%% Part 1a

h_force = @( force_mass, rotation_speed, load_distance  )  force_mass .* rotation_speed.^2 .* load_distance;


rpm = 0:1:350;  % RPM
    rpm_conversion_to_radians_per_second = 0.10472;
        angular_velocity = rpm .* rpm_conversion_to_radians_per_second;  % radians\s


% Characteristics of the dynamic load.
force_frequency = 300 * rpm_conversion_to_radians_per_second / ( 2 * pi );  % 5 Hz
h_force( dynamic_force.mass, 300*rpm_conversion_to_radians_per_second, dynamic_force.radius );  % 493.5 N


figure( 'Name', 'Load Force Versus Angular Velocity' ); ...
    plot( rpm, h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ) );  hold on;
    plot( 300, 494, 'Marker', '.', 'MarkerSize', 20 );
    line( [ 300, 300 ], [ 400 600 ], 'Color', 'r' );
    text( 305, 494, '493.5 N' );  grid on;
    xlabel( 'Angular Speed [RPM]' );  ylabel( 'Load Force [N]' );
    axis( [ -10 355  -5 705 ] );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1a', '-dpdf', '-r0' );
    %             end

return

%% Part 1b

total_mass = machine.mass + machine.load_mass;  % 110 kg

maximum_allowable_displacement = 1e-2;  % m

ks = (total_mass * 9.81) / maximum_allowable_displacement;  % 1.0791e5 N\m

wo = sqrt( ks / total_mass );  % 31.3 radians\s
    fo = wo / (2*pi);  % 4.98 Hz - Natural frequency of the mount.



%% Part 1c

frequency_set = 0.1:0.1:1e3;

r = frequency_set ./ fo;

damping_ratio = logspace( log10(0.001), log10(1), 6 );
    epsilon = damping_ratio;

h_transmissibility = @( epsilon, r )  sqrt( ( 1 + (2.*epsilon.*r).^2 ) ./ ( ( 1 - r.^2 ).^2 + ( 2.*epsilon.*r ).^2 ) );


figure( 'Name', 'Problem 1c - Transmission Loss' ); ...
    hold on;
    %
    for index = 1:1:numel( epsilon )
        plot( r, h_transmissibility( epsilon( index ), r ) );
    end
    %
    xlabel( 'Frequency Ratio [$\frac{f}{fo}$;  unitless]' );  ylabel( 'Force Transmissibility [$T_F$;  unitless]' );
        legend( '0.001', '0.004', '0.015', '0.063', '0.25', '1', 'Location', 'NorthEast' );
    axis( [ 0.1 1e2  1e-4 5e2 ] );
    grid on;
    set( gca, 'XScale', 'log', 'YScale', 'log' );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1c', '-dpdf', '-r0' );
    %             end

% return

%% Part 1d

rpm = 0:1:1e3;  % rotations per minute
    rpm_conversion_to_radians_per_second = 0.10472;  % radians\s
        angular_velocity = rpm .* rpm_conversion_to_radians_per_second;  % radians\s

frequency_set = angular_velocity / (2*pi);
    r = frequency_set ./ fo;


figure( ); ...
    hold on;
    %
    for index = 1:1:numel( epsilon )
        plot( rpm, h_transmissibility( epsilon( index ), r ) .* h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ) );
    end
    %
    line( [ 300, 300 ], [ 200 1e6 ], 'Color', 'r' );
        text( 305, 300, '5 Hz' );  grid on;
    line( [ 0 350 ], [ 100 100 ], 'Color', 'r' );
        text( 0, 205, '100 N' );  grid on;    
    xlabel( 'Angular Speed [RPM]' );  ylabel( 'Force Applied to Foundation [N]' );
        legend( '0.001', '0.004', '0.015', '0.063', '0.25', '1', 'Location', 'SouthEast' );
    axis( [ -10 355  0 1e6 ] );
    grid on;
    set( gca, 'YScale', 'log' );
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1d', '-dpdf', '-r0' );
    %             end



%% Part 1e

% Damping ratio is epsilon.
%
% Damping coefficient is C -> epsilon * 2 * sqrt( ks * m )


% The best damping ratio is 1 kg^0.5 m\s

C = 1 * 2*sqrt(ks*total_mass);  % 6,890.6

% ONLY the system mass


% 10 kg static for the displacement.



% The units of a viscous damping coefficient are Newton-seconds per meter (Ns/m) or kilograms per second (kg/s).


% Plot admittance as function of frequency.
m = total_mass;
k = [ ks  0 ];
dampings = [ 1  0 ];

f = frequency_set;

admittance = nDOF_direct_solution( m, k, dampings, f, 'admittance' );

% H = FRF( :, 1, 1 );  % Need to multiply this be frequency dependent forcing function.

figure( 'Name', 'Admittance' ); ...
    semilogy( rpm, abs( admittance ) );  grid on;
    xlabel( 'Angular Speed [RPM]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
    set( gca, 'YScale', 'log' );
    axis( [ -10 355  0 0.09 ] );
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1e', '-dpdf', '-r0' );
    %             end


% The admittance at 5 Hz (300 RPM) is 0.00152 m\N.



%% Part 1f

displacement = abs( admittance ) .* h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ).';


figure( 'Name', 'Displacement' ); ...
    plot( rpm, displacement );  grid on;
    xlabel( 'Angular Speed [RPM]' );  ylabel( 'Displacement [m]' );
    set( gca, 'YScale', 'log' );
    axis( [ -10 355  0 1e6 ] );
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1f', '-dpdf', '-r0' );
    %             end


% The displacement is 0.75 meters at 300 RPM.



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




