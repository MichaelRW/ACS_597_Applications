


%% Synopsis

% Problem 1 - Washing Machine Mount Design - 1 Degree of Freedom (DOF)



%% To Do

% 1c - Transmission loss is in dB.  It is the 10*log10 of transmissibility.
%       Check plots and make the context clear in the write-up.



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


rpm = 0:1:350;
    rpm_conversion_to_radians_per_second = 0.10472;
        angular_velocity = rpm .* rpm_conversion_to_radians_per_second;  % radians\s


% Characteristics of the dynamic load.
force_frequency = 300 * rpm_conversion_to_radians_per_second / ( 2 * pi );  % 5 Hz
h_force( dynamic_force.mass, 300*rpm_conversion_to_radians_per_second, dynamic_force.radius );  % 493.5 N


figure( 'Name', 'Problem 1a - Load Force Versus Angular Velocity' ); ...
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



%% Part 1b

maximum_allowable_displacement = 1e-2;  % m
    ks = ( machine.load_mass * 9.81) / maximum_allowable_displacement;  % 9,810 N\m


total_mass = machine.mass + machine.load_mass;  % 110 kg
    wo = sqrt( ks / total_mass );  % 9.4 radians\s
        fo = wo / (2*pi);  % 1.5 Hz - Natural frequency of the mount.



%% Part 1c

frequency_set = 0.1:0.01:1e3;
    r = frequency_set ./ fo;

damping_ratio = logspace( log10(0.001), log10(1), 6 );
    
h_transmissibility = @( epsilon, r )  sqrt( ( 1 + (2.*epsilon.*r).^2 ) ./ ( ( 1 - r.^2 ).^2 + ( 2.*epsilon.*r ).^2 ) );


figure( 'Name', 'Problem 1c - Transmission Loss' ); ...
    hold on;
    %
    for index = 1:1:numel( damping_ratio )
        plot( r, 10*log10( h_transmissibility( damping_ratio( index ), r ) ) );
    end
    %
    xlabel( 'Frequency Ratio [$\frac{f}{fo}$;  unitless]' );  ylabel( 'Tranmission Loss [dB]' );
        legend( '0.001', '0.004', '0.015', '0.063', '0.25', '1', 'Location', 'NorthEast' );
    axis( [ 0.1 665  -60 25 ] );
    grid on;  box on;
    set( gca, 'XScale', 'log' );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1c', '-dpdf', '-r0' );
    %             end



%% Part 1d

rpm = 0:1:350;  % rotations per minute
    rpm_conversion_to_radians_per_second = 0.10472;  % radians\s
        angular_velocity = rpm .* rpm_conversion_to_radians_per_second;  % radians\s

frequency_set = angular_velocity / (2*pi);
    r = frequency_set ./ fo;


figure( 'Name', 'Problem 1d - Applied Force to Foundation' ); ...
    hold on;
    %
    for index = 1:1:numel( damping_ratio )
        plot( rpm, h_transmissibility( damping_ratio( index ), r ) .* h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ) );
    end
    %
    line( [ 300, 300 ], [ 2 1e4 ], 'Color', 'k', 'LineStyle', '--' );
        text( 305, 2, '5 Hz' );  grid on;
    line( [ 0 350 ], [ 100 100 ], 'Color', 'k', 'LineStyle', '--' );
        text( 0, 205, '100 N' );  grid on;    
    xlabel( 'Angular Speed [RPM]' );  ylabel( 'Force Applied to Foundation [N]' );
        legend( '0.001', '0.004', '0.015', '0.063', '0.25', '1', 'Location', 'South' );
    axis( [ -10 355  1 2e4 ] );
    grid on;  box on;
    set( gca, 'YScale', 'log' );
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1d', '-dpdf', '-r0' );
    %             end



%% Part 1e

% The best damping ratio is 0.25.

C = 0.25 * 2*sqrt( ks * total_mass );  % 519.4 kg\s

% The units of a viscous damping coefficient are Newton-seconds per meter (Ns/m) or kilograms per second (kg/s).


m = total_mass;
k = [ ks  0 ];
dampings = [ C  0 ];
    admittance = nDOF_direct_solution( m, k, dampings, frequency_set, 'admittance' );


figure( 'Name', 'Problem 1e - Admittance' ); ...
    semilogy( rpm, abs( admittance ) );  grid on;
    xlabel( 'Angular Speed [RPM]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
    set( gca, 'YScale', 'log' );
    %
    line( [ 300, 300 ], [ 7e-6 2e-5 ], 'Color', 'k', 'LineStyle', '--' );
        text( 305, 1.2e-5, '5 Hz' );  grid on;
    axis( [ -10 355  6e-6 2.5e-4 ] );
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1e', '-dpdf', '-r0' );
    %             end

% return

%% Part 1f

displacement = abs( admittance ) .* h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ).';


figure( 'Name', 'Displacement' ); ...
    plot( rpm, displacement );  grid on;
    xlabel( 'Angular Speed [RPM]' );  ylabel( 'Displacement [m]' );
    line( [ 300, 300 ], [ 1e-3 1e-2 ], 'Color', 'k', 'LineStyle', '--' );
        text( 305, 1e-3, '5 Hz' );  grid on;
    set( gca, 'YScale', 'log' );
    axis( [ -10 355  1e-7 2e-2 ] );
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_1f', '-dpdf', '-r0' );
    %             end



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




