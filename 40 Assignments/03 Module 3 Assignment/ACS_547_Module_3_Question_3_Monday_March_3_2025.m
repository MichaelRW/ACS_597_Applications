


%% Synopsis

% Problem 3 - Washing Machine Dynamic Vibration Absorber (DVA) Design



%% Note(s)

% At 300 RPM, the load frequency, this system will isolated vibration more
% than the 2 DOF and 1 DOF systems.

% The DVA behaves like a Helmholtz resonator, working as a mechnical notch
% filter at a particular frequency.

% DVAs work well in systems that have one dominate mode.  A good
% example of this is building sway and its compensation.  Based on the
% design and construction of a building, it will typically have a single,
% dominate mode.

% DVA is useful if you can not get away from the resonance frequency;  can
% not tune a 2 DOF system when the resonance frequency is away from the
% forcing frequency.



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



%% Define Anonymous Functions

h_f_to_rpm = @( f )  ( f*2*pi ) / 0.10471;
h_w_to_rpm = @( w )  h_f_to_rpm( w*2*pi);
h_w_to_f = @( w )  w/(2*pi);
h_rpm_to_f = @( rpm )  ( rpm * 0.10471 ) / ( 2 * pi );



%% Washing Machine Information

machine.mass = 100;  % kg
machine.rotation_speed = 31.4;  % radians\s or 300 rotations per minute (RPM)
machine.load_mass = 10;  % kg
machine.static_displacement = 1e-2;  % m



%% Load Force Information

dynamic_force.mass = 1;  % kg
dynamic_force.radius = 0.5;  % m

% Maximum force applied to foundation is 100 N.



%% Dynamic Vibration Absorber (DVA)

maximum_allowable_displacement = 1e-2;  % m
    ks = ( machine.load_mass * 9.81) / maximum_allowable_displacement;  % 9,810 N\m


m = 110;  % kg;  washing machine mass and load mass
k = ks;  % N\m

wo = sqrt( k / m);  % 9.4 radians\s
    fo = wo/(2*pi);  % 1.5 Hz


% Add a Dynamic Vibration Absorber (DVA)
m_dva = 15;  % kg
% k_dva = 1800;  % N\m ->  7.7095e-05;  2.4411e-05;  0.23947 ->  f2 IS LESS THAN f1
k_dva = 1000;  % N\m ->  7.2494e-05;  2.2954e-05;  0.22518
% k_dva = 10;  % N\m ->  6.7593e-05;  2.1402e-05;  0.20995
% k_dva = 1;  % N\m ->  6.7552e-05;  2.1389e-05;  0.20983

wo_dva = sqrt( k_dva / m_dva );  % 8.2 radians\s
    fo_dva = wo_dva/(2*pi);  % 1.3 Hz

% return

%% Blocked Frequencies (Euqation 9.1 of the Notes)

m1 = m_dva;
k1 = 0;  % Spring above DVA mass;  does not exist;  set to zero.

m2 = m;
k2 = k;  %  N/m

k3 = k_dva;  % 1,800 N/m


w1 = sqrt( ( k1 + k3 ) / m1 );  % 10.95 radians\s;  smaller than the system resonant frequency.
    f1 = w1/(2*pi)  % 1.3 Hz


w2 = sqrt( ( k2 + k3 ) / m2 );  % 10.27 radians\s;  greater than the system resonant frequency.
    f2 = w2/(2*pi)  % 1.6 Hz

% return

%% Coupled Frequencies (Equation 9.15 of the Notes)

mu_4 = ( k3^2 / (m1*m2) );  % 1,963.6 (radians\s)^4
    mu = ( mu_4 )^0.25;  % 6.66 radians\s

w(1) = 0.5*( ( w1^2 + w2^2 )  +  sqrt( ( w1^2 - w2^2 )^2 + 4*mu_4 ) );
w(2) = 0.5*( ( w1^2 + w2^2 )  -  sqrt( ( w1^2 - w2^2 )^2 + 4*mu_4 ) );
    w = sqrt( w );


w_minus = w(2);              % 8.2 radians\s  (lower than 42.4 and 45.1 radians\s)
    f_minus = w(2)/(2*pi);  % 1.31 Hz
%
w_plus = w(1);                 % 12.6 radians\s (higher than 42.4 and 45.1 radians\s)
    f_plus = w(1)/(2*pi);    % 2.0 Hz
%
% Note(s):
%
%   The blocked frequencies are always between these two frequencies.

% return

%% Mode Shapes (Equations 9.18 and 9.19 of the Notes)

phi_plus = [ ...
    1; ...
    -(1/mu^2)*sqrt(m1/m2)*(w_plus^2 - w1^2), ...
    ];

phi_minus = [ ...
    1; ...
    (1/mu^2)*sqrt(m1/m2)*(w1^2 - w_minus^2), ...
    ];



%% Plot Admittance - Magnitude of Displacement over Force

F0 = 1;

w = 0:0.01:40*2*pi;
    f = w/(2*pi);


if ( 1 )

    m = [ 15 110 ];
    k = [  0  k_dva  k ];
    % k = [  0  100e3*(1+0.1*1j)  1.8e3*(1+0.1*1j)  ];
    dampings = [ 0 0 0 ];

    X2 = F0/(m(1)*m(2))*k(3)./((w.^2-w_plus^2).*(w.^2-w_minus^2));
    X2_original = F0./m(2)*1./(w2^2-w.^2);

else

    m = [ 110  15 ];
    k = [ k  k_dva  0 ];
    % dampings = [ 0 0 0 ];

    X2 = F0/(m(2)*m(1))*k(1)./((w.^2-w_plus^2).*(w.^2-w_minus^2));
    X2_original = F0./m(1)*1./(w2^2-w.^2);

end


figure( 'Name', 'Displacement with DVA' ); ...
    loglog( f, abs(X2_original) );  hold on;
    loglog( f, abs(X2) );  grid on;
        legend( 'Original Displacment', 'Displacement With DVA', 'Location', 'NorthWest' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Displacment [m]' );



%% Admittance

% Damping can be added by using a dampings vector or appending a complex
% value to the respective spring stiffness.  THESE ARE NOT INTERCHANGABLE.


% 350 RPM is 5.8 Hz


f = 0:0.01:5.8;  % No matrix singularity warning issued (set 0 to 0.1 if it occurs).

FRF = nDOF_direct_solution( m, k, dampings, f, 'admittance' );  % Symmetric matrix.
%
squeeze( FRF( numel( f ), :, :) );
%
% (1, 1) - Force on M1 and its associated displacement.
% (2, 2) - Force on M2 and its associated displacement.
%
% (1, 2) - Force on M1 and the displacement of M2.
% (2, 1) - Force on M2 and the displacment on M1.
%
%   These are interchangeable;  the same result.


figure( 'Name', 'Admittance of DVA' ); ...
    loglog( f, abs( FRF( :, 1, 1 ) ), 'Color', 'r' );  hold on;
    line( [ 5  5 ], [ 1e-5  1e-3 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
    text( 6, 3e-4, '5 Hz' );
        legend( 'Washing Machine', 'Load Frequency', 'Location', 'NorthWest' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );


temp2 = abs( FRF( :, 1, 1 ) );
    temp2( 501 )



%% Displacement

dynamic_force.mass = 1;  % kg
dynamic_force.radius = 0.5;  % m

h_force = @( force_mass, rotation_speed, load_distance  )  force_mass .* rotation_speed.^2 .* load_distance;

angular_velocity = f / (2*pi);
rpm = h_f_to_rpm( f );

temp = h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ).';


figure( 'Name', '' ); ...
    h1 = loglog( rpm, abs( FRF( :, 1, 1 ) ).*temp, 'Color', 'r' );  hold on;
    h2 = line( [ 300 300 ], [ 1e-6 1e-1 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
        text( 220, 1e-1, '5 Hz' );
        legend( [ h1  h2 ], 'Washing Machine', 'Load Frequency', 'Location', 'South' );
    xlabel( 'Rotation [RPM]' );  ylabel( 'Displacment [m]' );
    axis( [ 6 400  1e-8 1e1 ] );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_2d', '-dpdf', '-r0' );
    %             end


temp2 = abs( FRF( :, 1, 1 ) ).*temp;
    temp2( 501 )



%% Force Applied to the Ground

dynamic_force.mass = 1;  % kg
dynamic_force.radius = 0.5;  % m

h_force = @( force_mass, rotation_speed, load_distance  )  force_mass .* rotation_speed.^2 .* load_distance;

angular_velocity = f / (2*pi);
rpm = h_f_to_rpm( f );

temp = h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ).';


figure( 'Name', 'Force Applied to Ground by Raft' ); ...
    h1 = loglog( rpm, abs( FRF( :, 1, 1 ) ).*temp*k2, 'Color', 'b' );  hold on;
    h2 = line( [ 300 300 ], [ 1e-4 1e-1 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
    h3 = line( [ 6 400 ], [ 100 100 ], 'Color', 'r' );
        text( 220, 1e-1, '5 Hz' );
        legend( [ h1 h2 h3 ], 'Raft', 'Load Frequency', 'Maximum Floor Load', 'Location', 'South' );
    xlabel( 'Rotation [RPM]' );  ylabel( 'Force Applied to Ground [N]' );
    axis( [ 6 400  1e-5 1e3 ] );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_2e', '-dpdf', '-r0' );
    %             end


temp2 = abs( FRF( :, 1, 1 ) ).*temp*k2;
    temp2( 501 )



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


