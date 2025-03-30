

% To do:

% remove comment statements

% verify indices for required values



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

% DVA introduces two resonances for the original resonance.  These new
% resonances will produce greater motion at the new resonance frequencies.



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

set( groot, 'DefaultFigurePosition', [ 100  750    750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Define Anonymous Functions

h_rpm_to_f = @( rpm )  ( rpm * 0.10471 ) / ( 2 * pi );
h_rpm_to_w = @( rpm )  ( rpm * 0.10471 );

h_f_to_rpm = @( f )  ( f*2*pi ) / 0.10471;

h_w_to_f = @( w )  w/(2*pi);
h_w_to_rpm = @( w )  h_f_to_rpm( w*2*pi);


h_force = @( force_mass, rotation_speed, load_distance  )  force_mass .* rotation_speed.^2 .* load_distance;



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
k = ks;  % 9,810 N\m

wo = sqrt( k / m);  % 9.4 radians\s
    fo = wo/(2*pi);  % 1.5 Hz


m_dva = 15;  % kg
k_dva = 1e5;  % N\m

wo_dva = sqrt( k_dva / m_dva );  % 81.7 radians\s
    fo_dva = wo_dva/(2*pi);  % 13.0 Hz



%% Blocked Frequencies (Equation 9.1 of the Notes)

m1 = m_dva;
k1 = 0;  % Spring above DVA mass;  does not exist;  set to zero.

m2 = m;
k2 = k;

k3 = k_dva;


w1 = sqrt( ( k1 + k3 ) / m1 );  % 81.7 radians\s;  smaller than the system resonant frequency.
    f1 = w1/(2*pi);  % 13 Hz

w2 = sqrt( ( k2 + k3 ) / m2 );  % 31.6 radians\s;  greater than the system resonant frequency.
    f2 = w2/(2*pi);  % 5.0 Hz



%% Coupled Frequencies (Equation 9.15 of the Notes)

mu_4 = ( k3^2 / (m1*m2) );  % 6.1 (kg\s)^4
    mu = ( mu_4 )^0.25;  % 49.1 kg\s

w(1) = 0.5*( ( w1^2 + w2^2 )  +  sqrt( ( w1^2 - w2^2 )^2 + 4*mu_4 ) );
w(2) = 0.5*( ( w1^2 + w2^2 )  -  sqrt( ( w1^2 - w2^2 )^2 + 4*mu_4 ) );
    w = sqrt( w );


w_minus = w(2);  % 8.9 radians\s  (lower than 31.6 and 81.7 radians\s)
    f_minus = w(2)/(2*pi);  % 1.4 Hz
%
w_plus = w(1);  % 87.1 radians\s (greater than 31.6 and 81.7 radians\s)
    f_plus = w(1)/(2*pi);  % 13.9 Hz


% The blocked frequencies are always between the lower and higher coupled frequencies.



%% Mode Shapes (Equations 9.18 and 9.19 of the Notes)

phi_plus = [ ...
    1; ...
    -(1/mu^2)*sqrt(m1/m2)*(w_plus^2 - w1^2), ...
    ];
%
%   1
%   -0.14


phi_minus = [ ...
    1; ...
    (1/mu^2)*sqrt(m1/m2)*(w1^2 - w_minus^2), ...
    ];
%
%   1
%   0.99



%% Admittance Using Equations from DVA Lecture

MAXIMUM_RPM = 3e3;

F0 = 1;

w = 0:0.01:h_rpm_to_w( MAXIMUM_RPM );
    f = w/(2*pi);
        rpm = h_f_to_rpm( f );


m = [ 15 110 ];
k = [  0  k_dva.*(1+0.1*1i)  k2.*(1+0.2*1i)  ];


admittance_modified = F0 ./ ( m(1)*m(2) )*k(2) ./ (( w.^2 - w_plus^2 ) .* ( w.^2 - w_minus^2) );
admittance_original = F0 ./ m(2) * 1 ./ ( w2^2 - w.^2 );


figure( 'Name', 'Admittance - Equations from DVA Lecture' ); ...
    h1 = loglog( rpm, abs( admittance_original ) );  hold on;
    h2 = loglog( rpm, abs( admittance_modified ) );
    h3 = line( [ 300 300 ], [ 1e-6 1e-2 ], 'Color', 'k', 'LineStyle', '--' );
        text( 220, 2e-2, '300 RPM' );
    h4 = line( [ h_f_to_rpm( f_minus )  h_f_to_rpm( f_minus ) ], [ 1e-7 1e3 ], 'Color', [ 0.72, 0.27, 1.00 ], 'LineStyle', '-.' );
    
    h5 = line( [ h_f_to_rpm( f1 )  h_f_to_rpm( f1 ) ], [ 1e-7 1e3 ], 'Color', [ 0.47, 0.67, 0.19 ], 'LineStyle', '-.' );
    h6 = line( [ h_f_to_rpm( f2 )  h_f_to_rpm( f2 ) ], [ 1e-7 1e3 ], 'Color', [ 0.47, 0.67, 0.19 ], 'LineStyle', '--' );
    
    h7 = line( [ h_f_to_rpm( f_plus )  h_f_to_rpm( f_plus ) ], [ 1e-7 1e3 ], 'Color', [ 0.72, 0.27, 1.00 ], 'LineStyle', '--' );    grid on;
    
    axis( [ 5 3e3  1e-7 1 ] );
        legend( [ h1 h2 h3 h4 h5 h6 h7 ], 'Original System', 'DVA System', 'Operating RPM', ...
            '$w_-$', '$w_1$', ...
            '$w_2$', '$w_+$', ...
            'Location', 'NorthWest', 'Interpreter', 'Latex' );
    xlabel( 'Rotation [RPM]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );

return

%% Admittance Using nDOF Function

% Damping can be added by using a dampings vector or appending a complex
% value to the respective spring stiffness.  THESE ARE NOT INTERCHANGABLE.


rpm = 0:0.1:MAXIMUM_RPM;
    f = h_rpm_to_f( rpm );
        w = f ./ (2*pi);

FRF = nDOF_direct_solution( m, k, [ 0 0 0 ], f, 'admittance' );  % Symmetric matrix.
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


figure( 'Name', 'Admittance - nDOF Calculated' ); ...
    h1 = loglog( h_f_to_rpm( f ), abs( FRF( :, 1, 1 ) ), 'Color', 'r' );  hold on;
    h2 = line( [ 300 300 ], [ 1e-6 1e-2 ], 'Color', 'k', 'LineStyle', '--' );
        text( 220, 2e-2, '300 RPM' );
    h3 = line( [ h_f_to_rpm( f_minus )  h_f_to_rpm( f_minus ) ], [ 1e-7 1e3 ], 'Color', [ 0.72, 0.27, 1.00 ], 'LineStyle', '-.' );
    
    h4 = line( [ h_f_to_rpm( f1 )  h_f_to_rpm( f1 ) ], [ 1e-7 1e3 ], 'Color', [ 0.47, 0.67, 0.19 ], 'LineStyle', '-.' );
    h5 = line( [ h_f_to_rpm( f2 )  h_f_to_rpm( f2 ) ], [ 1e-7 1e3 ], 'Color', [ 0.47, 0.67, 0.19 ], 'LineStyle', '--' );
    
    h6 = line( [ h_f_to_rpm( f_plus )  h_f_to_rpm( f_plus ) ], [ 1e-7 1e3 ], 'Color', [ 0.72, 0.27, 1.00 ], 'LineStyle', '--' );    grid on;
    
    axis( [ 2e1 MAXIMUM_RPM  1e-7 1 ] );
        legend( [ h1 h2 h3 h4 h5 h6 ], 'DVA System', 'Operating RPM', ...
            '$w_-$', '$w_1$', ...
            '$w_2$', '$w_+$', ...
            'Location', 'NorthWest', 'Interpreter', 'Latex' );
    xlabel( 'Rotation [RPM]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );


% temp2 = abs( FRF( :, 1, 1 ) );
    % round( temp2( 3001 ), 3, 'significant' )



%% Displacement

dynamic_force.mass = 1;  % kg
dynamic_force.radius = 0.5;  % m


figure( 'Name', 'Displacement' ); ...
    h1 = loglog( rpm, abs( FRF( :, 1, 1 ) ).*h_force( dynamic_force.mass, w, dynamic_force.radius ).', 'Color', 'r' );  hold on;
    h2 = line( [ 300 300 ], [ 1e-7 1e-3 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
        text( 220, 2e-3, '300 RPM' );
    h4 = line( [ h_f_to_rpm( f_minus )  h_f_to_rpm( f_minus ) ], [ 1e-8 1e3 ], 'Color', [ 0.72, 0.27, 1.00 ], 'LineStyle', '-.' );
    
    h5 = line( [ h_f_to_rpm( f1 )  h_f_to_rpm( f1 ) ], [ 1e-8 1e3 ], 'Color', [ 0.47, 0.67, 0.19 ], 'LineStyle', '-.' );
    h6 = line( [ h_f_to_rpm( f2 )  h_f_to_rpm( f2 ) ], [ 1e-8 1e3 ], 'Color', [ 0.47, 0.67, 0.19 ], 'LineStyle', '--' );
    
    h7 = line( [ h_f_to_rpm( f_plus )  h_f_to_rpm( f_plus ) ], [ 1e-8 1e3 ], 'Color', [ 0.72, 0.27, 1.00 ], 'LineStyle', '--' );    grid on;
        
        legend( [ h1 h2 h3 h4 h5 h6 ], 'DVA System', 'Operating RPM', ...
            '$w_-$', '$w_1$', ...
            '$w_2$', '$w_+$', ...
            'Location', 'NorthWest', 'Interpreter', 'Latex' );

    xlabel( 'Rotation [RPM]' );  ylabel( 'Displacment [m]' );
    axis( [ 2e1 MAXIMUM_RPM  1e-8 1 ] );


% temp2 = abs( FRF( :, 1, 1 ) ).*h_force( dynamic_force.mass, w, dynamic_force.radius ).';
    % round( temp2( 3001 ), 3, 'significant' )



%% Ground Force

dynamic_force.mass = 1;  % kg
dynamic_force.radius = 0.5;  % m


figure( 'Name', 'Ground Force' ); ...
    h1 = loglog( rpm, abs( FRF( :, 1, 1 ) ).*h_force( dynamic_force.mass, w, dynamic_force.radius ).'.*k2, 'Color', 'k' );  hold on;
    h2 = line( [ 300 300 ], [ 1e-4 1e-0 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
        text( 220, 2e-0, '300 RPM' );
    h3 = line( [ 6 MAXIMUM_RPM ], [ 100 100 ], 'Color', 'r' );
        legend( [ h1 h2 h3 ], 'Ground Force', 'Operating RPM', 'Maximum Floor Load', 'Location', 'SouthEast' );
    xlabel( 'Rotation [RPM]' );  ylabel( ' Force [N]' );
    axis( [ 2e1 MAXIMUM_RPM  1e-5 1e3 ] );


temp2 = abs( FRF( :, 1, 1 ) ).*h_force( dynamic_force.mass, w, dynamic_force.radius ).'*k2;
    % round( temp2( 3001 ), 3, 'significant' )



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


