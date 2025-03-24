


%% Synopsis

% Problem 2 - Washing Machine Two-stage Mount Design

% See Lecture 14, Monday, March 3, 2025



%% Note(s):

% Allowing the top and bottom masses, the washing machine and the raft,
% respectively, to each be attached to ground independently.


% With the force applied to the washing machine and the displacement of the
% washing machine, use FSF( :, 1, 1 ).

% However, if the force is applied elsewhere, than FSF( :, 2, 2, ) might
% need to be used.



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

% return

%% Problem 2a

% See report.



%% Problem 2b - Blocked Frequencies of Resonance

% The frequency of the forcing function is 5 Hz and the amplitude is 493.5 N.
fo = 5;  % Hz
    wo = 2*pi*fo;  % 31.4 radians\s


m1 = 100 + 10;  % kg - Total mass of the washing machine and the load.
    k1 = 1e2;  % N\m
    c1 = 5;

m2 = 100;  % kg - UPPER LIMIT
    k2 = 1e2;  % N\m
    c2 = 5;

% m2 = m1;  k2 = k1;
% m1 = 100;  m2 = 10;


k3 = 1e3;  % Admittance:  9.31e-6 m\N  KEEP
% k3 = 1e2;  % Admittance:  9.23e-6 m\N
%
% Note(s):
%
%   a.) Strong coupling produces a high w+.


w1 = sqrt( ( k1 + k3 ) / m1 );  %  radians\s;  washing machine
    f1 = w1 / (2*pi);
        rpm1 = f1*(2*pi)/0.10471;

w2 = sqrt( ( k2 + k3 ) / m2 );  %  radians\s;  raft
    f2 = w2 / (2*pi);
        rpm2 = f2*(2*pi)/0.10471;

% Note:  If k3 is set to zero, then these frequencies represent the
% resonance frequencies for each mass on their own.



%% Problem 2b - Coupled Frequencies

mu_4 = ( k3^2 / (m1*m2) );  %  (radians\s)^4
    mu = ( mu_4 )^0.25;

w(1) = 0.5*( ( w1^2 + w2^2 )  +  sqrt( ( w1^2 - w2^2 )^2 + 4*mu_4 ) );
w(2) = 0.5*( ( w1^2 + w2^2 )  -  sqrt( ( w1^2 - w2^2 )^2 + 4*mu_4 ) );
    w = sqrt( w );


w_minus = w(2);              %  0.97 radians\s  (lower than w1, 1, and w2, 1.05 radians\s)
    f_minus = w(2)/(2*pi);  %  0.1545 Hz
        rpm_minus = h_f_to_rpm( f_minus );
%
w_plus = w(1);                 %  1.08 radians\s  (higher than w1, 1, and w2, 1.05 radians\s)
    f_plus = w(1)/(2*pi);    %  0.17124 Hz
        rpm_plus = h_f_to_rpm( f_plus );
%
% Note(s):
%
%   The blocked frequencies, w1 and w2, are always between these two frequencies.


clear w mu_4;



%% Problem 2b - Mode Shapes

phi_plus = [ ...
    1; ...
    -(1/mu^2)*sqrt(m1/m2)*(w_plus^2 - w1^2), ...
    ];

phi_minus = [ ...
    1; ...
    (1/mu^2)*sqrt(m1/m2)*(w1^2 - w_minus^2), ...
    ];


% Note(s):
%
% When m1 = m2 and k1 = k3,
%
%   a.)  At w_minus, m1 and m2 move the same amount and are in-phase.
%   b.)  At w_plus, m1 and m2 move the same amount, but are out-of-phase.

% In general, for a 2 DOF system there will be 2 modes (in-phase and out-of-phase).


% With different initial conditions (no forcing), the behaviour is a weighted sum of the two modes.

% return

%% Problem 2c

% masses = [ m2 m1 ];
% stiffnesses = [ k2  k3  k1 ];
% dampings = [ c2  0  c1 ];

masses = [ m1 m2 ];
stiffnesses = [ k1  k3  k2 ];
dampings = [ c1  0  c2 ];

rpm = 0:0.1:500;  % rotations per minute
    rpm_conversion_to_radians_per_second = 0.10472;  % radians\s
        angular_velocity = rpm .* rpm_conversion_to_radians_per_second;
            frequencies = angular_velocity / (2*pi);

FRF = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );


admittance = zeros( numel( frequencies ), 1 );

for index = 1:1:numel( frequencies )
    temp = diag( squeeze( FRF( index, 1, 1 ) ) );  % Check indexing to ensure force on washing machine.
        admittance( index ) = abs( temp );
end
%
clear temp temp2;

washing_machine.admittance = admittance;


admittance = zeros( numel( frequencies ), 1 );

for index = 1:1:numel( frequencies )
    temp = diag( squeeze( FRF( index, 1, 2 ) ) );  % Check indexing to ensure force on washing machine.
        admittance( index ) = abs( temp );
end
%
clear temp temp2;

raft.admittance = admittance;


figure( 'Name', 'Admittance - Magnitude' ); ...
    h1 = loglog( rpm, washing_machine.admittance, 'Color', 'r' );  hold on;
    h2 = loglog( rpm, raft.admittance, 'Color', 'b' );
    h3 = line( [ 300 300 ], [ 10^-9 10^-2 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
    %
    h4 = line( [ rpm_minus rpm_minus], [ 1e-8 1e1 ], 'Color', 'g', 'LineStyle', '-.' );
    %
    line( [ rpm1 rpm1 ], [ 1e-8 1e1 ], 'Color', 'm', 'LineStyle', '--' );
    line( [ rpm2 rpm2 ], [ 1e-8 1e1 ], 'Color', 'm', 'LineStyle', '--' );
    %
    h5 = line( [ rpm_plus rpm_plus ], [ 1e-8 1e1 ], 'Color', 'g', 'LineStyle', '-' );
    %
    legend( [ h1 h2 h3 h4 h5 ], 'Washing Machine', 'Raft', 'Load Frequency', 'w-', 'w+', 'Location', 'NorthEast' );
    %
    xlabel( 'Rotation [RPM]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
    axis( [ 6 400  1e-8 1e1 ] );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_2c', '-dpdf', '-r0' );
    %             end


washing_machine.admittance( 3001 )
washing_machine.admittance( 3001 ) < 9.99049512983748e-06

( 9.99049512983748e-06 - washing_machine.admittance( 3001 ) ) / washing_machine.admittance( 3001 ) * 100



admittance = zeros( numel( frequencies ), 1 );

for index = 1:1:numel( frequencies )
    temp = diag( squeeze( FRF( index, 1, 1 ) ) );  % Check indexing to ensure force on washing machine.
        admittance( index ) = unwrap( angle ( temp ) ) * 180 / pi;
end
%
clear temp temp2;

washing_machine.admittance_phase = admittance;


admittance = zeros( numel( frequencies ), 1 );

for index = 1:1:numel( frequencies )
    temp = diag( squeeze( FRF( index, 1, 2 ) ) );  % Check indexing to ensure force on washing machine.
        admittance( index ) = unwrap( angle ( temp ) ) * 180 / pi;
end
%
clear temp temp2;

raft.admittance_phase = admittance;

figure( 'Name', 'Admittance - Phase' ); ...
    h1 = semilogx( rpm, washing_machine.admittance_phase, 'Color', 'r' );  hold on;
    h2 = semilogx( rpm, raft.admittance_phase   , 'Color', 'b', 'LineStyle', '--' );  grid on;
    h3 = line( [ 300 300 ], [ -200 200 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
    %
    h4 = line( [ rpm_minus rpm_minus], [ -200 200 ], 'Color', 'g', 'LineStyle', '-.' );
    %
    line( [ rpm1 rpm1 ], [ -200 200 ], 'Color', 'm', 'LineStyle', '--' );
    line( [ rpm2 rpm2 ], [ -200 200 ], 'Color', 'm', 'LineStyle', '--' );
    %
    h5 = line( [ rpm_plus rpm_plus ], [ -200 200 ], 'Color', 'g', 'LineStyle', '-' );
    %
    legend( [ h1 h2 h3 h4 h5 ], 'Washing Machine', 'Raft', 'Load Frequency', 'w-', 'w+', 'Location', 'NorthEast' );
    %
    xlabel( 'Rotation [RPM]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
    xlim( [ 6 400 ] );
    % axis( [ 6 400  1e-8 1e1 ] );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_2c', '-dpdf', '-r0' );
    %             end

% return

%% Problem 2d

% dynamic_force.mass = 1;  % kg
% dynamic_force.radius = 0.5;  % m
% 
% 
% % rpm = 0:1:350;  % rotations per minute
% %     rpm_conversion_to_radians_per_second = 0.10472;  % radians\s
% %         angular_velocity = rpm .* rpm_conversion_to_radians_per_second;
% %             frequencies = angular_velocity / (2*pi);
% 
% h_force = @( force_mass, rotation_speed, load_distance  )  force_mass .* rotation_speed.^2 .* load_distance;
% 
% 
% temp = h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ).';
% 
% 
% figure( 'Name', '' ); ...
%     h1 = loglog( rpm, washing_machine.admittance.*temp );  hold on;
%     h2 = loglog( rpm, raft.admittance.*temp );
%     h3 = line( [ 300 300 ], [ 10^-9 10^-2 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
%         legend( [ h1 h2 h3 ], 'Washing Machine', 'Raft', 'Load Frequency', 'Location', 'SouthWest' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Displacment [m]' );
%     % axis( [ 0 max(frequencies)  1e-9 1e-2 ] );
%     % 
%     % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
%     %     pos = get( gcf, 'Position' );
%     %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
%     %             if ( PRINT_FIGURES == 1 )
%     %                 print(gcf, 'Homework_3_Part_2c', '-dpdf', '-r0' );
%     %             end

% return

%% Problem 2e

% return

%% Problem 2f

% return

%% Problem 2g

% return

%% 2 DOF Example

% masses = [ 50  1 ];
% stiffnesses = [ 100e3  1800  0 ];
% 
% dampings = [ 0  0  0 ];
% % dampings = [ 0  (1 + 0.1*1j)  (1 + 0.1*1j) ];
% 
% frequencies = 0:0.01:40;
% 
% FRF = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );  % 4001-by-2-by-2
% 
% 
% admittance = zeros( numel( frequencies ), 1 );
% 
% for index = 1:1:numel( frequencies )
%     temp = diag( squeeze( FRF( index, :, : ) ) );
%         temp2 = temp(1) + temp(2)*1j;
%         admittance( index ) = abs( temp2 );
% end
% 
% clear temp1 temp2;
% 
% 
% % figure( ); ...
% %     semilogy( frequencies, admittance );  grid on;
% %     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
% 
% 
% figure( 'Name', 'Admittance' ); ...
%     semilogy( frequencies, abs( FRF( :, 1, 1 ) ) );  hold on;
%     semilogy( frequencies, abs( FRF( :, 1, 2 ) ) );
%     semilogy( frequencies, abs( FRF( :, 2, 1 ) ) );
%     semilogy( frequencies, abs( FRF( :, 2, 2 ) ) );  grid on;
%         legend( 'M1 Movement with F1', 'M1 Movement with F2', 'M2 Movement with F1', 'M2 Movement with F2', 'Location', 'North' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
%     axis( [ 5.5 10.5  1e-7 1e2 ] );


% figure( 'Name', 'Dynamic Stiffness (Reciprocal of Admittance)' ); ...
%     semilogy( frequencies, 1 ./ abs( FRF( :, 1, 1 ) ) );  hold on;
%     semilogy( frequencies, 1 ./ abs( FRF( :, 1, 2 ) ) );
%     semilogy( frequencies, 1 ./ abs( FRF( :, 2, 1 ) ) );
%     semilogy( frequencies, 1 ./ abs( FRF( :, 2, 2 ) ) );  grid on;
%         legend( 'M1 Movement with F1', 'M1 Movement with F2', 'M2 Movement with F1', 'M2 Movement with F2', 'Location', 'North' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );



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


