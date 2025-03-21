


%% Synopsis

% Problem 2 - Washing Machine Two-stage Mount Design



%% Note(s):

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



%% Problem 2a

% Use Visio to draw diagram of the two-stage mount.  Include dash-pots.



%% Problem 2b

% The frequency of the forcing function is 5 Hz and the amplitude is 494.4 N.
fo = 5;  % Hz
    wo = 2*pi*fo;  % 31.4 radians\s


m1 = 100 + 10;  % kg - Total mass of washing machine and its load.
    k1 = 5e3;  % N\m
    c1 = 10;

m2 = 100;  % kg - UPPER LIMIT
    k2 = 6e3;  % N\m
    c2 = 10;

k3 = 100e3;


w1 = sqrt( ( k1 + k3 ) / m1 );  % 30.9 radians\s

w2 = sqrt( ( k2 + k3 ) / m2 );  % 32.6 radians\s

% return

%% Problem 2c

masses = [ m2 m1 ];
stiffnesses = [ k2  k3  k1 ];
dampings = [ c2  0  c1 ];

% frequencies = 0:0.01:40;

rpm = 0:1:500;  % rotations per minute
    rpm_conversion_to_radians_per_second = 0.10472;  % radians\s
        angular_velocity = rpm .* rpm_conversion_to_radians_per_second;
            frequencies = angular_velocity / (2*pi);

FRF = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );


admittance = zeros( numel( frequencies ), 1 );

for index = 1:1:numel( frequencies )
    temp = diag( squeeze( FRF( index, 2, 2 ) ) );
        admittance( index ) = abs( temp );
end
%
clear temp temp2;

washing_machine.admittance = admittance;


admittance = zeros( numel( frequencies ), 1 );

for index = 1:1:numel( frequencies )
    temp = diag( squeeze( FRF( index, 2, 1 ) ) );
        admittance( index ) = abs( temp );
end
%
clear temp temp2;

raft.admittance = admittance;


figure( 'Name', '' ); ...
    h1 = loglog( rpm, washing_machine.admittance );  hold on;
    h2 = loglog( rpm, raft.admittance );
    h3 = line( [ 300 300 ], [ 10^-9 10^-2 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
        legend( [ h1 h2 h3 ], 'Washing Machine', 'Raft', 'Load Frequency', 'Location', 'SouthWest' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
    % axis( [ 0 max(frequencies)  1e-9 1e-2 ] );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_2c', '-dpdf', '-r0' );
    %             end


% The admittance at 5 Hz (300 RPM) is 7.29e-7 m\N, or 2,085 times smaller
% than it is for the first-order system.



%% Problem 2d

dynamic_force.mass = 1;  % kg
dynamic_force.radius = 0.5;  % m


% rpm = 0:1:350;  % rotations per minute
%     rpm_conversion_to_radians_per_second = 0.10472;  % radians\s
%         angular_velocity = rpm .* rpm_conversion_to_radians_per_second;
%             frequencies = angular_velocity / (2*pi);

h_force = @( force_mass, rotation_speed, load_distance  )  force_mass .* rotation_speed.^2 .* load_distance;


temp = h_force( dynamic_force.mass, angular_velocity, dynamic_force.radius ).';


figure( 'Name', '' ); ...
    h1 = loglog( rpm, washing_machine.admittance.*temp );  hold on;
    h2 = loglog( rpm, raft.admittance.*temp );
    h3 = line( [ 300 300 ], [ 10^-9 10^-2 ], 'Color', 'k', 'LineStyle', '--' );  grid on;
        legend( [ h1 h2 h3 ], 'Washing Machine', 'Raft', 'Load Frequency', 'Location', 'SouthWest' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Displacment [m]' );
    % axis( [ 0 max(frequencies)  1e-9 1e-2 ] );
    % 
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             if ( PRINT_FIGURES == 1 )
    %                 print(gcf, 'Homework_3_Part_2c', '-dpdf', '-r0' );
    %             end

return

%% Problem 2e

return

%% Problem 2f

return

%% Problem 2g

return

%% 1 DOF Example

masses = 50;
stiffnesses = [ 100e3  0 ];
dampings = [ 0  0 ];
    % dampings = [ 0  (1 + 0.1*1j) ];

wo = sqrt( 100e3 / 50 );
    fo = wo/(2*pi);  % 7.1 Hz

frequencies = 0:0.1:100;


admittance = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );  % 4001-by-2-by-2
%
figure( ); ...
    loglog( frequencies, abs( admittance ) );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );


impedance_original_function = nDOF_direct_solution_org( masses, stiffnesses, dampings, frequencies, 'impedance' );  % 4001-by-2-by-2
impedance_updated_function = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'impedance' );  % 4001-by-2-by-2
%
figure( ); ...
    loglog( frequencies, abs( impedance_original_function ) );  hold on;
    loglog( frequencies, abs( impedance_updated_function ) );  grid on;
        legend( 'Original Function', 'Updated Function', 'Location', 'North' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Impedance [$\frac{N\cdots}{m}$]' );



%% 2 DOF Example

masses = [ 50  1 ];
stiffnesses = [ 100e3  1800  0 ];

dampings = [ 0  0  0 ];
% dampings = [ 0  (1 + 0.1*1j)  (1 + 0.1*1j) ];

frequencies = 0:0.01:40;

FRF = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );  % 4001-by-2-by-2


admittance = zeros( numel( frequencies ), 1 );

for index = 1:1:numel( frequencies )
    temp = diag( squeeze( FRF( index, :, : ) ) );
        temp2 = temp(1) + temp(2)*1j;
        admittance( index ) = abs( temp2 );
end

clear temp1 temp2;


% figure( ); ...
%     semilogy( frequencies, admittance );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );


figure( 'Name', 'Admittance' ); ...
    semilogy( frequencies, abs( FRF( :, 1, 1 ) ) );  hold on;
    semilogy( frequencies, abs( FRF( :, 1, 2 ) ) );
    semilogy( frequencies, abs( FRF( :, 2, 1 ) ) );
    semilogy( frequencies, abs( FRF( :, 2, 2 ) ) );  grid on;
        legend( 'M1 Movement with F1', 'M1 Movement with F2', 'M2 Movement with F1', 'M2 Movement with F2', 'Location', 'North' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
    
figure( 'Name', 'Dynamic Stiffness' ); ...
    semilogy( frequencies, 1 ./ abs( FRF( :, 1, 1 ) ) );  hold on;
    semilogy( frequencies, 1 ./ abs( FRF( :, 1, 2 ) ) );
    semilogy( frequencies, 1 ./ abs( FRF( :, 2, 1 ) ) );
    semilogy( frequencies, 1 ./ abs( FRF( :, 2, 2 ) ) );  grid on;
        legend( 'M1 Movement with F1', 'M1 Movement with F2', 'M2 Movement with F1', 'M2 Movement with F2', 'Location', 'North' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );



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


