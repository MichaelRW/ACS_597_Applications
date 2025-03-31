


%% Synopsis

% Problem 4 - Washing Machine Abuse Testing



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

h_unit_impulse = @( washing_machine_mass, wd, wo, epsilon, t )  ( 1./( washing_machine_mass.*wd) )  .*  exp( -wo.*epsilon.*t ) .* sin( wd.*t );  % mvo = 1 kg\(m s)



%% Impulse Signal

time_step = 1e-4;  % s
net_time = 10;  % s
    time_indices = 0:time_step:( net_time - time_step );

impulse = zeros( size( time_indices ) );
    impulse( 4/time_step ) = 1./time_step;
        brick_impulse = 5 .* impulse;



%% 1 DOF System Impulse Function

ks = 9810;  % N\m
    wo = sqrt( ks / 110 );  % 9.4 radians\s
        fo = wo / (2*pi);  % 1.5 Hz

epsilon = 0.25;
    C = epsilon * 2*sqrt( ks * 110 );  % 519.4 kg\s

wd = wo * sqrt( 1 - epsilon^2 );  % 9.1 radians\s
    fd = wd / (2*pi);  % 1.46 Hz



%% 1 DOF System Impulse Response

impulse_response = h_unit_impulse( 100, wd, wo, epsilon, time_indices );

impulse_response_at_4 = conv( impulse, impulse_response ) * time_step;
    time_indices_at_4 = ( 0:1:( numel( impulse_response_at_4 ) - 1 ) ) .* time_step;
%
% figure( 'Name', '1 DOF System Impulse Response' ); ...
%     yyaxis left; ...
%         plot( time_indices, impulse );  grid on;
%         ylabel( 'Force [N]' );
%         ylim( 'auto' );
%     yyaxis right; ...
%         plot( time_indices_at_4, impulse_response_at_4 );  hold on;
%         ylabel( 'Admittance [$\frac{m}{N}$]' );
%         ylim( 'auto' );
%     xlabel( 'Time [s]' );
%     xlim( [ 0 10 ] );



%% Combined Response

sinusoid_input = 493.5 * sin( 2 * pi * 5 * time_indices );
    sinusoid_input( 1:1:( 1 / time_step ) ) = 0;
    sinusoid_input( ( 8 / time_step ):1:( 10 / time_step ) ) = 0;


joint_signal = sinusoid_input + brick_impulse;


response_to_combined_forcing = conv( impulse_response, joint_signal ) .* time_step;
    time_indices_sinusoidal_forcing = ( 0:1:( numel( response_to_combined_forcing ) - 1 ) ) .* time_step;
%
figure( 'Name', '1 DOF Combined Forcing Response' ); ...
    h1 = plot( time_indices_sinusoidal_forcing, response_to_combined_forcing );  hold on;
    %
    h2 = line( [ 1 1 ], [ -2e-2 2e-2 ], 'Color', [ 0.47, 0.67, 0.19 ], 'LineStyle', '-.' );
    
    h3 = line( [ 4 4 ], [ -2e-2 2e-2 ], 'Color', [ 0.47, 0.67, 0.19 ], 'LineStyle', '--' );
    h4 = line( [ 6.3 6.3 ], [ -2e-2 2e-2 ], 'Color', 'r', 'LineStyle', '--' );
    % 
    h5 = line( [ 8 8 ], [ -2e-2 2e-2 ], 'Color', 'r', 'LineStyle', '-.' );  grid on;

    h6 = line( [4 8 ], [ 5.42e-3 5.42e-3 ], 'LineStyle', ':', 'Color', 'k' );
        
        legend( [ h1 h2 h3 h4 h5 h6 ], 'Response', 'Start Sinusoidal Forcing', ...
            'Brick Inserted', 'End of Brick Response', ...
            'Stop Sinusoidal Forcing', 'Impulse Response Measure', ...
            'Location', 'NorthEast', 'Interpreter', 'Latex' );

    xlabel( 'Time [s]' );  ylabel( 'Displacement [m]' );
    axis( [ 0 12    -2e-2 2e-2  ] );

return

%% Acceleration

velocity = movingslope( response_to_combined_forcing );
velocity2 = gradient( response_to_combined_forcing );
%
figure( 'Name', 'Velocity' ); ...
    plot( time_indices_sinusoidal_forcing, velocity );  hold on;
    plot( time_indices_sinusoidal_forcing, velocity2, 'LineStyle', '--' );  grid on;
    xlabel( 'Time [s]' );  ylabel( 'Velocity [$\frac{m}{s}$]' );
    % xlim( [ 0 2.5 ] );


acceleration = movingslope( velocity );
acceleration2 = gradient( velocity2 );
%
figure( 'Name', 'Acceleration' ); ...
    plot( time_indices_sinusoidal_forcing, acceleration );  hold on;
    plot( time_indices_sinusoidal_forcing, acceleration2, 'LineStyle', '--' );  grid on;
    xlabel( 'Time [s]' );  ylabel( 'Acceleration [$\frac{m}{s^2}$]' );
    % xlim( [ 0 2.5 ] );



%% Verify

close all;

t = time_indices_sinusoidal_forcing;
    x = sin( 2 * pi * 5 * t );

time_step_verify = t(2) - t(1);

figure( ); ...
    plot( t, x );  hold on;
    % plot( t(1:1:end-1), diff( x ) / time_step_verify );
    plot( t, gradient( x ) / time_step_verify );
    grid on;
    % axis( [ 0 1  -1.2 1.2 ] );



%% Plot Difference

% figure( 'Name', '1 DOF Combined Forcing Response' ); ...  % Figure 7
%     plot( time_indices_sinusoidal_forcing, response_to_combined_forcing - response_to_sinusoidal_forcing );  grid on;
%     xlabel( 'Time [s]' );  ylabel( 'Displacement [m]' );
%     % xlim( [ 0 2.5 ] );



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


