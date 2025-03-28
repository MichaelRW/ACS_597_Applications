


%% Synopsis

% Problem 4 - Washing Machine Abuse Testing



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



%% Impulse Signal

time_step = 1e-2;  % s
net_time = 5;  % s
    time_indices = 0:time_step:( net_time - time_step );

impulse = zeros( size( time_indices ) );
    impulse( 1/time_step ) = 1./time_step;
        impulse = impulse .* 5;


figure( 'Name', 'Unit Impulse Signal' ); ...  % Figure 1
    stem( time_indices, impulse, 'Marker', '.' );  grid on;
        legend( 'Unit Impulse Signal', 'Location', 'NorthEast' );
    xlabel( 'Time [s]' );  ylabel( 'Force [N]' );
    xlim( [ 0 net_time ] );
    ylim( [ -600 600 ] );



%% 1 DOF System Impulse Function

m = 5;  % kg
m = 10;  % kg

% Initial conditions of 1 DOF system.
xo = 0;  % m
vo = 1 / m;  % m\s

wd = 31.4;  % radians\s or 5 Hz or 300 RPM


ks = 9810;  % N\m
    wo = sqrt( ks / 110 );  % 9.4 radians\s
        fo = wo / (2*pi);  % 1.5 Hz


epsilon = 0.25;
    C = epsilon * 2*sqrt( ks * 110 );  % 519.4 kg\s


h_unit_impulse = @( total_mass, wd, wo, epsilon, t )  ( 1./( total_mass.*wd) )  .*  exp( -wo.*epsilon.*t ) .* sin( wd.*t );  % mvo = 1 kg\(m s)



%% 1 DOF System Impulse Response

impulse_response = h_unit_impulse( 110, wd, wo, epsilon, time_indices );
%
figure( 'Name', '1 DOF System Impulse Response Function' ); ...  % Figure 2
    plot( time_indices, impulse_response );  grid on;
    xlabel( 'Time [s]' );  ylabel( 'Amplitude [WU]' );


response_to_brick_perturbation = conv( impulse_response, impulse );
    time_indices_perturbation = ( 0:1:( numel( response_to_brick_perturbation ) - 1 ) ) .* time_step;
%
figure( 'Name', '1 DOF Impulse Response' ); ...  % Figure 3
    plot( time_indices_perturbation, response_to_brick_perturbation );  grid on;
    xlabel( 'Time [s]' );  ylabel( 'Displacement [m]' );



%% Sinusoidal Impulse Response

% last_index = (1/5) / time_step;
last_index = ( 1.1 * 0.2 ) / time_step;

sinusoid_input = 493.5/2 * sin( 2 * pi * 5 * time_indices - 1.1*pi );
    sinusoid_input( 1:1:last_index ) = 0;
%
figure( 'Name', 'Sinusoidal Input' ); ... % Figure 4
    plot( time_indices, sinusoid_input );  grid on;
    xlabel( 'Time [s]' );  ylabel( 'Force [N]' );
    xlim( [ 0 2.5 ] );
    ylim( [ -600 600 ] );


response_to_sinusoidal_forcing = conv( impulse_response, sinusoid_input );
    time_indices_sinusoidal_forcing = ( 0:1:( numel( response_to_brick_perturbation ) - 1 ) ) .* time_step;
%
figure( 'Name', '1 DOF Sinusoidal Forcing Response' ); ...  % Figure 5
    plot( time_indices_sinusoidal_forcing, response_to_sinusoidal_forcing );  grid on;
    xlabel( 'Time [s]' );  ylabel( 'Displacement [m]' );
    xlim( [ 0 2.5 ] );



%% Combined Response

joint_signal = sinusoid_input + impulse;
%
figure( 'Name', 'Combined Signal' ); ... % Figure 6
    plot( time_indices, joint_signal );  hold on;
    plot( time_indices, sinusoid_input, 'LineStyle', '--' );
    plot( time_indices, impulse, 'LineStyle', '--' );
    xlabel( 'Time [s]' );  ylabel( 'Force [N]' );
    xlim( [ 0 2.5 ] );
    ylim( [ -600 600 ] );



response_to_combined_forcing = conv( impulse_response, joint_signal );
    time_indices_sinusoidal_forcing = ( 0:1:( numel( response_to_brick_perturbation ) - 1 ) ) .* time_step;
%
figure( 'Name', '1 DOF Combined Forcing Response' ); ...  % Figure 7
    plot( time_indices_sinusoidal_forcing, response_to_combined_forcing );  grid on;
    xlabel( 'Time [s]' );  ylabel( 'Displacement [m]' );
    xlim( [ 0 2.5 ] );
        



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


