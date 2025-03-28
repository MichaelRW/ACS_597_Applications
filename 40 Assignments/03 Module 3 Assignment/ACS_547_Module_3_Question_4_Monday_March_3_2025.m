


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



%% Impulsive Perturbation

time_step = 1e-3;  % s
net_time = 4;  % s
    time_indices = 0:time_step:( net_time - time_step );

impulse = zeros( size( time_indices ) );
    impulse( 1/time_step ) = 1;


figure( 'Name', 'Impulse Signal' ); ...
    stem( time_indices, impulse, 'Marker', '.' );  grid on;
        legend( 'Impulse Perturbation', 'Location', 'NorthEast' );
    xlabel( 'Time [s]' );  ylabel( 'Amplitude [WU]' );
    axis( [ 0 net_time -0.1 1.5 ] );



%% System Impulse Response

m = 5;  % kg

% Initial conditions of 1 DOF system.
xo = 0;  % m
vo = 1 / m;  % m\s

wd = 31.4;  % radians\s or 5 Hz or 300 RPM


ks = 9810;  % N\m
    wo = sqrt( ks / 110 );  % 9.4 radians\s
        fo = wo / (2*pi);  % 1.5 Hz


epsilon = 0.25;
    C = epsilon * 2*sqrt( ks * 110 );  % 519.4 kg\s


h_unit_impulse = @( total_mass, wd, wo, epsilon, t )  ( 1./( total_mass.*wd) )  .*  exp( -wo.*epsilon.*t ) .* sin( wd.*t );



%%

impulse_response = h_unit_impulse( 110, wd, wo, epsilon, time_indices );
%
figure( ); ...
    plot( time_indices, impulse_response );  grid on;


response_to_brick_perturbation = conv( impulse_response, 5*impulse );
    time_indices_perturbation = ( 0:1:( numel( response_to_brick_perturbation ) - 1 ) ) .* time_step;
%
figure( ); ...
    plot( time_indices_perturbation, response_to_brick_perturbation );  grid on;
    xlim( [ 0 net_time ] );



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


