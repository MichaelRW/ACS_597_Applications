


%% Synopsis

% Problem 2 - Simulation of Sources from Monopoles



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

set( groot, 'DefaultFigurePosition', [ 100  750    750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format LongG;

pause( 1 );



%% Anonymous Function Definitions



%% Strikeforce

c = 343;  % m\s
rho0 = 1.21;  % kg/m^3

f = 1e3;
    lambda = c / f;  % 0.343 m

fractions = [ 0.5  0.75  1 1.25  1.5  2  2.5  3  3.5  5  10 ];



xyz_sources = [ 0, 0, 0 ];
Q_sources = 1 + 1i*1;



%% Problem 2a

x = lambda .* fractions;
    xyz_receivers = [ x.'  zeros( numel(x) , 1 )  zeros( numel(x) , 1 ) ];

p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, f, rho0, c );

% figure( ); ...
%     loglog( x, abs(p) );  grid on;
%     xlabel( 'Distance [m]' );  ylabel( 'Pressure [dB]' );



%% Problem 2b

% Concentric Circle 1
N = 180;
    n = 0:1:N;

temp = exp( 1i.*(2*pi)/N*n);

    
temp2 = temp * fractions( 1 );
polarplot( angle(temp2), abs(temp2) );
hold on;
for index = 2:1:numel( fractions )
    temp2 = temp * fractions( index );
    polarplot( angle(temp2), abs(temp2), 'Color', 'k' );
end
grid on;


polarplot3d



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


