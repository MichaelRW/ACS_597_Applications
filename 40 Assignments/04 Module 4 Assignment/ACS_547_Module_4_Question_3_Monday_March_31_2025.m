


%% Synopsis

% Problem 3 - Coherent and Incoherent Speakers

% See Lecture 22 on Monday, April 7, 2025
%   ""D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\22 Monday, April 7, 2025\Lecture 22 - Coherence effects - filled.pptx""



%% Note(s)

% Distance (source and receivers) are in units of meters.

% Source can be real-valued or complex-valued.  Complex-valued sources have magnitude and phase information.

% Sources have units of volume velocity, m^3/s.

% Complex pressures are in units of Pascals.



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

set( groot, 'DefaultFigurePosition', [ 1.7e3  775    750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'norma' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Define Anonymous Function(s)

h_volume_velocity = @( rho, c, k, watts )  sqrt( ( watts*8*pi ) / (rho*c*k^2 ) );



%% Constants and Parameters

c = 343;  % m/s
rho0 = 1.21;  % kg/m^3

frequency = 125;  % Hz
    wavelengths = c ./ frequency;  % 2.74 m

k = (2*pi*frequency) ./ c;  % 2.29 1/m



%% Problem 2a

d = 1e-1;  % m


% Calculate volume velocity (see Lecture 19, slide 35).
volume_velocity_0_5_watts = h_volume_velocity ( rho0, c, k, 0.5 );
volume_velocity_1_0_watts = h_volume_velocity ( rho0, c, k, 1 );


xyz_sources = [ ... 
    +d, 0, 0; ...
    0, 0, 0; ...
    -d, 0, 0; ...
    ];  % m
%
Q_sources = [ ...
    +volume_velocity_0_5_watts; ...
    -volume_velocity_1_0_watts; ...
    +volume_velocity_0_5_watts; ...
    ];  % m^3/s

xyz_receiver = [ 20  tand(30)*20  1 ];
% xyz_receiver = [ 0  0  1 ];

p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receiver, frequency, rho0, c );
    L = 10*log10( abs(p).^2 );  % -17.3 dB



%% Problem 2b

xyz_sources = [ ... 
    +d, 0, 0; ...
    0, 0, 0; ...
    -d, 0, 0; ...
    ];  % m
%
Q_sources = [ ...
    +volume_velocity_0_5_watts; ...
    -volume_velocity_1_0_watts; ...
    +volume_velocity_0_5_watts; ...
    ];  % m^3/s


% Speaker 1
p1 = sum_of_monopoles( xyz_sources(1, :), Q_sources(1), xyz_receiver, frequency, rho0, c );
    L1 = 10*log10( abs( p1 ) );  % -6.0 dB

% Speaker 2
p2 = sum_of_monopoles( xyz_sources(2, :), Q_sources(2), xyz_receiver, frequency, rho0, c );
    L2 = 10*log10( abs( p2 ) );  % -4.5 dB

% Speaker 3
p3 = sum_of_monopoles( xyz_sources(3, :), Q_sources(3), xyz_receiver, frequency, rho0, c );
    L3 = 10*log10( abs( p3 ) );  % -6.1 dB


L_incoherent_sources = 10*log10(  10^(L1/10) + 10^(L2/10) + 10^(L3/10 ) );  % -0.71 dB


    
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

p = 10 + 1i*5;  % Pascals;  sinusoid

p_mag = abs( p );  % 11.18 Pascals

p_rms = p_mag / sqrt(2);  % 7.91 Pascals RMS

p_dB_SPL = 20*log10( p_rms / 20e-6 );  % 111.94 dB SPL Z


p_dB_SPL_verify = convert_complex_pressure_to_dB_SPL( p );




%% URLs


