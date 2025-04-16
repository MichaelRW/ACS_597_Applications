


%% Synopsis

% Problem 2 - Active Control of Enclosure Opening

% See Lecture 22
%   ""D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\22 Monday, April 7, 2025\Lecture 22 - Coherence effects - filled.pptx""



%% To Do

% ph



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

h_sound_power_ratio = @( k, d )  ( 4/5 ) .* ( k .* d ).^4;

h_pressure_ratio = @( theta, sound_level_ratio )  5 .* cos( theta ).^4 * sound_level_ratio;



%% Constants and Parameters

c = 343;  % m/s
rho0 = 1.21;  % kg/m^3

frequencies = [ 63  125  250  500 ];  % Hz
    wavelengths = c ./ frequencies;  % 5.44 m, 2.74 m, 1.37 m, and 0.69 m

k = (2*pi*frequencies) ./ c;  % 1.15 1/m, 2.29 1/m, 4.58 1/m, and 9.16 1/m


color_map = slanCM( 'ColorBlind', 4 );



%% Problem 2a

d = 205e-3;  % m;  distance between each of the two speakers and the center of the square


figure( 'Name', 'Sound Power Level Reductions' ); ...

    plot( frequencies, 10*log10( h_sound_power_ratio( k, d ) ) );  grid on;
        xlabel( 'Octave Band Center Frequency [Hz]' );  ylabel( 'Sound Power Level Reduction [dB]' );
    xticks( [ 63, 125, 250, 500 ] );  xticklabels( { '63 Hz', '125 Hz', '250 Hz', '500 Hz' } );
    axis( [ 50 550 -30 12 ] );



%% Problem 2b

theta_degrees = 0:0.5:90;
    theta_radians = theta_degrees.*pi./180;


figure( 'Name', 'Sound Power Level Reductions Versus Theta' ); ...

    plot( theta_degrees, 10*log10( h_pressure_ratio( theta_radians, h_sound_power_ratio( k(1), d ) ) ) );  hold on;
    plot( theta_degrees, 10*log10( h_pressure_ratio( theta_radians, h_sound_power_ratio( k(2), d ) ) ) );
    plot( theta_degrees, 10*log10( h_pressure_ratio( theta_radians, h_sound_power_ratio( k(3), d ) ) ) );
    plot( theta_degrees, 10*log10( h_pressure_ratio( theta_radians, h_sound_power_ratio( k(4), d ) ) ) );  grid on;
        legend( '63 Hz', '125 Hz', '250 Hz', '500 Hz', 'Location', 'SouthWest' )
    xlabel( 'Octave Band Center Frequency [Hz]' );  ylabel( 'Sound Power Level Reduction [dB]' );
    axis( [ -5 95 -1e2 25 ] );


    
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


