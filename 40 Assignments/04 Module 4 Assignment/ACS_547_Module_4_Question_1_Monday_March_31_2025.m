


%% Synopsis

% Problem 1 - Simulation of Sources from Monopoles

% See Lecture 20
%   "D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\20 Monday, March 31, 2025\Lecture 20 - Dipoles and quadrupoles - Filled.pptx"



%% To Do

% Confirm relationship of pressure dependence on distance for a given frequency.

% https://www.google.com/search?q=wiki+wave+number&rlz=1C1UEAD_enCA1080CA1080&oq=wiki+wave+number&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIICAEQABgWGB4yBwgCEAAY7wUyCggDEAAYgAQYogQyBggEEEUYPNIBCDI0NDNqMGo3qAIAsAIA&sourceid=chrome&ie=UTF-8
% https://en.wikipedia.org/wiki/Wavenumber



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



%% Constants and Parameters

c = 343;  % m/s
rho0 = 1.21;  % kg/m^3

frequency = 1e3;  % Hz
    wavelength = c / frequency;  % 0.343 m

k = (2*pi*frequency) ./ c;  % 18.3 1/m (alternative calculation:  2pi / wavelength)
        
radius_fractions = [ 0.5  2.4  4.8  9.6 ];


r = radius_fractions;  % m (unity multiplier)

theta = ( 0:0.01:2*pi ).';

xyz_receivers = nan( numel( r ), numel( theta), 3 );  % m

for radius_index = 1:1:numel( r )
    x = r( radius_index )*cos( theta );  y = r( radius_index )*sin( theta );  z = zeros( size( x ) );
        xyz_receivers( radius_index, :, :, : ) = [ x(:) y(:) z(:) ];
end


color_map = slanCM( 'ColorBlind', 4 );



%% Problem 2a - Pressure Versus Distance for a Monopole Source

% Monopole sound pressure verus distance has a decay of 6 dB per doubling of distance.

% xyz_sources = [ ...
%     0, 0, 0; ...
%     ];
% %
% Q_sources = 1;  % m^3/s
% 
% 
% x = 0.01:0.01:10;
%     y = zeros( length( x ), 1 );  z = zeros( length( x ), 1 );
%         xyz_receivers = [ x(:) y(:) z(:) ];
% 
% 
% figure( 'Name', 'Monopole Source - Pressure Magnitude Versus Distance' ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, frequency, rho0, c );
%         L = 10*log10( abs(p).^2 );
%                 plot( x, L, 'Color', [ color_map( 1, : ) 0.8 ] );  grid on;
%     legend( '1 kHz', 'Location', 'NorthEast' );
%     xlabel( 'Distance [m]' );  ylabel( 'Pressure Magnitude [$dB$]' );
%     set( gca, 'XScale', 'log' );


% Note(s):

%   With a higher frequency, the wave number, "k" will be larger.  Pressure will be higher at the same distance.
%   For a doubling in frequency, the pressure will be 20 dB larger.



%% Problem 2b - Directivity Pattern for a Monopole Source

% xyz_sources = [ ... 
%     0, 0, 0; ...
%     ];  % m
% %
% Q_sources = 1.4;  % m^3/s
% 
% 
% xyz_receivers = [ cos( theta )  sin( theta )  zeros( size( cos( theta ) ) ) ];
% 
% figure( 'Name', 'Monopole Source - Directivity Pattern' ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, frequency, rho0, c );
%         monopole.L = 10*log10( abs(p).^2 );
% 
%     polarplot( theta, monopole.L, 'Color', [ color_map( 1, : ) 0.8 ] );
%         tick_marks = monopole.L( 1 );
% 
%     monopole.label = sprintf( '1 kHz, Q = 1.4 %s @ %1g m, %3.1f dB', '$\frac{m^3}{s}$', string( [ 1  round( tick_marks, 1 ) ] ) );
% 
%     rlim( [ 0  1e2 ] );
%         rticklabels( sprintf( '1 kHz, @ %1g m, %3.1f dB', string( [ 1  round( tick_marks, 1 ) ] ) ) );  rtickangle( 45 );
%         rticks( tick_marks );


% figure( 'Name', 'Monopole Source - Pressure Versus Angle' ); ...
% 
%     plot( theta.*180./pi, monopole.L, 'Color', [ color_map( 1, : ) 0.8 ] );  grid on;
%     legend( monopole.label, 'Location', 'North', 'Interpreter', 'Latex' );
%     xlabel( 'Angle [Degree]' );  ylabel( 'Pressure at 1 m [dB]' );
%     axis( [ -10 370  0 80 ] );



%% Problem 2b - Directivity Pattern for a Dipole Source

% d = wavelength / 8;  % 0.043 m
%     k*d;  % 0.79
% 
% xyz_sources = [ ... 
%     -d, 0, 0; ...
%     +d, 0, 0; ...
%     ];  % m
% %
% Q_sources = [ ...
%     +1.4; ...
%     -1.4; ...
%     ];  % m^3/s
% 
% 
% xyz_receivers = [ cos( theta )  sin( theta )  zeros( size( cos( theta ) ) ) ];
% 
% figure( 'Name', 'Dipole Source - Directivity Pattern' ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, frequency, rho0, c );
%         dipole.L = 10*log10( abs(p).^2 );
% 
%     polarplot( theta, dipole.L, 'Color', [ color_map( 2, : ) 0.8 ] );
%         tick_marks = max( dipole.L );
% 
%     dipole.label = sprintf( '1 kHz, Q = 1.4 %s @ %1g m, %3.1f dB', '$\frac{m^3}{s}$', string( [ 1  round( tick_marks, 1 ) ] ) );
% 
%     rlim( [ 0  80 ] );
%         rticklabels( sprintf( '1 kHz, @ %1g m, %3.1f dB', string( [ 1  round( tick_marks, 1 ) ] ) ) );  rtickangle( 45 );
%         rticks( tick_marks );


% figure( 'Name', 'Dipole Source - Pressure Versus Angle' ); ...
% 
%     plot( theta.*180./pi, dipole.L, 'Color', [ color_map( 2, : ) 0.8 ] );  grid on;
%     legend( dipole.label, 'Location', 'North', 'Interpreter', 'Latex' );
%     xlabel( 'Angle [Degree]' );  ylabel( 'Pressure at 1 m [dB]' );
%     axis( [ -10 370  0 80 ] );



%% Problem 2b - Lateral Quadrupole Directivity Pattern

% d = 1e-1;  % m
% 
% xyz_sources = [ ... 
%     +d, +d, 0; ...
%     -d, +d, 0; ...
%     -d, -d, 0; ...
%     +d, -d, 0; ...
%     ];  % m
% %
% Q_sources = [ ...
%     +1.4; ...
%     -1.4; ...
%     +1.4; ...
%     -1.4; ...
%     ];  % m^3/s
% 
% 
% xyz_receivers = [ cos( theta )  sin( theta )  zeros( size( cos( theta ) ) ) ];
% 
% figure( 'Name', 'Lateral Quadrupole Source - Directivity Pattern' ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, frequency, rho0, c );
%         lateral_quadrupole.L = 10*log10( abs(p).^2 );
% 
%     polarplot( theta, lateral_quadrupole.L, 'Color', [ color_map( 3, : ) 0.8 ] );
%         tick_marks = max( lateral_quadrupole.L );
% 
%     lateral_quadrupole.label = sprintf( '1 kHz, Q = 1.4 %s @ %1g m, %3.1f dB', '$\frac{m^3}{s}$', string( [ 1  round( tick_marks, 1 ) ] ) );
% 
%     rlim( [ 0  80 ] );
%         rticklabels( sprintf( '1 kHz, @ %1g m, %3.1f dB', string( [ 1  round( tick_marks, 1 ) ] ) ) );  rtickangle( 45 );
%         rticks( tick_marks );


% figure( 'Name', 'Lateral Quadrupole Source - Pressure Versus Angle' ); ...
% 
%     plot( theta.*180./pi, lateral_quadrupole.L, 'Color', [ color_map( 3, : ) 0.8 ] );  grid on;
%     legend( lateral_quadrupole.label, 'Location', 'North', 'Interpreter', 'Latex' );
%     xlabel( 'Angle [Degree]' );  ylabel( 'Pressure at 1 m [dB]' );
%     axis( [ -10 370  0 100 ] );



%% Problem 2b - Linear Quadrupole Directivity Pattern

% d = wavelength / 8;  % m
% 
% xyz_sources = [ ... 
%     +d, 0, 0; ...
%     0, 0, 0; ...
%     -d, 0, 0; ...
%     ];  % m
% %
% Q_sources = [ ...
%     -1.4; ...
%     +2.8; ...
%     -1.4; ...
%     ];  % m^3/s
% 
% 
% xyz_receivers = [ cos( theta )  sin( theta )  zeros( size( cos( theta ) ) ) ];
% 
% figure( 'Name', 'Linear Quadrupole Source - Directivity Pattern' ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, frequency, rho0, c );
%         linear_quadrupole.L = 10*log10( abs(p).^2 );
% 
%     polarplot( theta, linear_quadrupole.L, 'Color', [ color_map( 4, : ) 0.8 ] );
%         tick_marks = max( linear_quadrupole.L );
% 
%     linear_quadrupole.label = sprintf( '1 kHz, Q = 1.4 %s @ %1g m, %3.1f dB', '$\frac{m^3}{s}$', string( [ 1  round( tick_marks, 1 ) ] ) );
% 
%     rlim( [ 0  80 ] );
%         rticklabels( sprintf( '1 kHz, @ %1g m, %3.1f dB', string( [ 1  round( tick_marks, 1 ) ] ) ) );  rtickangle( 45 );
%         rticks( tick_marks );


% figure( 'Name', 'Linear Quadrupole Source - Pressure Versus Angle' ); ...
% 
%     plot( theta.*180./pi, linear_quadrupole.L, 'Color', [ color_map( 4, : ) 0.8 ] );  grid on;
%     legend( linear_quadrupole.label, 'Location', 'North', 'Interpreter', 'Latex' );
%     xlabel( 'Angle [Degree]' );  ylabel( 'Pressure at 1 m [dB]' );
%     axis( [ -10 370  0 100 ] );



%% Problem 2b - Combined Plot

% figure( 'Name', 'Combined Directivity Patterns - Polar Plot' ); ...
% 
%     h_1 = polarplot( theta, monopole.L, 'Color', [ color_map( 1, : ) 0.8 ] );  hold on;
%     h_2 = polarplot( theta, dipole.L, 'Color', [ color_map( 2, : ) 0.8 ] );
%     h_3 = polarplot( theta, lateral_quadrupole.L, 'Color', [ color_map( 3, : ) 0.8 ] );
%     h_4 = polarplot( theta, linear_quadrupole.L, 'Color', [ color_map( 4, : ) 0.8 ] );
%         legend( [ h_1, h_2, h_3, h_4 ], { 'Monopole', 'Dipole', 'Lateral Quadrupole', 'Linear Quadrupole' }, 'Location', 'EastOutside' );
%     %
%     rlim( [ 0  80 ] );


% figure( 'Name', 'Combined Directivity Patterns - Pressure Versus Angle' ); ...
% 
%     h_1 = plot( theta.*180./pi, monopole.L, 'Color', [ color_map( 1, : ) 0.8 ] );  hold on;
%     h_2 = plot( theta.*180./pi, dipole.L, 'Color', [ color_map( 2, : ) 0.8 ] );
%     h_3 = plot( theta.*180./pi, lateral_quadrupole.L, 'Color', [ color_map( 3, : ) 0.8 ] );
%     h_4 = plot( theta.*180./pi, linear_quadrupole.L, 'Color', [ color_map( 4, : ) 0.8 ] );  grid on;
%          legend( linear_quadrupole.label, 'Location', 'North', 'Interpreter', 'Latex' );
%     xlabel( 'Angle [Degree]' );  ylabel( 'Pressure at 1 m [dB]' );
%     axis( [ -10 370  0 90 ] );



%% Problem 2c - Finite Line Source - 4 Sources Uniformly Separation - Smaller than a Wavelength

% See Lecture 21 - Distributed Sources (D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\21 Wednesday, April 2, 2025)


% d = wavelength / 8;  % Source separation distance.
% 
% xyz_sources = [ ... 
%     -4*d, 0, 0; ...
%     -3*d, 0, 0; ...
%     -2*d, 0, 0; ...
%     -1*d, 0, 0; ...
%     0, 0, 0; ...
%     1*d, 0, 0; ...
%     2*d, 0, 0; ...
%     3*d, 0, 0; ...
%     4*d, 0, 0; ...
%     ];  % m
% %
% Q_sources = [ ...
%     +1.4 + 1i*0; ...
%     +1.4 + 1i*0; ...
%     +1.4 + 1i*0; ...
%     +1.4 + 1i*0; ...
%     +1.4 + 1i*0; ...
%     +1.4 + 1i*0; ...
%     +1.4 + 1i*0; ...
%     +1.4 + 1i*0; ...
%     +1.4 + 1i*0; ...
%     ];  % m^3/s
% 
% 
% xyz_receivers = [ cos( theta )  sin( theta )  zeros( size( cos( theta ) ) ) ];
% 
% array_length = ( size( xyz_sources, 1 ) - 1 )*d;
% 
% 
% figure( 'Name', 'Finite Line Source - 4 Sources Uniformly Spaced' ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, 1e2, rho0, c );
%         L = 10*log10( abs(p).^2 );
%             polarplot( theta, L, 'Color', [ color_map( 1, : ) 0.8 ] ); hold on;
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, 2e3, rho0, c );
%         L = 10*log10( abs(p).^2 );
%             polarplot( theta, L, 'Color', [ color_map( 2, : ) 0.8 ] );
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, 5e3, rho0, c );
%         L = 10*log10( abs(p).^2 );
%             polarplot( theta, L, 'Color', [ color_map( 3, : ) 0.8 ] ); grid on;
% 
%     legend( sprintf( 'kL = %3.1f', (2*pi*1e2)/c*array_length ), sprintf( 'kL = %3.1f', (2*pi*2e3)/c*array_length ), sprintf( 'kL = %3.1f', (2*pi*5e3)/c*array_length ), 'Location', 'EastOutside' );
%     rlim( [ 0  1e2 ] );



%% Problem 2c - Finite Line Source - 4 Sources Uniformly Separation - Greater than a Wavelength

% See Lecture 21 - Distributed Sources (D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\21 Wednesday, April 2, 2025)


% d = wavelength;
% 
% xyz_sources = [ ... 
%     -3*d, 0, 0; ...
%     -2*d, 0, 0; ...
%     -1*d, 0, 0; ...
%     0, 0, 0; ...
%     1*d, 0, 0; ...
%     2*d, 0, 0; ...
%     3*d, 0, 0; ...
%     ];  % m
% %
% Q_sources = [ ...
%     +1 + 1i*0; ...
%     +1 + 1i*0; ...
%     +1 + 1i*0; ...
%     +1 + 1i*0; ...
%     +1 + 1i*0; ...
%     +1 + 1i*0; ...
%     +1 + 1i*0; ...
%     ];  % m^3/s
% 
% 
% xyz_receivers = [ cos( theta )  sin( theta )  zeros( size( cos( theta ) ) ) ];
% 
% array_length = ( size( xyz_sources, 1 ) - 1 )*d;
% 
% 
% figure( 'Name', 'Finite Line Source - 4 Sources Uniformly Spaced' ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, 1e2, rho0, c );
%         L = 10*log10( abs(p).^2 );
%             polarplot( theta, L, 'Color', [ color_map( 1, : ) 0.8 ] ); hold on;
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, 2e3, rho0, c );
%         L = 10*log10( abs(p).^2 );
%             polarplot( theta, L, 'Color', [ color_map( 2, : ) 0.8 ] );
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, 5e3, rho0, c );
%         L = 10*log10( abs(p).^2 );
%             polarplot( theta, L, 'Color', [ color_map( 3, : ) 0.8 ] ); grid on;
% 
%     legend( sprintf( 'kL = %3.1f', (2*pi*1e2)/c*array_length ), sprintf( 'kL = %3.1f', (2*pi*2e3)/c*array_length ), sprintf( 'kL = %3.1f', (2*pi*5e3)/c*array_length ), 'Location', 'EastOutside' );
%     rlim( [ 0  1e2 ] );


    
%% Problem 2d - Baffled Piston

% See Lecture 21 - Distributed Sources (D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\21 Wednesday, April 2, 2025)

% See:  https://demonstrations.wolfram.com/SunflowerSeedArrangements/


a = 5*wavelength;


N = 1e2;
    n = 0:1:( N - 1 );

phi = ( 1 + sqrt(5) ) / 2;

theta_n = (2*pi)/(phi^2) .* n;
r_n = a.*( sqrt( n./N ) );

x_n = r_n.*cos(theta_n);  y_n = r_n.*sin(theta_n);
    xyz_sources = [ x_n(:)  y_n(:)  zeros( N, 1 ) ];
    %
    % figure( 'Name', 'Distribution of Sources on Baffled Piston' ); ...
    %     scatter( x_n, y_n, 40.*ones( size(x_n) ), 'Marker', '*' );  grid on;
    %     xlabel( 'X-axis [m]' );  ylabel( 'Y-axis [m]' );
    %     daspect( [ 1 1 1 ] );
%
Q_sources = ones( size( xyz_sources, 1 ), 1 );

z_start = a/10;
    z_receivers = linspace( 0.1, 20.25*wavelength, 1e2 );
        z_receivers( 60:end ) = [ ];
            xyz_receivers = [ zeros( numel( z_receivers ), 1 )  zeros( numel( z_receivers ), 1 )  z_receivers(:) ];


m = 1:1:( a / wavelength );
    r_null = ( (a/wavelength)^2 - m.^2) ./ (2*m/wavelength);

figure( 'Name', 'Baffled Piston - ka = 8' ); ...

    p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, frequency, rho0, c );
        L = 10*log10( abs(p).^2 );

    h1 = plot( z_receivers, L );  hold on;
    h2 = line( [ r_null(1) r_null(1) ], [ 25 75 ], 'Color', [ 1, 0, 0, 0.5 ] );
    line( [ r_null(2) r_null(2) ], [ 25 75 ], 'Color', [ 1, 0, 0, 0.5 ] );
    line( [ r_null(3) r_null(3) ], [ 25 75 ], 'Color', [ 1, 0, 0, 0.5 ] );
    line( [ r_null(4) r_null(4) ], [ 25 75 ], 'Color', [ 1, 0, 0, 0.5 ] );  grid on;
        legend( [ h1 h2 ], { 'Sound Pressure on Z-axis', 'Theoretical Null Locations' }, 'Location', 'South' );
    axis( [ 0 4.5  15 78 ] );
    xlabel( 'Distance above the XY-plane [m]' );  ylabel( 'Sound Pressure [dB]' );


    
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

% https://www.mathworks.com/matlabcentral/answers/1771015-how-to-customize-a-colormap

% https://www.mathworks.com/help/matlab/ref/colororder.html#mw_a3ac0dbb-c8d8-48c9-977c-1bfa84840b08

% https://www.mathworks.com/matlabcentral/answers/493751-increase-levels-on-colorbar

% https://www.mathworks.com/matlabcentral/answers/216283-get-current-axes-from-multiple-figures

% https://www.mathworks.com/matlabcentral/answers/288557-how-to-join-every-row-of-an-array-of-strings-in-a-single-string-for-time-plotting-purposes

% https://www.mathworks.com/matlabcentral/answers/352290-set-colorbar-ticklabels-and-tickmarks

% https://www.mathworks.com/matlabcentral/answers/152426-sprintf-d-x-prints-out-exponential-notation-instead-of-decimal-notation


