


%% Synopsis

% Problem 1 - Simulation of Sources from Monopoles



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


OFFSET = 3;
            cmap = slanCM( 'Purples', numel( r ) + OFFSET );
                cmap(1:OFFSET, :) = [ ];
                    cmap = flipud( cmap );



%% Problem 2a - Pressure Versus Distance for a Monopole Source

% Monopole sound pressure verus distance has a decay of 6 dB per doubling of distance.

xyz_sources = [ ...
    0, 0, 0; ...
    ];
%
Q_sources = 1;  % m^3/s


x = 0.01:0.01:10;
    y = zeros( length( x ), 1 );  z = zeros( length( x ), 1 );
        xyz_receivers = [ x(:) y(:) z(:) ];


% figure( 'Name', 'Monople Source - Pressure Magnitude Versus Distance' ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, frequency, rho0, c );
%         L = 10*log10( abs(p).^2 );
%                 plot( x, L, 'Color', 'b' );  grid on;
%     legend( '1 kHz', 'Location', 'NorthEast' );
%     xlabel( 'Distance [m]' );  ylabel( 'Pressure Magnitude [$dB$]' );
%     set( gca, 'XScale', 'log' );


% Note(s):

%   With a higher frequency, the wave number, "k" will be larger.  Pressure will be higher at the same distance.
%   For a doubling in frequency, the pressure will be 20 dB larger.



%% Problem 2b - Directivity Pattern for a Monopole Source

xyz_sources = [ ... 
    0, 0, 0; ...
    ];  % m
%
Q_sources = 1.4;  % m^3/s


xyz_receivers = [ cos( theta )  sin( theta )  zeros( size( cos( theta ) ) ) ];

h = figure( 'Name', 'Monople Source - Directivity Pattern' ); ...

    p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, frequency, rho0, c );
        monopole.L = 10*log10( abs(p).^2 );

    polarplot( theta, monopole.L, 'Color', cmap( 1, : ) );
        tick_marks = monopole.L( 1 );


    colormap( gca( h ), flipud( cmap ) );
    % 
    h_colorbar = colorbar( );
        h_colorbar.Ticks = linspace( 0, 1, 4 );
        h_colorbar.TickLabels = num2cell( round( fliplr( tick_marks ) ) );

    temp2 = [ round( tick_marks, 1 )  1 ];
        temp3 = string( temp2 );

    monopole.label = sprintf( '%1g m, %3.1f dB', string( [ 1  round( tick_marks, 1 ) ] ) );

    h_colorbar.TickLabels = flipud( monopole.label );
    h_colorbar.Label.String = 'Sound Pressure [dB]';

    rlim( [ 0  80 ] );


figure( 'Name', 'Monople Source - Pressure Versus Angle' ); ...
    
    plot( theta.*180./pi, monopole.L, 'Color', cmap( 1, : ) );  grid on;
    legend( '1 kHz Monopole, Q = 1.4 $\frac{m^3}{s}$', 'Location', 'East', 'Interpreter', 'Latex' );
    xlabel( 'Angle [Degree]' );  ylabel( 'Pressure at 1 m [dB]' );
    axis( [ -10 370  0 80 ] );

return

%% Problem 2b - Dipole Directivity Pattern

xyz_sources = [ ... 
    0, 0, 0; ...
    0.1715, 0, 0; ...
    ];  % m
%
Q_sources = [ ...
    1 + 1i; ...
    -1 - 1i; ...
    ];  % m^3/s


% for frequency_index = 2  % 1 kHz
% 
%     h = figure( 'Name', 'Dipole - Spherical Spreading' ); ...
% 
%     OFFSET = 3;
%     cmap = slanCM( color_maps{ frequency_index }, numel( r ) + OFFSET );
%     cmap(1:OFFSET, :) = [ ];
%     cmap = flipud( cmap );
% 
%     tick_marks = [ ];
% 
%     for radius_index = 4
% 
%         p = sum_of_monopoles( xyz_sources, Q_sources, squeeze( xyz_receivers( radius_index, :, : ) ), f( frequency_index ), rho0, c );
%         L = 10*log10( abs(p).^2 );
% 
%         % polarplot( theta, r( radius_index ).*ones( size( theta ) ), 'Color', cmap( radius_index, : ) ); hold on;
%         polarplot( theta, L ); hold on;
%         tick_marks = [ L(1) tick_marks ];
% 
%     end
% 
% 
%     % colormap( gca( h ), flipud( cmap ) );
%     %
%     % h_colorbar = colorbar( );
%     %     h_colorbar.Ticks = linspace( 0, 1, numel( r ) );
%     %     h_colorbar.TickLabels = num2cell( round( fliplr( tick_marks ) ) );
%     %
%     % temp2 = [ round( fliplr( tick_marks ) ).'  r.' ];
%     %     temp3 = string( temp2 );
%     %
%     % for i = 1:1:size( temp3, 1 )
%     %     new_labels( i, : ) = sprintf( '%1g m, %i dB', fliplr( temp2( i, : ) ) );
%     % end
%     %
%     % h_colorbar.TickLabels = flipud( new_labels );
%     % h_colorbar.Label.String = 'Sound Pressure [dB]';
% 
%     % rlim( [ 0, 1.2*max( r ) ] );
%     % %
%     % rticks( r );
%     %     rticklabels( { 'r = 0.5', 'r = 2.4', 'r = 4.8', 'r = 9.6' } );
%     %         set( gca, 'RTickLabelRotation', 45 );
% 
% end


% h = [ ];
% 
% figure( 'Name', 'Dipole - Pressure Versus Angle' ); ...
%     plot( nan, nan );  hold on;
%     %
%     for index = 1:1:numel( f )
%         h_plot = plot( theta.*180./pi, L( index, : ), 'Color', color_set( index, : ) );  grid on;
%                     h = [ h h_plot ];
%     end
%     %
%     % legend( h, f_set, 'Location', 'EastOutside' );
%     xlabel( 'Angle [Degree]' );  ylabel( 'Pressure at 9.6 m [$dB$]' );
%     % axis( [ -10 370  0 60 ] );

% return

%% Problem 2b - Lateral Quadrupole Directivity Pattern

d = 1e-1;

xyz_sources = [ ... 
    +d, +d, 0; ...
    -d, +d, 0; ...
    -d, -d, 0; ...
    +d, -d, 0; ...
    ];  % m
%
Q_sources = [ ...
    +1; ...
    -1; ...
    +1; ...
    -1; ...
    ];  % m^3/s


% for frequency_index = 2  % 1 kHz
% 
%         h = figure( 'Name', 'Lateral Quadrupole Directivity Pattern' ); ...
% 
%         OFFSET = 3;
%         cmap = slanCM( color_maps{ frequency_index }, numel( r ) + OFFSET );
%         cmap(1:OFFSET, :) = [ ];
%         cmap = flipud( cmap );
% 
%         tick_marks = [ ];
% 
%         for radius_index = 4
% 
%             p = sum_of_monopoles( xyz_sources, Q_sources, squeeze( xyz_receivers( radius_index, :, : ) ), f( frequency_index ), rho0, c );
%             L = 10*log10( abs(p).^2 );
% 
%             polarplot( theta, L ); hold on;
%             tick_marks = [ L(1) tick_marks ];
% 
%         end
% 
% end

% return

%% Problem 2b - Linear Quadrupole Directivity Pattern

d = 1e-1;

xyz_sources = [ ... 
    +d, 0, 0; ...
    0, 0, 0; ...
    -d, 0, 0; ...
    ];  % m
%
Q_sources = [ ...
    -1; ...
    +2; ...
    -1; ...
    ];  % m^3/s


% for frequency_index = 2  % 1 kHz
% 
%     h = figure( 'Name', 'Linear Quadrupole Directivity Pattern' ); ...
% 
%     OFFSET = 3;
%     cmap = slanCM( color_maps{ frequency_index }, numel( r ) + OFFSET );
%     cmap(1:OFFSET, :) = [ ];
%     cmap = flipud( cmap );
% 
%     tick_marks = [ ];
% 
%     for radius_index = 4
% 
%         p = sum_of_monopoles( xyz_sources, Q_sources, squeeze( xyz_receivers( radius_index, :, : ) ), f( frequency_index ), rho0, c );
%         L = 10*log10( abs(p).^2 );
% 
%         polarplot( theta, L ); hold on;
%         tick_marks = [ L(1) tick_marks ];
% 
%     end
% 
% end



%% Problem 2c - Finite Line Source - 4 Sources Uniformly Separation - Smaller than a Wavelength

% See Lecture 21 - Distributed Sources (D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\21 Wednesday, April 2, 2025)


% d = 0.343 / ( 8 - 1 );  % Source separation distance.
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
%     +1 + 1i*0; ...
%     +1 + 1i*0; ...
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
% array_length = ( size( xyz_sources, 1 ) - 1 )*d;
% 
% h = figure( 'Name', sprintf( 'Finite Line Source - 4 Sources Uniformly Spaced - kL = %2.1f, L = %2.1f', k*array_length, array_length ) ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, squeeze( xyz_receivers( 4, :, : ) ), f( 2 ), rho0, c );
%         L = 10*log10( abs(p).^2 );
% 
%     polarplot( theta, L ); hold on;



%% Problem 2c - Finite Line Source - 4 Sources Uniformly Separation - Greater than a Wavelength

% See Lecture 21 - Distributed Sources (D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\21 Wednesday, April 2, 2025)


% d = 0.343;  % Source separation distance.
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
% array_length = ( size( xyz_sources, 1 ) - 1 )*d;
% 
% h = figure( 'Name', sprintf( 'Finite Line Source - 4 Sources Uniformly Spaced - kL = %2.1f, L = %2.1f', k*array_length, array_length ) ); ...
% 
%     p = sum_of_monopoles( xyz_sources, Q_sources, squeeze( xyz_receivers( 4, :, : ) ), f( 2 ), rho0, c );
%         L = 10*log10( abs(p).^2 );
% 
%     polarplot( theta, L ); hold on;


    
%% Problem 2d - Baffled Piston

% See Lecture 21 - Distributed Sources (D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\21 Wednesday, April 2, 2025)


d = 0.343;  % Source separation distance.

xyz_sources = [ ... 
    -3*d, 0, 0; ...
    -2*d, 0, 0; ...
    -1*d, 0, 0; ...
    0, 0, 0; ...
    1*d, 0, 0; ...
    2*d, 0, 0; ...
    3*d, 0, 0; ...
    ];  % m
%
Q_sources = [ ...
    +1 + 1i*0; ...
    +1 + 1i*0; ...
    +1 + 1i*0; ...
    +1 + 1i*0; ...
    +1 + 1i*0; ...
    +1 + 1i*0; ...
    +1 + 1i*0; ...
    ];  % m^3/s


array_length = ( size( xyz_sources, 1 ) - 1 )*d;

h = figure( 'Name', sprintf( 'Finite Line Source - 4 Sources Uniformly Spaced - kL = %2.1f, L = %2.1f', k*array_length, array_length ) ); ...

    p = sum_of_monopoles( xyz_sources, Q_sources, squeeze( xyz_receivers( 4, :, : ) ), f( 2 ), rho0, c );
        L = 10*log10( abs(p).^2 );

    polarplot( theta, L ); hold on;


    
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


































