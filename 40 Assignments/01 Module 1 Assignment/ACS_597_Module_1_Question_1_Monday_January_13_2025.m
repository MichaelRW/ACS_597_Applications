


%% Synopsis

% Question 1 - Cut-on Frequencies in Ducts and Pipes



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../40 Assignments/00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Define Constants and Anonymous Functions

c_air = 343;  % The speed of sound in air (meters per second).
c_water = 1500;  % The speed of sound in water (meters per second).

gamma = 1.4;  % The ratio of specific heats [unitless].
R = 287;  % The gas constant [Joules per ( kilogram * Kelvin)].


h_f_cut_on_rectangular_duct = @( c, L )  0.5 .* c ./ L;
%
% c - The speed of sound.
% L - The largest cross-section dimension of the rectangular duct.


h_f_cut_on_circular_duct = @( c, d )  0.568 .* c ./ d;
%
% c - The speec of sound.
% L - The diameter of the circular duct.


h_speed_of_sound_in_air = @( gamma, R, temperature_Kelvin)  sqrt( gamma .* R .* temperature_Kelvin );



%% Problem 1a

% The cross-sectional dimensions for the rectangular duct are:  Lx = 12 cm and Ly = 20 cm.

% The largest dimension is Ly = 20 cm or 0.2 m.

% The cut-on frequency is,
h_f_cut_on_rectangular_duct( c_air, 0.2 );  % 857.5 Hz (shown in class 858 Hz)
    fprintf( 1, '\n Problem 1a:  The lowest cut-on frequency for the rectangular pipe with air is %3.1f Hz.\n', h_f_cut_on_rectangular_duct( c_air, 0.2 ) );



%% Problem 1b

% The cross-sectional dimensions for the rectangular duct are:  Lx = 12 cm and Ly = 20 cm.

% The cross-sectional area of the rectangular duct is 12 cm * 20 cm = 240 cm^2 or 0.024 m^2.
rectangular_duct_cross_sectional_area = 0.12 * 0.20;  % 0.024 m^2

% The diameter of the circulat pipe is,
circular_duct_diameter = sqrt( 0.024 / pi ) * 2;  % 0.17481 meters
%
% Check:
    % pi * ( circular_duct_diameter / 2 )^2  CHECKED


% The cut-on frequency for the circular duct is,
h_f_cut_on_circular_duct( c_air, circular_duct_diameter );  % 1,114.5 Hz
    fprintf( 1, '\n Problem 1b:  The lowest cut-on frequency for the circular pipe (of equal area) with air is %3.1f Hz.\n', h_f_cut_on_circular_duct( c_air, circular_duct_diameter ) );



%% Problem 1c

% The cut-on frequency for the circular duct with water is,
h_f_cut_on_circular_duct( c_water, circular_duct_diameter );  % 4,873.9 Hz
    fprintf( 1, '\n Problem 1c:  The lowest cut-on frequency for the circular pipe (of equal area) with water is %3.1f Hz.\n', h_f_cut_on_circular_duct( c_water, circular_duct_diameter ) );

% The cut-on frequency should be higher because it is proportional to the
% speed of sound in a given medium.



%% Problem 1d

fprintf( 1, '\n Problem 1d:  See the figure.\n' );

temperature_range_celsius = 0:0.1:500;  % Celsius
    temperature_range_kelvin = temperature_range_celsius + 273.15;  % Kelvin


FONT_SIZE = 14;

figure( ); ...
    plot( temperature_range_celsius, h_f_cut_on_circular_duct( h_speed_of_sound_in_air( gamma, R, temperature_range_kelvin ), 0.05 ) ./ 1e3 );  grid on;
        legend( 'Duct Diameter = 5.0 cm', 'Location', 'East', 'FontSize', FONT_SIZE, 'Interpreter', 'Latex' );
        set( gca, 'FontSize', FONT_SIZE );
    %
    xlabel( 'Temperature [Celsius]', 'FontSize', FONT_SIZE );
        % xl = get( gca, 'xlabel' );    pxl = get( xl, 'position' );  pxl( 2 ) = 1.1 * pxl( 2 );
        %     set( xl, 'position', pxl );
    %
    ylabel( 'Lowest Cut-on Frequency [kHz]', 'FontSize', FONT_SIZE );
        % yl = get( gca, 'ylabel' );  pyl = get( yl, 'position' );  pyl( 1 ) = 1.2 * pyl( 1 );
        %     set( yl, 'position', pyl );
    %
    caption = sprintf( 'Lowest Cut-on Frequency for a Circular Pipe with Air Flow Versus Air Temperature\n' );
        title( caption, 'FontSize', FONT_SIZE );
    %
    ylim( [ 3  7 ] );



%% Problem 1e

fprintf( 1, '\n Problem 1e:  See Section Problem 1e of the Matlab script for the answers.\n\n' );


% Question:  Are cut-on frequencies higher for a circular or rectangular duct for a given cross-sectional area?

% The lowest cut-on frequency is higher for a circular duct than for a
% rectangular duct for a given cross-sectional area.

% For the dimensions given in class, the rectangular duct is not square.
% This produces a larger dimension and thus a smaller, lowest cut-on
% frequency.

% If the rectangular duct is square dimensions on the order of the circular
% duct diameter with the same cross-sectional area, the the cut-on
% frequencies are approximately equal.


% Question:  What about in air versus water?

% The lowest cut-on frequency is larger with water than air.  This due to
% the fact that the cut-on frequency is proportional to the speed of sound
% and the speed of sound in water is greater than it is in air.


% Question:  What about cold versus hot air?

% For a circular pipe, the cut-on frequency is higher in warm air than cold
% air.



%% Clean-up

if ( ~isempty( findobj( 'Type', 'figure' ) ) )
    monitors = get( 0, 'MonitorPositions' );
        if ( size( monitors, 1 ) == 1 )
            autoArrangeFigures( 2, 2, 1 );
        elseif ( 1 < size( monitors, 1 ) )
            autoArrangeFigures( 2, 2, 1 );
        end
end

if ( PRINT_FIGURES == 1 )
    saveas( gcf, 'Cut-on Frequency Versus Temperature - Sunday, January 19, 2025.pdf' );
end

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );


