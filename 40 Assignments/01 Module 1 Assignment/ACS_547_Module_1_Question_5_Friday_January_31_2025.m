


%% Synopsis

% Question 5 - Intake Duct Silencer



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



%% Problem 5a

% The considered dip frequency 940 Hz.

% The additional attentuatio required is 14 dB.

liner_thickness = 0.0381;  % meters
    h = sqrt( ( pi*( 0.1016 - 2*liner_thickness  )^2 ) / 4 );  % 0.0113 meters

m = 1 + liner_thickness/h;  % 2.7

% As noted in the discussion, there are two cases to be considered here.
%
% Figure 8.37 considers the combined effect of the expansion and the liner.
% In Problem 4, the effect of the expansion was calculated, so including it
% here would include its effect twice.  Therefore, with an m of 1, there is
% not additional attenuation with the linear.
%
% The second case is to consider the m value of 2.7, which provides
% additional attenuation.  From Figure 8.37, the total attenuation of 
% the lining is 10 dB.
attenuation_rate = 10 / 0.127;  % 78.4  dB per meter



%% Problem 5b

% For a circular duct.

liner_thickness = 0.0381;  % meters
half_diameter_of_open_orifice = 0.0127;  % meters

h = sqrt(pi) * half_diameter_of_open_orifice / 2;  % 0.0113 meters



%% Problem 5c

% The liner thickness ratio is,
liner_thickness / h;  % 3.3852 unitless

% The normalized frequency is,
( 2 * h ) / ( 343 / 940 );  % 0.06169 unitless



%% Problem 5d

% The approximate resistivity parameter is 4.



%% Problem 5e

% Calculate the flow resistivity.

rho0 = 1.21;  % Density of air (kg per cubic-meter).
c = 343;  % Speed of sound in air (meters per second).

R1 = 16 * rho0*c / liner_thickness;  % 1.74e5 kg / m^3*s


%% Clean-up

if ( ~isempty( findobj( 'Type', 'figure' ) ) )
    monitors = get( 0, 'MonitorPositions' );
        if ( size( monitors, 1 ) == 1 )
            autoArrangeFigures( 2, 2, 1 );
        elseif ( 1 < size( monitors, 1 ) )
            autoArrangeFigures( 2, 2, 1 );
        end
end


fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


