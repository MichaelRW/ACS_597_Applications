


%% Synopsis

% Acoustic Horn Example

% From ACS 597 - Noise Control Applications, Lecture 3 (Wednesday, January 22, 2025)

% Material from lectures 2 and 3.



%% To Do

% Focus on interpretation;  "big picture" understanding.

% Do not repeat the prompt\question.

% Code should have "meat" of understanding;  avoid over-commenting.

% Code should use intuitive variable names.



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../40 Assignments/00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.5 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Constants

rho0 = 1.21;  % Density of air (kilograms per cubic-meter).
c = 343;  % Speed of sound in air (meters per second).



%% Define Shape

number_of_segments = 20;
    segment_lengths = ones( number_of_segments, 1 ) .* 0.1;  % 20 segments, each 0.1 m long (2 meters total).
    segment_diameters = ones( number_of_segments, 1 ) .* 0.01;  % 1 cm cross-section.



%% Calculation

frequency_set = 1:1:5e3;  % Hertz

horn_amplification = horn_amplification( frequency_set, segment_diameters, segment_lengths, rho0, c );



%% Calculate Peak Placment for Plot

f_base = c / 4;
    peak_set = f_base .* ( 1:1:58 );



%% Plot

x = repmat( peak_set.', 1, 2 ).';
y = repmat( [ -60; 60 ], 1, size( x, 2 ) );

figure( ); ...
    plot( frequency_set, horn_amplification );  hold on;
    line( x, y, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 0.8 );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Horn Amplification Profile' );
    %
    Ax = gca;
        Ax.XAxis.TickLabelInterpreter = 'latex';
        Ax.YAxis.TickLabelInterpreter = 'latex';
    %
    axis( [ -50  5e3+50  -65 65 ] );



%% Interpretation

% The total length of this horn is 2 meters.

% At a 2 meter wavelength, the corresponding frequency is 171.5 Hz.
f = c / 2;  % 171.5 Hz

% The peaks occur at frequencies of integer multiples of the half-wavelength.

% The broad shape (gradually increasing) comes from the end connection.
%
% ?






    


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



%% Reference(s)


