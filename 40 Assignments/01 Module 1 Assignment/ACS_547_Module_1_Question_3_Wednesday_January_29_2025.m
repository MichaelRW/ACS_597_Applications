


%% Synopsis

% Bugle Recorder

% From ACS 597 - Noise Control Applications, Lecture 3 (Wednesday, January 22, 2025)

% Material from lectures 2 and 3.

% Assumptions:  source is high-impedence;  load is Z_open;  no mean flow;



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

rho0 = 1.21;  % Density of air (kg per cubic-meter).
c = 343;  % Speed of sound in air (meters per second).



%% Define Shape

L_mouth_piece = 0.09;  % Meters

pipe.inner_diameter = 0.009;  % Meters
pipe.thickness = 0.004;  % Meters

% The recorder is unflanged.

hole_diameter = 0.006;  % Meters



%% Part a

% Determine the length of the recorder to produce 523 Hz.

% The total length of the recorder, including the 0.09 meter long mouthpiece, is L.

a = 0.009 / 2;  % Meters
    L_o = 0.61*a;  % Slide 18 of Lecture 2 slide set.

f = 523;  % Hz
    k = 2*pi*f/c;  % The wave number for the respective frequency.    

S = pi/4*(0.009)^2;  % squared-meters


test_lengths = 0:0.0001:0.1;
    test_lengths = test_lengths + 0.09;


nLengths = length( test_lengths );
    A = zeros( nLengths, 1 );


for iLength = 1:1:nLengths

    L = test_lengths(iLength);

    T_total = [ 1 0; 0 1 ];

    L_e = L + L_o;
        Z = 1j * rho0 * c / S * tan( k* L_e );

    T = [ ...
    cos(k*L),                           1j*rho0*c/S*sin(k*L); ...
    1j*S/(rho0*c)*sin(k*L),      cos(k*L) ...
    ];


    T_total = T * T_total;
        T11 = T_total(1, 1);  T12 = T_total(1, 2);

    A( iLength ) = -10*log10( abs( T11 + T12 / Z )^2 );

end


figure( ); ...
    plot( test_lengths * 1e3, A );  grid on;
    xlabel( 'Total Recorder Length [mm]' );  ylabel( 'Amplitude [dB]' );
    title( 'Amplification Versus Recorder Length' );



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


