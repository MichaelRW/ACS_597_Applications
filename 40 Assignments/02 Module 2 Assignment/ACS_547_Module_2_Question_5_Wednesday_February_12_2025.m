


%% Synopsis

% Slide 8 - Noise Reduction and Transmission Loss

% Volume of the enclosure is much bigger than the machine.  Diffuse sound field in the enclosure.



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



%% Define Machine

machine.area = 3;  % m^2
machine.absorption = 0.07;  % Sabine
machine.D = 2;  % Unitless

machine.distance = 10;  % m



%% Data

octave_band_frequencies = [ 250  500  1000  2000  4000 ].';  % Hz
Lw = [ 105  115  106  108  119 ].';  % dB re: 1 pW
    % [ octave_band_frequencies  Lw ]

% figure( ); ...
%     h1 = stem( octave_band_frequencies, Lw, 'Marker', '.', 'MarkerSize', 12, 'Color', 'r' );  hold on;
%     h2 = line( [ 2e2 5e3 ], [ 30 30 ] );  grid on;
%         legend( [ h1 h2 ], 'Current Sound Pressure Levels', 'Target Sound Pressure Level', 'Location', 'North' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Sound Pressure Level [dB re: 20e-6 Pa]' );
%     title( 'Sound Power Level Versus Octave Band Center Frequency' );
%     %
%     axis( [ 150 6e3  0 140 ] );
%     set( gca, 'XScale', 'log' );



%% Per Octave Band Insertion Loss

Lp_10_meters = Lw  +  10*log10( machine.D /( 4 * pi * machine.distance ) );  % dB re: 20e-6 Pa
    % [ octave_band_frequencies  Lp_10_meters ]
%
%  The value of R is infinite.  The machine is outside in open air.

octave_band_IL = Lp_10_meters - 30;
    % [ octave_band_frequencies  octave_band_IL ]



%% Anonymous Function for Insertion Loss

h_IL_large = @( Sw, alpha_w, Si, alpha_i, TL )  10*log10(  1  +  (Sw*alpha_w  +  Si*alpha_i)./(Sw + Si)*10^(TL/10)  );



%% Find Values of TL and Aborption that will Meet the Target Insertion Loss - Ground Reflecting

% Assumption(s):
%
%   1.)  The ground is a hard reflecting survice
%   2.)  The enclosure is a cube.
%   3.)  There is no noise transmission through the ground.

enclosure.length = 2;  % m
enclosure.width = 2;  % m
enclosure.height = 2;  % m
    enclosure

IL_estimates = h_IL_large( )



%% Find Values of TL and Aborption that will Meet the Target Insertion Loss - Ground with Cover

% Assumption(s):
%
%   1.)  The ground is covered with the absorption material.
%   2.)  The enclosure is a cube.



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


