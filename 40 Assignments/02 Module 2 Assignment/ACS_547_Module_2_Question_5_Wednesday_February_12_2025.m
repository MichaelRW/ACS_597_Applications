 


%% Synopsis

% Problem 5 - Large Enclosure Design



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



%% Define Machine

machine.area = 3;  % m^2
machine.absorption = 0.07;  % m^2 or Sabin
machine.D = 1;  % Unitless - In air.

machine.distance = 10;  % m



%% Measurement Data

octave_band_frequencies = [ 250  500  1000  2000  4000 ].';  % Hz
Lw = [ 105  115  106  108  119 ].';  % dB re: 1 pW

% figure( ); ...
%     h1 = stem( octave_band_frequencies, Lw, 'Marker', '.', 'MarkerSize', 12, 'Color', 'r' );  hold on;
%     h2 = line( [ 2e2 5e3 ], [ 30 30 ] );  grid on;
%         legend( [ h1 h2 ], 'Current Sound Pressure Level', 'Target Sound Pressure Level', 'Location', 'North' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Sound Pressure Level [dB re: 20e-6 Pa]' );
%     title( 'Sound Power Level Versus Octave Band Center Frequency' );
%     %
%     axis( [ 150 6e3  0 140 ] );
%     set( gca, 'XScale', 'log' );



%% Per Octave Band Insertion Loss

Lp_10_meters = Lw  +  10*log10( machine.D /( 4 * pi * machine.distance^2 ) );  % dB re: 20e-6 Pa
%
%  The value of R is infinite.  The machine is outside in open air.

octave_band_IL = Lp_10_meters - 30;



%% Define Anonymous Function for Insertion Loss

h_IL_large = @( Sw, alpha_w, Si, alpha_i, TL )  10*log10(  1  +  (Sw*alpha_w  +  Si*alpha_i)./(Sw + Si)*10^(TL/10)  );



%% Find Values of TL and Aborption that will Meet the Target Insertion Loss - Ground Reflecting

% Assumption(s):
%
%   1.)  The enclosure is a cube.
%   2.)  The machine sits on the ground.
%   2.)  There is no noise transmission through the ground.

enclosure.dimension = 2;  % m
    enclosure.area = 5 * enclosure.dimension^2;  % 20 m^2



% https://www.cecoenviro.com/wp-content/uploads/2023/12/Acoustic-Enclosures-8pp-A4-web.pdf
% https://www.controlnoise.com/product/acoustic-enclosures/


% Volume of the enclosure is much bigger than the machine.  Diffuse sound field in the enclosure.


switch ( 4 )

    case 1

        % https://www.controlnoise.com/wp-content/uploads/2022/02/Acoustic-Enclosures-Datasheet.pdf

        % 250 Hz - QBV-2
        alpha_w = 0.27;  % From specification sheet.
        TL = 22;  % From specification sheet.
            IL_estimates(1) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 500 Hz - QBV-2
        alpha_w = 0.96;  % From specification sheet.
        TL = 28;  % From specification sheet.
            IL_estimates(2) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 1 kHz - QBV-2
        % alpha_w = 1.13;  % From specification sheet.
        alpha_w = 0.99;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 40;  % From specification sheet.
            IL_estimates(3) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 2 kHz - QBV-2
        % alpha_w = 1.08;  % From specification sheet.
        alpha_w = 0.99;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 56;  % From specification sheet.
            IL_estimates(4) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 4 kHz - QBV-2
        alpha_w = 0.99;  % From specification sheet.
        TL = 61;  % From specification sheet.
            IL_estimates(5) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );


    case 2

        % https://www.controlnoise.com/wp-content/uploads/2022/02/Acoustic-Enclosures-Datasheet.pdf

        % Note:  Abosrption values are carried over from QBV-2.

        % 250 Hz - QBV-3
        alpha_w = 0.96;  % From specification sheet.
        TL = 25;  % From specification sheet.
            IL_estimates(1) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 500 Hz - QBV-3
        alpha_w = 0.87;  % From specification sheet.
        TL = 33;  % From specification sheet.
            IL_estimates(2) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 1 kHz - QBV-3
        % alpha_w = 1.13;  % From specification sheet.
        alpha_w = 0.66;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 46;  % From specification sheet.
            IL_estimates(3) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 2 kHz - QBV-3
        % alpha_w = 1.08;  % From specification sheet.
        alpha_w = 0.47;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 53;  % From specification sheet.
            IL_estimates(4) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 4 kHz - QBV-3
        alpha_w = 0.30;  % From specification sheet.
        TL = 58;  % From specification sheet.
            IL_estimates(5) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );


    case 3

        % HTL-4 from specification sheet shown in class.

        % 250 Hz
        % alpha_w = 1.13;  % From specification sheet.
        alpha_w = 0.99;
        TL = 39;  % From specification sheet.
            IL_estimates(1) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 500 Hz
        % alpha_w = 1.12;  % From specification sheet.
        alpha_w = 0.99;
        TL = 59;  % From specification sheet.
            IL_estimates(2) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 1 kHz
        % alpha_w = 1.09;  % From specification sheet.
        alpha_w = 0.99;
        TL = 68;  % From specification sheet.
            IL_estimates(3) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 2 kHz
        % alpha_w = 1.03;  % From specification sheet.
        alpha_w = 0.99;
        TL = 67;  % From specification sheet.
            IL_estimates(4) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 4 kHz
        alpha_w = 0.91;  % From specification sheet.
        TL = 72;  % From specification sheet.
            IL_estimates(5) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );


    case 4

        % https://www.acousticsworld.com/machine-acoustic-enclosures/

        % 250 Hz
        % alpha_w = 1.13;  % From specification sheet.
        alpha_w = 0.99;
        TL = 34;  % From specification sheet.
            IL_estimates(1) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 500 Hz
        % alpha_w = 1.12;  % From specification sheet.
        alpha_w = 0.99;
        TL = 48;  % From specification sheet.
            IL_estimates(2) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 1 kHz
        % alpha_w = 1.09;  % From specification sheet.
        alpha_w = 0.99;
        TL = 61;  % From specification sheet.
            IL_estimates(3) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 2 kHz
        % alpha_w = 1.03;  % From specification sheet.
        alpha_w = 0.99;
        TL = 66;  % From specification sheet.
            IL_estimates(4) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 4 kHz
        alpha_w = 0.91;  % From specification sheet.
        TL = 70;  % From specification sheet.
            IL_estimates(5) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );

end


[ octave_band_IL    IL_estimates.'    (octave_band_IL - IL_estimates.') ]


% What is the most restrictive case?



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

    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.45 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Q4d TL for 75 AOI', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex


