 


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
set( 0, 'DefaultLineLineWidth', 1.2 );
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

% Volume of the enclosure is much bigger than the machine  Diffuse sound field in the enclosure.


switch ( 1 )

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



%% Calculate the Differences

[ octave_band_IL    IL_estimates.'    (octave_band_IL - IL_estimates.') ]

data = [ ...
    250      44    28.0    19.7   10.6    5.6; ...
    500      54    26.7    22.2    6.6    -4.4; ...
    1000    45    5.6    1.4    -15.4    -22.4; ...
    2000    47    -8.4    -2.2    -18.4    -19.4; ...
    4000    58    -2.4    5.7     -11.0    -13.0 ];

figure( ); ...
    h1 = plot( data(:, 1 ), data( :, 3 ), 'LineStyle', '-', 'Color', 'k', 'Marker', '+', 'MarkerSize', 8 );  hold on;
    h2 = plot( data(:, 1 ), data( :, 4 ), 'LineStyle', '--', 'Color', 'k', 'Marker', '+', 'MarkerSize', 8 );
    h3 = plot( data(:, 1 ), data( :, 5 ), 'LineStyle', ':', 'Color', 'k', 'Marker', '+', 'MarkerSize', 8 );
    h4 = plot( data(:, 1 ), data( :, 6 ), 'LineStyle', '-.', 'Color', 'k', 'Marker', '+', 'MarkerSize', 8 );  grid on;
    line( [200 5500], [0 0 ], 'Color', 'r' );
    legend( [ h1 h2 h3 h4 ], { 'Case 1 - Control Noise QBV 2', 'Case 2 - Control Noise QBV 3', 'Case 3 - Acoustics World HTL (100MM)', 'Case 4 - HTL 4 (class table)' }, 'Location', 'SouthOutside' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Insertion Loss Difference [dB]' );
    %
    axis( [ 150 6e3  -30 30 ] );
    set( gca, 'XScale', 'log' );
    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.5 ] );
        pos = get( gcf, 'Position' );
            set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Q5 IL Plot', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex



%% Clean-up

% if ( ~isempty( findobj( 'Type', 'figure' ) ) )
%     monitors = get( 0, 'MonitorPositions' );
%         if ( size( monitors, 1 ) == 1 )
%             autoArrangeFigures( 2, 2, 1 );
%         elseif ( 1 < size( monitors, 1 ) )
%             autoArrangeFigures( 2, 2, 1 );
%         end
% end


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


