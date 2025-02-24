


%% Synopsis

% Problem 4 - Panel Tramission Loss



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



%% Define Panel

panel.length = 80e-2;  % m
panel.E = 200e9;  % Pa
panel.density = 7800;  % kg/m^3
panel.v = 0.29;  % Poisson's Ratio (unitless)
panel.thickness = 1.2e-3;  % m
panel.eta = 0.001;  % Loss factor (unitless)

c = 343;  % m/s
rho0 = 1.21;  % kg/m^3



%% Measured Panel Data

octave_band_frequencies = [ 63  125  250  500  1000  2000  4000  8000 ].';  % Hz
TL = [ 9  14  21  27  32  37  43  42 ].';  % dB



%% Problem 4a - Infinite, Rigid Panel Model with Normal Incidence

D = ( panel.E * panel.thickness.^3 ) / ( 12 * ( 1 - panel.v^2 ) );  % 31.4

ms = panel.density * panel.thickness;  % 9.4 kg/m^2

wo = pi^2 / panel.length * sqrt( D / ms );  % 22.6 radians/s
    s = wo^2 * ms;  % 4,785.9 kg radians / m^2s^2

fo = wo / (2*pi);  % 4 Hz


% Define an Anonymous function for the rigid panel with normal incidence.
h_tau_infinite_rigid_panel = @( f, wo, ms, s, rho0, c, eta)  4 ./ ( ( (2*pi.*f*ms - s./(2*pi.*f)) ./ (rho0 * c) ).^2  +  ( (wo*ms*eta) ./ (rho0*c) + 2 ).^2 );



%% Problem 4b - Infinite, Flexible Panel Model with Random Incidence

% The panel has bending waves.


% Part (i.)

% The critical frequency.
 critical_frequency = c^2 / (2*pi) * sqrt( ms / D );  % 10.22 kHz

% Verify the critical frequency using a 90 degree angle of incidence.
critical_frequency_verify_1 =  1./(2*pi) .* sqrt( ms / D ) .* ( c / sind( 90 )).^2;  % 10.22 kHz
        
% Verify the critical frequency using the properties of the panel.
critical_frequency_verify_2 = c^2 / ( 1.8 * panel.thickness * sqrt( panel.E / ( panel.density * ( 1 - panel.v^2) ) ) );  % 10.3 kHz


% Coincidence frequency for a 75 degree angle of incidence.
phi = 75;
    h_coincidence_frequency =  @( ms, D, c, phi )  1./(2*pi) * sqrt( ms / D ) .* ( c ./ sind( phi )).^2;
        h_coincidence_frequency( ms, D, c, phi );  % 10,949 Hz


% Define an Anonymous function for the flexible panel with random incidence.        
h_tau_infinite_flexible_panel = @( f, rho0, c, phi, D, eta )  ( 2*rho0.*c*secd(phi)).^2 ./ ( (2*rho0.*c*secd(phi) + D*eta*(2*pi.*f./c).^4./(2*pi.*f)*sind(phi)^4).^2  +  ...
    (2*pi.*f*ms - D*(2*pi.*f./c).^4./(2*pi.*f)*sind(phi)^4).^2 );



% Part (ii.) - Transmission loss for a 75 degree angle of incidence.
f = 0.1:1:100e3;

figure( ); ...
    plot( f, -10*log10( h_tau_infinite_flexible_panel( f, rho0, c, phi, D, panel.eta ) ), 'LineStyle', '-', 'Marker', 'none', 'Color', [0.00, 0.45, 0.74] );  grid on;
    text( 12e3, 5, sprintf( 'Coincidence\nFrequency\nof 10,494 Hz' ) );
    xlabel( 'Frequency [Hz] ' );  ylabel( 'Transmission Loss [dB]' );
    set( gca, 'XScale', 'log' );
    axis( [ 1 200e3 -5 80 ] );
    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.3 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Q4 TL for 75 AOI', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex



% Part (iii.)
f = 0.1:1:100e3;

eta = panel.eta;  phi_set = 0:10:90;  t_set = [ ];

figure( ); ...
    hold on;
    for phi = phi_set
        plot( f, -10*log10( h_tau_infinite_flexible_panel( f, rho0, c, phi, D, panel.eta ) ), 'LineStyle', '-', 'Marker', 'none', 'Color', [0.00, 0.45, 0.74] );
            t_set = [ t_set;  h_tau_infinite_flexible_panel( f, rho0, c, phi, D, panel.eta ) ];
    end
    %
    grid on;  box on;
    xlabel( 'Frequency [Hz] ' );  ylabel( 'Transmission Loss [dB]' );
    set( gca, 'XScale', 'log' );
    axis( [ 1 200e3 -5 80 ] );
    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.3 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Q4iii TL for 75 AOI', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex


tau_d = nanmean( t_set .* sind( 2*phi_set ).', 1 );



%% Combined Transmission Loss Plot

figure( ); ...
    h1 = stem( octave_band_frequencies, TL, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k' );  hold on;
    h2 = plot( f, -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta ) ), 'LineStyle', '-', 'Color', 'k' );
    h3 = plot( f, -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta ) ./(200*panel.eta) * ( 4*panel.length / ( panel.length^2 * critical_frequency ) ) ), 'LineStyle', '--', 'Marker', 'none', 'Color', 'k' );
    h4 = plot( f, -10*log10( tau_d ), 'LineStyle', '-', 'Marker', 'none', 'Color', [0.00, 0.45, 0.74] );
    h5 = line( [ 16e3 16e3 ], [ 0 65 ], 'Color', 'r' );  grid on;  % 16 kHz Demarcation
    legend( ...
        [ h1, h2, h3, h4, h5 ], ...
        'Target Transmission Loss', ...
        'Infinite Rigid Panel', ...
        'Infinite Flexible Panel with Diffuse Incidence', ...
        'Finite Flexible Panel Model', ...
        '16 kHz Demarcation', ...
        'Location', 'SouthOutside' );
            
    xlabel( 'Frequency [Hz] ' );  ylabel( 'Transmission Loss [dB]' );
    title( 'Measured Panel Transmission Losses' );
    set( gca, 'XScale', 'log' );
    axis( [ 1 200e3 -5 80 ] );
    %
    % Textheight:  744 pt. and Textwidth:  493 pt. from LaTex document
    %
    % set( gcf, 'units', 'point', 'pos', [ 200 200    493*0.8 744*0.45 ] );
    %     pos = get( gcf, 'Position' );
    %         set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', 'PaperSize', [pos(3), pos(4)] );
    %             print(gcf, 'Q4d TL for 75 AOI', '-dpdf', '-r0' );
%
% https://tex.stackexchange.com/questions/179382/best-practices-for-using-matlab-images-in-latex



%% Plot Data and Model - Different Side Materials

% h_tau_infinite_rigid_panel_side_materials = @( f, wo, ms, s, rho0, c, eta, n )  (4*n) ./ ( ( (2*pi.*f*ms - s./(2*pi.*f)) ./ (rho0 * c) ).^2  +  ( (wo*ms*eta) ./ (rho0*c) + n + 1 ).^2 );

% f = 1e-2:1e-2:20e3;
% 
% phi = 15;
% eta = panel.eta;
% 
% 
% figure( ); ...
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel_side_materials( f, wo, ms, s, rho0, c, panel.eta, 1 ) ), 'LineStyle', '-' );  hold on;
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel_side_materials( f, wo, ms, s, rho0, c, panel.eta, 1/3600 ) ), 'LineStyle', '-' );
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel_side_materials( f, wo, ms, s, rho0, c, panel.eta, 3600 ) ), 'LineStyle', '-' );  grid on;
%         %
%         legend( ...
%             'Same Fluid', ...
%             'Water to Air', ...
%             'Air to Water', ...
%             'Location', 'North' );
%     %
%     xlabel( 'Frequency [$\frac{\omega}{\omega_o}$] ' );  ylabel( 'Transmission Loss [dB]' );
%     title( 'Measured Panel Transmission Losses' );
%     set( gca, 'XScale', 'log' );
%     % axis( [ 40 12e3 -5 45] );



%% Change is Stiffness

% figure( ); ...
%     stem( octave_band_frequencies ./ (wo / (2*pi) ), TL, 'LineWidth', 0.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r' );  hold on;
%     %
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta ) ), 'LineStyle', '-' );
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s*100, rho0, c, panel.eta ) ), 'LineStyle', '--' );
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s*1e-2, rho0, c, panel.eta ) ), 'LineStyle', '--' );
%         %
%         legend( ...
%             'Target TL Values', ...
%             'Infinite Rigid Panel with Normal Incidence Sound', ...
%             'Infinite Rigid Panel with Normal Incidence Sound (s * 100)', ...
%             'Infinite Rigid Panel with Normal Incidence Sound (s / 100)', ...
%             'Location', 'North' );
%     %
%     xlabel( 'Frequency [$\frac{\omega}{\omega_o}$] ' );  ylabel( 'Transmission Loss [dB]' );
%     title( 'Measured Panel Transmission Losses - Change in Stiffness' );
%     set( gca, 'XScale', 'log' );
%     % axis( [ 40 12e3 -5 45] );



%% Change in Mass

% figure( ); ...
%     stem( octave_band_frequencies ./ (wo / (2*pi) ), TL, 'LineWidth', 0.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r' );  hold on;
%     %
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms*100, s, rho0, c, panel.eta ) ), 'LineStyle', '-' );
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta ) ), 'LineStyle', '-' );
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms*1e-2, s, rho0, c, panel.eta ) ), 'LineStyle', '-' );
%         %
%         legend( ...
%             'Target TL Values', ...
%             'Infinite Rigid Panel with Normal Incidence Sound', ...
%             'Infinite Rigid Panel with Normal Incidence Sound (s * 100)', ...
%             'Infinite Rigid Panel with Normal Incidence Sound (s / 100)', ...
%             'Location', 'North' );
%     %
%     xlabel( 'Frequency [$\frac{\omega}{\omega_o}$] ' );  ylabel( 'Transmission Loss [dB]' );
%     title( 'Measured Panel Transmission Losses - Change in Mass' );
%     set( gca, 'XScale', 'log' );
%     % axis( [ 40 12e3 -5 45] );



%% Change in Loss Factor

% figure( ); ...
%     stem( octave_band_frequencies ./ (wo / (2*pi) ), TL, 'LineWidth', 0.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r' );  hold on;
%     %
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta ) ), 'LineStyle', '-' );
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta*1e2 ) ), 'LineStyle', ':' );
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta*1e-2 ) ), 'LineStyle', ':' );
%         %
%         legend( ...
%             'Target TL Values', ...
%             'Infinite Rigid Panel with Normal Incidence Sound', ...
%             'Infinite Rigid Panel with Normal Incidence Sound (eta * 100)', ...
%             'Infinite Rigid Panel with Normal Incidence Sound (eta / 100)', ...
%             'Location', 'North' );
%     %
%     xlabel( 'Frequency [$\frac{\omega}{\omega_o}$] ' );  ylabel( 'Transmission Loss [dB]' );
%     title( 'Measured Panel Transmission Losses - Change in Loss Factor' );
%     set( gca, 'XScale', 'log', 'YScale', 'log' );
%     % axis( [ 40 12e3 -5 45] );



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


