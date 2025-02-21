


%% Synopsis

% Slide 8 - Noise Reduction and Transmission Loss



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



%% Parameters

c = 343;  % m/s
rho0 = 1.21;  % kg

panel.length = 80e-2;  % meters
%
panel.E = 200e9;  % Pascals
panel.density = 7800;  % kg / m^3
panel.v = 0.29;  % Poisson's Ratio (unitless)
panel.thickness = 1.2e-3;  % m
panel.eta = 0.001;  % Loss factor (unitless)



%% Panel Data

octave_band_frequencies = [ 63  125  250  500  1000  2000  4000  8000 ].';  % Hz
TL = [ 9  14  21  27  32  37  43  42 ].';  % dB


% figure( ); ...
%     stem( octave_band_frequencies, TL, 'LineWidth', 0.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r' );  grid on;
%     xlabel( 'Frequency [Hz] ' );  ylabel( 'Transmission Loss [dB]' );
%     title( 'Measured Panel Transmission Losses' );
%     set( gca, 'XScale', 'log' );
%     axis( [ 40 12e3 -5 45] );



%% Problem 4a - Infinite, Rigid Panel Model with Normal Incidence

D = ( panel.E * panel.thickness.^3 ) / ( 12 * ( 1 - panel.v^2 ) );  % 31.4

% From lecture 9 on Wednesday, February 12, 2025, the equivalent bending
% moment of the panel is a half-wavelength.
%
wavelength = 2 * panel.length;  % 1.6 meters

ms = panel.density * panel.thickness;

wo = pi^2 / panel.length * sqrt( D / ms );
    s = wo^2 * ms;

h_tau_infinite_rigid_panel = @( f, wo, ms, s, rho0, c, eta)  4 ./ ( ( (2*pi.*f*ms - s./(2*pi.*f)) ./ (rho0 * c) ).^2  +  ( (wo*ms*eta) ./ (rho0*c) + 2 ).^2 );



%% Problem 4b - Infinite, Flexible Panel Model with Random Incidence

h_tau_term1 = @( rho0, c, phi )  ( 2*rho0*c*secd(phi) ).^2;  % Checked

h_tau_term2 = @( rho0, c, phi, D, eta, f )  ( 2*rho0*c*secd(phi) + D*eta*(2*pi*f/c).^4/(2*pi*f) * sind(phi).^4 ).^2;

h_tau_term3 = @( f, ms, D, phi )  ( 2*pi*f*ms  -  D*(2*pi*f/c).^4/(2*pi*f) * sind(phi).^4).^2;

% h_tau_infinite_flexible_panel = @( f, rho0, c, phi, D, eta )  (2*rho0*c*secd(phi)).^2 ./ ( (2*rho0*c*secd(phi) + D*eta*(2*pi*f/c).^4/(2*pi*f)*sind(phi)^4).^2  +  ...
%     (2*pi*f*ms - D*(2*pi*f/c).^4/(2*pi*f)*sind(phi)^4).^2 );

% return
 
%% Plot Data and Model

f = 1e-2:1e-2:20e3;

phi = 15;
eta = panel.eta;


% figure( ); ...
%     stem( octave_band_frequencies ./ (wo / (2*pi) ), TL, 'LineWidth', 0.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r' );  hold on;
%     %
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta ) ), 'LineStyle', '-' );
%     %
%     plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_term1(rho0, c, phi) ./ ( h_tau_term2( rho0, c, phi, D, eta, f ) + h_tau_term3( f, ms ,D, phi) ) ) );  grid on;
%     % plot( f, -10*log10( h_tau_infinite_flexible_panel( f, rho0, c, 75, D, panel.eta ) ) );  grid on;
%         %
%         legend( ...
%             'Target TL Values', ...
%             'Infinite Rigid Panel with Normal Incidence Sound', ...
%             'Infinite Flexible Panel with Random Incidence Sound', ...
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

figure( ); ...
    stem( octave_band_frequencies ./ (wo / (2*pi) ), TL, 'LineWidth', 0.5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r' );  hold on;
    %
    plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta ) ), 'LineStyle', '-' );
    plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta*1e2 ) ), 'LineStyle', ':' );
    plot( f ./ (wo / (2*pi) ), -10*log10( h_tau_infinite_rigid_panel( f, wo, ms, s, rho0, c, panel.eta*1e-2 ) ), 'LineStyle', ':' );
        %
        legend( ...
            'Target TL Values', ...
            'Infinite Rigid Panel with Normal Incidence Sound', ...
            'Infinite Rigid Panel with Normal Incidence Sound (eta * 100)', ...
            'Infinite Rigid Panel with Normal Incidence Sound (eta / 100)', ...
            'Location', 'North' );
    %
    xlabel( 'Frequency [$\frac{\omega}{\omega_o}$] ' );  ylabel( 'Transmission Loss [dB]' );
    title( 'Measured Panel Transmission Losses - Change in Loss Factor' );
    set( gca, 'XScale', 'log', 'YScale', 'log' );
    % axis( [ 40 12e3 -5 45] );



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


