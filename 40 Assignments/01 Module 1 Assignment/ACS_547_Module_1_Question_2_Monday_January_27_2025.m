


%% Synopsis

% Question 2 - Muffler Design Comparison



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.5 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Constants

rho0 = 1.21;  % Air density (kg per m^3).
c = 343;  % Speed of sound in air (meters per second).

frequency_set = 0:1:5e3;  % Hertz



%% Dimensions

convert.inches_to_meters = 0.0254;
convert.foot_to_meters = 0.3048;

dimensions.inlet_diameter_meters = 2 * convert.inches_to_meters;  % 0.0508 meters
dimensions.inlet_length_meters = 6 * convert.foot_to_meters;  % 1.82 meters

dimensions.muffler_diameter_meters = 10 * convert.inches_to_meters;  % 0.254 meters
dimensions.muffler_length_meters = 18 * convert.inches_to_meters;  % 0.4572 meters

dimensions.outlet_diameter_meters = 2 * convert.inches_to_meters;  % 0.0508 meters
dimensions.outlet_length_meters = 1 * convert.foot_to_meters;  % 0.3048 meters

outlet_flanged = false;

dimensions.overhang = 3 *convert.inches_to_meters;  % 0.0762 meters

segment_diameters = [ ...
    dimensions.outlet_diameter_meters, ...
    dimensions.muffler_diameter_meters, ...
    dimensions.inlet_diameter_meters, ...
    ].';
%
h_area_from_diameter = @( d )  pi .* d.^2 ./ 4;
%
segment_areas = h_area_from_diameter( segment_diameters );

segment_lengths = [ ...
    dimensions.outlet_length_meters, ...
    dimensions.muffler_length_meters, ...
    dimensions.inlet_length_meters, ...
    dimensions.overhang, ...
    ].';



%% Part a - Simple Expansion Chamber

nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );

    T_total = [ 1 0; 0 1 ];

    T1 = duct_segment_transfer_matrix( f, rho0, c, 0.3048, 0.0020268 );  % Duct - Outlet
    T2 = duct_segment_transfer_matrix( f, rho0, c, 0.4572, 0.050671 );  % Duct
    T3 = duct_segment_transfer_matrix( f, rho0, c, 1.8288, 0.0020268 );  % Duct - Inlet

    T_net = T3 * T2 * T1 * T_total;
    % T_net = T_inlet *  T_total;  % Zero transmission loss for a straight duct.

    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
        TL( frequency_index ) = 10 * log10( abs( ( T11  +   0.0020268*T12/(rho0*c)  +  (rho0*c)*T21/0.0020268  +  T22 ) / 2 )^2 );
        %
        % The transmission loss calculation does not require a load impedance.

end

TL_parta = TL;  % The maximum peak value should be about 22 (21.952) dB.
%
% Expected behaviour:
%
%   1.)  0 dB at 0 Hz.
%   2.)  The transmission loss of a straight duct section is zero;  energy out equals energy in.

max( TL_parta );  % 22 dB



%% Part b - Double-tuned Expansion Chamber

annulus_area_squared_meters = pi/4 * ( 0.254^2 - 0.0508^2 );

L_o = 0;  % Assume that the Lo extension is neglible.


nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );
    
    T_total = [ 1 0; 0 1 ];

    T1 = duct_segment_transfer_matrix( f, rho0, c, (0.3048 + 0.0762), 0.0020268 );
    T3 = duct_segment_transfer_matrix( f, rho0, c, (0.4572 - 2*0.0762), 0.050671 );
    T5 = duct_segment_transfer_matrix( f, rho0, c, (1.8288 + 0.0762), 0.0020268 );

    k = 2*pi*f/c;
        Z_A = -1j*rho0*c/annulus_area_squared_meters*cot(k * ( 0.0762 + L_o ) );
            T2 = [ 1  0;  1/Z_A  1 ];
                T4 = T2;

    T_net = T5 * T4 * T3 * T2 * T1 * T_total;


    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  0.0020268*T12/(rho0*c)  +  (rho0*c)*T21/0.0020268  +  T22 ) / 2 )^2 );
        %
        % The transmission loss calculation does not require a load impedance.

end

TL_partb = TL;
%
% Expected behaviour:
%
%   1.)  0 dB at 0 Hz.
%   2.)  0 dB at same locations as a simple expansion chamber.
%   3.)  Peaks at 1,125 Hz and 3,376 Hz;  

% Frequency at which the quarter-wavelength is 0.0762 meters.
% 343 / ( 4 * 0.0762 );  % 1,125 Hz.

% Also work at three-quarter-wavelength.
% 3 * 1125;  % 3,375 Hz



%% Part c - Cascaded, Double-tuned Expansion Chamber

annulus_area_squared_meters = pi/4 * ( 0.254^2 - 0.0508^2 );

L_o = 0;  % Assume that the Lo extension is neglible.


nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );
        k = 2*pi*f/c;
            Z_A = -1j*rho0*c/annulus_area_squared_meters*cot(k * ( 0.0762 + L_o ) );
    
    T_total = [ 1 0; 0 1 ];

    T1 = duct_segment_transfer_matrix( f, rho0, c, (0.3048 + 0.0762), 0.0020268 );  % Duct - Outlet
    T2 = [ 1  0;  1/Z_A  1 ];  % Straight Side Branch
    T3 = duct_segment_transfer_matrix( f, rho0, c, (0.2286 - 2*0.0762), 0.050671 );  % Duct
    T4 = T2;  % Straight Side Branch
    T5 = duct_segment_transfer_matrix( f, rho0, c, 2*0.0762, 0.050671 );  % Duct
    T6 = T2;  % Straight Side Branch
    T7 = T3;  % Duct
    T8 = T2;  % Straight Side Branch
    T9 = duct_segment_transfer_matrix( f, rho0, c, (1.8288 + 0.0762), 0.0020268 );  % Duct - Inlet

    T_net = T9 * T8 * T7 * T6 * T5 * T4 * T3 * T2 * T1 * T_total;

    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  0.0020268*T12/(rho0*c)  +  (rho0*c)*T21/0.0020268  +  T22 ) / 2 )^2 );

end

TL_partc = TL;



%% Part d - Cascaded, Double-tuned Expansion Chamber

nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );
        k = 2*pi*f/c;
            Z_A = -1j*rho0*c/annulus_area_squared_meters*cot(k * ( 0.0762 + L_o ) );
    
    T_total = [ 1 0; 0 1 ];

    T1 = duct_segment_transfer_matrix( f, rho0, c, (0.3048 + 0.0762), 0.0020268 );  % Duct - Outlet
    T2 = [ 1  0;  1/Z_A  1 ];  % Straight Side Branch
    T3 = duct_segment_transfer_matrix( f, rho0, c, 0.0076, 0.050671 );  % Duct
    T4 = T2;  % Straight Side Branch
    T5 = duct_segment_transfer_matrix( f, rho0, c, 2*0.0762, 0.050671 );  % Duct
    T6 = T2;  % Straight Side Branch
    T7 = duct_segment_transfer_matrix( f, rho0, c, (0.29718 - 2*0.0762), 0.050671 );  % Duct
    T8 = T2;  % Straight Side Branch
    T9 = duct_segment_transfer_matrix( f, rho0, c, (1.8288 + 0.0762), 0.0020268 );  % Duct - Inlet

    T_net = T9 * T8 * T7 * T6 * T5 * T4 * T3 * T2 * T1 * T_total;

    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  0.0020268*T12/(rho0*c)  +  (rho0*c)*T21/0.0020268  +  T22 ) / 2 )^2 );

end

TL_partd = TL;



%% Plot Transmission Loss Profiles

Y_LIMITS = [ -5  320 ];

h_figure_1 = figure( ); ...
    plot( frequency_set, TL_parta, 'LineWidth', 1.0, 'LineStyle', ':', 'Color', 'r' );  hold on;
    plot( frequency_set, TL_partb, 'LineWidth', 0.9, 'LineStyle', '--', 'Color', 'b' );
    plot( frequency_set, TL_partc, 'LineWidth', 0.9, 'LineStyle', '-', 'Color', 'k' );  grid on;
        legend( ...
            'Simple Expansion Chamber', ...
            'Double-tuned Expansion Chamber', ...
            'Cascaded Double-tuned Expansion Chamber', ...
            'Location', 'SouthOutside' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Transmission Loss [dB]' );
    title( 'Transmission Loss Profiles' );
    %
    Ax = gca;
        Ax.XAxis.TickLabelInterpreter = 'latex';
        Ax.YAxis.TickLabelInterpreter = 'latex';
    %
    axis( [ -50  5e3+50  Y_LIMITS ] );


Y_LIMITS = [ -5  315 ];

h_figure_2 = figure( ); ...
    plot( frequency_set, TL_partc, 'LineWidth', 1.0, 'LineStyle', '-', 'Color', 'k' );  hold on;
    plot( frequency_set, TL_partd, 'LineWidth', 0.6, 'LineStyle', '--', 'Color', 'b' );  grid on;
        legend( ...
            'Cascaded Double-tuned Expansion Chamber', ...
            'Modified Cascaded Double-tuned Expansion Chamber', ...
            'Location', 'SouthOutside' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Transmission Loss [dB]' );
    title( 'Transmission Loss Profiles - Cascaded and Modified Double-tuned Cascaded Systems' );
    %
    Ax = gca;
        Ax.XAxis.TickLabelInterpreter = 'latex';
        Ax.YAxis.TickLabelInterpreter = 'latex';
    %
    axis( [ -50  5e3+50  Y_LIMITS ] );



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
        exportgraphics( h_figure_1, 'Assignment 1 - Question 2 Figure All TL Profiles.pdf', 'Append', true );
        exportgraphics( h_figure_2, 'Assignment 1 - Question 2 Figure Comparison TL Plot For Cascaded Systems.pdf', 'Append', true );
end

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );


