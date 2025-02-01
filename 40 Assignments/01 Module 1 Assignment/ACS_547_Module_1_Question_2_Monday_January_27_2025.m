


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

rho0 = 1.21;  % Ratio of specific heats (unitless).
c = 343;  % Speed of sound in air (meters per second).

frequency_set = 0:1:5e3;  % Hertz



%% Dimensions.

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


    % Right-to-left - 1: duct;  2: muffler;  3: duct.
    %
    T_outlet = duct_segment_transfer_matrix( f, rho0, c, segment_lengths( 1 ), segment_areas ( 1 ) );
        T_outlet_muffler_connection = [ 1  0;  0  segment_areas( 1 ) / segment_areas( 2 ) ];
    %
    T_muffler = duct_segment_transfer_matrix( f, rho0, c, segment_lengths( 2 ), segment_areas( 2 ) );
        T_muffler_inlet_connection = duct_connection_transfer_matrix( segment_areas( 2 ), segment_areas( 3 ) );
    %    
    T_inlet = duct_segment_transfer_matrix( f, rho0, c, segment_lengths( 3 ), segment_areas( 3 ) );


    T_net = T_inlet * T_muffler_inlet_connection * T_muffler * T_outlet_muffler_connection * T_outlet * T_total;
        T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);

    Z = open_end_impedance( f, rho0, c, segment_lengths( 1 ), segment_areas( 1 ), outlet_flanged );
        TL( frequency_index ) = 10 * log10( 0.25 * abs( ( T11  +  segment_areas(3)*T12/(rho0*c)  +  (rho0*c)*T21/segment_areas(1)  +  T22 ) / 2 )^2 );

end

TL_parta = TL;



%% Part b - Double-tuned Expansion Chamber

annulus_area_squared_meters = pi/4 * ( segment_diameters(2)^2 - segment_diameters(1)^2 );

epsilon = segment_diameters(3) / segment_diameters(2);  % 0.2
%
% Using Ji (2005):
L_o = segment_diameters(3) * ( 0.8216 - 0.0644*epsilon - 0.694*epsilon^2 );


nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );
    
    T_total = [ 1 0; 0 1 ];

    
    % Right-to-left - 1: duct;  2: branch;  3: muffler;  4: branch;  5: duct.
    %
    T_outlet = duct_segment_transfer_matrix( f, rho0, c, segment_lengths( 1 ), segment_areas ( 1 ) );
    %
    k = 2*pi*f/c;
        Z_A = -1j*rho0*c/annulus_area_squared_meters*cot(k * ( dimensions.overhang + L_o ) );
            T_branch_1 = [ 1  0;  1/Z_A  1 ];
    %
    T_muffler = duct_segment_transfer_matrix( f, rho0, c, segment_lengths(2) - 2*segment_lengths(4), segment_areas( 2 ) );
    %
    T_branch_2 = [ 1  0;  1/Z_A  1 ];
    %    
    T_inlet = duct_segment_transfer_matrix( f, rho0, c, segment_areas( 3 ), segment_areas( 3 ) );


    T_net = T_inlet * T_branch_2 * T_muffler * T_branch_1 * T_outlet * T_total;
        T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);

    Z = open_end_impedance( f, rho0, c, segment_lengths( 1 ), segment_areas( 1 ), outlet_flanged );
        TL( frequency_index ) = 10 * log10( 0.25 * abs( ( T11  +  segment_areas(3)*T12/(rho0*c)  +  (rho0*c)*T21/segment_areas(1)  +  T22 ) / 2 )^2 );

end

TL_partb = TL;    



%% Part c - Cascaded, Double-tuned Expansion Chamber

nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );
    
    T_total = [ 1 0; 0 1 ];

    
    % Right-to-left - 1: duct;  2: branch;  3: muffler;  4: branch;  5;  branch;  6: muffler;  7: branch;  8: duct.
    %
    T_outlet = duct_segment_transfer_matrix( f, rho0, c, segment_lengths( 1 ), segment_areas ( 1 ) );
    %
    k = 2*pi*f/c;
        Z_A = -1j*rho0*c/annulus_area_squared_meters*cot(k * ( dimensions.overhang + L_o ) );
            T_branch_1 = [ 1  0;  1/Z_A  1 ];
    %
    T_muffler_1 = duct_segment_transfer_matrix( f, rho0, c, segment_lengths(4), segment_areas( 2 ) );
    %
    T_branch_2 = [ 1  0;  1/Z_A  1 ];
    %
    T_branch_3 = [ 1  0;  1/Z_A  1 ];
    %
    T_muffler_2 = duct_segment_transfer_matrix( f, rho0, c, segment_lengths(4), segment_areas( 2 ) );
    %
    T_branch_4 = [ 1  0;  1/Z_A  1 ];
    %
    T_inlet = duct_segment_transfer_matrix( f, rho0, c, segment_areas( 3 ), segment_areas( 3 ) );


    T_net = T_inlet * T_branch_4 * T_muffler_2 * T_branch_3 * T_branch_2 * T_muffler_1 * T_branch_1 * T_outlet * T_total;
        T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);

    Z = open_end_impedance( f, rho0, c, segment_lengths( 1 ), segment_areas( 1 ), outlet_flanged );
        TL( frequency_index ) = 10 * log10( 0.25 * abs( ( T11  +  segment_areas(3)*T12/(rho0*c)  +  (rho0*c)*T21/segment_areas(1)  +  T22 ) / 2 )^2 );

end

TL_partc = TL;



%% Part d - Cascaded, Double-tuned Expansion Chamber

nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );

    T_total = [ 1 0; 0 1 ];


    % Right-to-left - 1: duct;  2: branch;  3: muffler;  4: branch;  5;  branch;  6: muffler;  7: branch;  8: duct.
    %
    T_outlet = duct_segment_transfer_matrix( f, rho0, c, segment_lengths( 1 ), segment_areas ( 1 ) );
    %
    k = 2*pi*f/c; 
        Z_A = -1j*rho0*c/annulus_area_squared_meters*cot(k * ( 0.0508 + L_o ) );  % Changed to 2 inches.
            T_branch_1 = [ 1  0;  1/Z_A  1 ];
    %
    T_muffler_1 = duct_segment_transfer_matrix( f, rho0, c, 0.1016, segment_areas( 2 ) );  % Changed to 4 inches.
    %
    T_branch_2 = [ 1  0;  1/Z_A  1 ];
    %
    Z3 = -1j*rho0*c/annulus_area_squared_meters*cot(k * ( 0.1016 + L_o ) );  % Changed to 4 inches.
        T_branch_3 = [ 1  0;  1/Z3  1 ];
    %
    T_muffler_2 = duct_segment_transfer_matrix( f, rho0, c, 0.0508, segment_areas( 2 ) );  % Changed to 2 inches.
    %
    T_branch_4 = [ 1  0;  1/Z_A  1 ];
    %
    T_inlet = duct_segment_transfer_matrix( f, rho0, c, segment_areas( 3 ), segment_areas( 3 ) );


    T_net = T_inlet * T_branch_4 * T_muffler_2 * T_branch_3 * T_branch_2 * T_muffler_1 * T_branch_1 * T_outlet * T_total;


    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);

    Z = open_end_impedance( f, rho0, c, segment_lengths( 1 ), segment_areas( 1 ), outlet_flanged );
        TL( frequency_index ) = 10 * log10( 0.25 * abs( ( T11  +  segment_areas(3)*T12/(rho0*c)  +  (rho0*c)*T21/segment_areas(1)  +  T22 ) / 2 )^2 );

end

TL_partd = TL;



%% Plot all Transmission Loss Profiles

Y_LIMITS = [ -20  280 ];

figure( ); ...
    plot( frequency_set, TL_parta, 'LineWidth', 0.6, 'Color', 'b' );  hold on;
    plot( frequency_set, TL_partb, 'LineWidth', 0.8, 'LineStyle', '--', 'Color', 'b' );
    plot( frequency_set, TL_partc, 'LineWidth', 1.0, 'LineStyle', '-.', 'Color', 'b' );
    plot( frequency_set, TL_partd, 'LineWidth', 0.6, 'LineStyle', '-.', 'Color', 'b' );  grid on;
        legend( ...
            'Simple Expansion Chamber', ...
            'Double-tuned Expansion Chamber', ...
            'Cascaded Double-tuned Expansion Chamber', ...
            'Modified Cascaded Double-tuned Expansion Chamber', ...
            'Location', 'SouthOutside' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Transmission Loss Profiles' );
    %
    Ax = gca;
        Ax.XAxis.TickLabelInterpreter = 'latex';
        Ax.YAxis.TickLabelInterpreter = 'latex';
    %
    axis( [ -50  5e3+50  Y_LIMITS ] );
    %
    if ( PRINT_FIGURES == 1 )
        exportgraphics( gcf, 'Figure TL All Profiles.pdf', 'Append', true );
    end



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



%% Delete

% % Plot transmission loss profile.
% zero_set = [ 0  3385 ];
% 
% x = repmat( zero_set.', 1, 2 ).';
% y = repmat( [ -60; 60 ], 1, size( x, 2 ) );
% 
% figure( ); ...
%     plot( frequency_set, TL );  grid on;
%     line( x, y, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 0.8 );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
%     title( 'Muffler Transmission Loss Profile' );
%     %
%     Ax = gca;
%         Ax.XAxis.TickLabelInterpreter = 'latex';
%         Ax.YAxis.TickLabelInterpreter = 'latex';
%     %
%     axis( [ -50  5e3+50  Y_LIMITS ] );
%     %
%     if ( PRINT_FIGURES == 1 )
%         exportgraphics( gcf, 'Figure Problem 2b.pdf', 'Append', true );
%     end





