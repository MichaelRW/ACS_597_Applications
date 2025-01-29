


%% Synopsis

% Acoustic Horn Example

% From ACS 597 - Noise Control Applications, Lecture 3 (Wednesday, January 22, 2025)

% Material from lectures 2 and 3.


% Assumptions:
%
%   No resistive terms (i.e., R_A is zero).
%   No flow.



%% To Do

% Focus on interpretation;  "big picture" understanding.

% Do not repeat the prompt\question.

% Code should have "meat" of understanding;  avoid over-commenting.

% Code should use intuitive variable names.



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

frequency_set = 0:0.1:5e3;  % Hertz



%% Measurements

% Conversions Factors
convert.inches_to_meters = 0.0254;
convert.foot_to_meters = 0.3048;


% Inlet
dimensions.inlet_diameter_meters = 2 * convert.inches_to_meters;  % 0.0508 meters
dimensions.inlet_length_meters = 6 * convert.foot_to_meters;  % 1.82 meters

% Muffler
dimensions.muffler_diameter_meters = 10 * convert.inches_to_meters;  % 0.254 meters
dimensions.muffler_length_meters = 18 * convert.inches_to_meters;  % 0.4572 meters

% Outlet
dimensions.outlet_diameter_meters = 2 * convert.inches_to_meters;  % 0.0508 meters
dimensions.outlet_length_meters = 1 * convert.foot_to_meters;  % 0.3048 meters
%
outlet_flanged = false;



%% Part a - Simple Expansion Chamber

h_area_from_diameter = @( d )  pi .* d.^2 ./ 4;

segment_diameters = [ ...
    dimensions.outlet_diameter_meters, ...
    dimensions.muffler_diameter_meters, ...
    dimensions.inlet_diameter_meters, ...
    ];

segment_areas = h_area_from_diameter( segment_diameters );

segment_lengths = [ ...
    dimensions.outlet_length_meters, ...
    dimensions.muffler_length_meters, ...
    dimensions.inlet_length_meters, ...
    ];



f = 100;  % Hz

T_total = [ 1 0; 0 1 ];  % Start with the identity matrix.
 
T_outlet = duct_segment_transfer_matrix( f, rho0, c, segment_lengths( 1 ), segment_areas ( 1 ) );

T_outlet_muffler_connection = [ 1  0;  0  segment_areas( 1 ) / segment_areas( 2 ) ];

T_muffler = duct_segment_transfer_matrix( f, rho0, c, segment_areas( 2 ), segment_areas( 2 ) );

T_muffler_inlet_connection = [ 1 0;  0 segment_areas( 2 ) / segment_areas( 3 ) ];
    duct_connection_transfer_matrix( segment_areas( 2 ), segment_areas( 3 ) );

T_inlet = duct_segment_transfer_matrix( f, rho0, c, segment_areas( 3 ), segment_areas( 3 ) );


T_net = T_inlet * T_muffler_inlet_connection * T_muffler * T_outlet_muffler_connection * T_outlet * T_total;


% Z = open_end_impedance( f, rho0, c, segment_lengths( 1 ), segment_areas( 1 ), outlet_flanged )
% 
% transmission_loss = 10 * log10( abs( ( T_net(1, 1)  +  segment_areas(1)*T_net(1, 2)/(rho0*c)  +  (rho0*c)*T_net(2, 1)/segment_areas(1)  +  T_net(2, 2) ) / 2 )^2 )



nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );


for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );
    
    T_total = [ 1 0; 0 1 ];  % Start with the identity matrix.

    
    T_outlet = duct_segment_transfer_matrix( f, rho0, c, segment_lengths( 1 ), segment_areas ( 1 ) );
        T_outlet_muffler_connection = [ 1  0;  0  segment_areas( 1 ) / segment_areas( 2 ) ];
    %
    T_muffler = duct_segment_transfer_matrix( f, rho0, c, segment_areas( 2 ), segment_areas( 2 ) );
        T_muffler_inlet_connection = duct_connection_transfer_matrix( segment_areas( 2 ), segment_areas( 3 ) );
    %    
    T_inlet = duct_segment_transfer_matrix( f, rho0, c, segment_areas( 3 ), segment_areas( 3 ) );
    
    T_net = T_inlet * T_muffler_inlet_connection * T_muffler * T_outlet_muffler_connection * T_outlet * T_total;


    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);

    Z = open_end_impedance( f, rho0, c, segment_lengths( 1 ), segment_areas( 1 ), outlet_flanged );
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  segment_areas(1)*T12/(rho0*c)  +  (rho0*c)*T21/segment_areas(1)  +  T22 ) / 2 )^2 );

end  % End:  for f = frequency_set



%% Plot Tranmission Loss Profile

close all;  clc;


peak_set = [ 0  3385 ]

x = repmat( peak_set.', 1, 2 ).';
y = repmat( [ -60; 60 ], 1, size( x, 2 ) );

figure( ); ...
    plot( frequency_set, TL );  grid on;
    line( x, y, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 0.8 );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Muffler Transmission Loss Profile' );
    %
    Ax = gca;
        Ax.XAxis.TickLabelInterpreter = 'latex';
        Ax.YAxis.TickLabelInterpreter = 'latex';
    %
    axis( [ -50  5e3+50  -10 60 ] );

return

%% Cylindrical Duct with Constant Diameter

horn_amplification = horn_amplification( frequency_set, segment_diameters, segment_lengths, rho0, c );

f_base = c / 4;
    peak_set = f_base .* ( 1:1:58 );

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

return

%% Part b - Double-tuned Expansion Chamber

return

%% Part c - Cascaded, Double-tuned Expansion Chamber

return

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


