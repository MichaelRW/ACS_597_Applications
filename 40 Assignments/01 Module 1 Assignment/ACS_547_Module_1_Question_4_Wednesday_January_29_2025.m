


%% Synopsis

% Question 4 - Intake Duct



%% Environment

% close all; clear; clc;
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

h_area = @( diameter)  pi * diameter^2 / 4;



%% Define Shape

% Source
duct_1.diameter_meters = 0.0254;  % meters
    duct_1.area = h_area( duct_1.diameter_meters );  % squared-meters
duct_1.length_meters = 0.1524;  % meters
    

% Outlet
duct_2.diameter_meters = 0.1016;  % 4 inches
    duct_2.area = h_area( duct_2.diameter_meters );
duct_2.length_meters = 0.127;  % 5 inches

% Flanged.

flow_rate_cubic_meters_per_second = 1.04772 / 60;  % or 37 cubic-feet per minute



%% Part a - Mach Numbers

% Calculate the Mach number of the flow in both pip sections.

duct_1.Mach = -1.0 * flow_rate_cubic_meters_per_second / ( (pi * duct_1.diameter_meters^2 / 4 ) * c );  % 0.100 unitless
duct_2.Mach = -1.0 * flow_rate_cubic_meters_per_second / ( (pi * duct_2.diameter_meters^2 / 4 ) * c );  % 0.00628 unitless



%% Part b - No Flow in Intake System

outlet_flanged = true;  % Flanged end.


frequency_set = 0:1:2.5e3;
    nFreq = length( frequency_set );
        TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );

    T_total = [ 1 0; 0 1 ];

    T1 = duct_segment_transfer_matrix( f, rho0, c, duct_2.length_meters, duct_2.area );
    T2 = [ 1  0; 0  duct_2.area/duct_1.area ];
    T3 = duct_segment_transfer_matrix( f, rho0, c, duct_1.length_meters, duct_1.area );

    T_net = T3 * T2 * T1 * T_total;

    Z = open_end_impedance( f, rho0, c, duct_2.length_meters, duct_2.area, outlet_flanged );

    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  duct_1.area*T12/(rho0*c)  +  (rho0*c)*T21/duct_2.area  +  T22 ) / 2 )^2 );

end


TL_part_b = TL;


figure( ); ...
    plot( frequency_set, TL_part_b );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Transmission Loss [dB]' );
    title( 'Transmission Loss Profile for Instake System with No Flow' );

return

%% Part c

% Flow present (use Mach numbers).

outlet_flanged = true;  % Flanged end.


frequency_set = 0:0.1:2.5e3;
    nFreq = length( frequency_set );
        TL = zeros( nFreq, 1 );


for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );


    T_total = [ 1 0; 0 1 ];

    T_outlet = duct_segment_transfer_matrix_flow( f, rho0, c, duct_2.length_meters, duct_2.area, duct_2.Mach );

    T_expansion = duct_expansion_connection_transfer_matrix( rho0, c, duct_2.area, duct_1.area, duct_1.Mach );
    
    T_inlet = duct_segment_transfer_matrix_flow( f, rho0, c, duct_1.length_meters, duct_1.area, duct_1.Mach );

    T_net = T_inlet * T_expansion * T_outlet * T_total;
    
    
    Z = open_end_impedance( f, rho0, c, duct_2.length_meters, duct_2.area, outlet_flanged );
    
    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  duct_1.area*T12/(rho0*c)  +  (rho0*c)*T21/duct_2.area  +  T22 ) / 2 )^2 );
 

end  % End:  for f = frequency_set

TL_part_c = TL;

% return

%% Part d

% Flow present (use Mach numbers).
% Helmholtz resonator in place (between lefthand duct and expansion).

outlet_flanged = true;  % Flanged end.


% Resonance
w_o = 2*pi*136.6;  % Estimated from plot.


% helmholtz_diameter_cavity = 
% helmholtz_diameter_neck = 1e-3;
%     helmholtz_L01 = 0.82 * ( 1 - 1.33*(helmholtz_diameter_neck/helmholtz_diameter_cavity ) );
%     %
%     epsilon = helmholtz_diameter_cavity / duct_1.length_meters;
%     % helmholtz_L02 = 
% 
% helmholtz_volume = 1e-3;
% 
% % keyboard
% 
% helmholtz_L_neck = 1e-3;
% Q = 2;
% 
% R_A = rho0*c / Q * sqrt( L_e / ( pi*helmholtz_diameter_neck^2/4 * helmholtz_volume ) );



frequency_set = 0:0.1:2.5e3;
    nFreq = length( frequency_set );
        TL = zeros( nFreq, 1 );


for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );


    T_total = [ 1 0; 0 1 ];

    T_outlet = duct_segment_transfer_matrix_flow( f, rho0, c, duct_2.length_meters, duct_2.area, duct_2.Mach );

    T_expansion = duct_expansion_connection_transfer_matrix( rho0, c, duct_2.area, duct_1.area, duct_1.Mach );


    % Z_A = 1j*rho0*2*pi*f(frequency_index) * L_e / ( pi*helmholtz_neck_diameter^2/4 )  -  1j*rho0*c^2/(helmholtz_volume*2*pi*f(frequency_index))  +  R_A;
        % T_Helmholtz = [ 1  0;  1/Z_A  1 ];
    

    T_inlet = duct_segment_transfer_matrix_flow( f, rho0, c, duct_1.length_meters, duct_1.area, duct_1.Mach );


    % T_net = T_inlet * T_Helmholtz * T_expansion * T_outlet * T_total;
    T_net = T_inlet * T_expansion * T_outlet * T_total;
    

    Z = open_end_impedance( f, rho0, c, duct_2.length_meters, duct_2.area, outlet_flanged );
    
    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  duct_1.area*T12/(rho0*c)  +  (rho0*c)*T21/duct_2.area  +  T22 ) / 2 )^2 );
 

end  % End:  for f = frequency_set

TL_part_d = TL;

% return

%% Plot

Y_LIMITS = [ 0  50 ];

figure( ); ...
    plot( frequency_set, TL_part_b );  hold on;
    plot( frequency_set, TL_part_c, '-' );
    plot( frequency_set, TL_part_d, '-.' );  grid on;
    legend( ...
        'No Flow - Part b', ...
        'Flow - Part c', ...
        'Flow and Resonator - Part d', ...
        'Location', 'SouthOutside' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Transmission Loss Profiles' );
    %
    Ax = gca;
        Ax.XAxis.TickLabelInterpreter = 'latex';
        Ax.YAxis.TickLabelInterpreter = 'latex';
    %
    % axis( [ -50  5e3+50  Y_LIMITS ] );
    %
    if ( PRINT_FIGURES == 1 )
        exportgraphics( gcf, 'Figure TL All Profiles.pdf', 'Append', true );
    end

return

%% Part c

return


%% Part d

return

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





























