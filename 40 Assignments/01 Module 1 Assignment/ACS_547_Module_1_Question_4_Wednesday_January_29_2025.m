


%% Synopsis

% Problem 4, Intake Duct

% From ACS 597 - Noise Control Applications, Lecture 3 (Wednesday, January 22, 2025)

% Material from lectures 2 and 3.

% Assumptions:  source is high-impedence;  load is Z_open;  no mean flow;



%% To Do

% CHECK:  equation for volumetric flow rate;  

% Focus on interpretation;  "big picture" understanding.

% Do not repeat the prompt\question.

% Code should have "meat" of understanding;  avoid over-commenting.

% Code should use intuitive variable names.



%% Environment

% close all; clear; clc;
clear; clc;
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
duct_1.diameter_meters = 0.0254;  % 1 inch
duct_1.length_meters = 0.1524;  % 6 inches
duct_1.area = h_area( duct_1.diameter_meters );

% Outlet
duct_2.diameter_meters = 0.1016;  % 4 inches
duct_2.length_meters = 0.127;  % 5 inches
duct_2.area = h_area( duct_2.diameter_meters );
%
% Flanged.

%

flow_rate_cubic_meters_per_second = 1.04772 / 60;  % or 37 cubic-feet per minute

% return

%% Part a

% Calculate the Mach number of the flow in both pip sections.

duct_1.Mach = -1.0 * flow_rate_cubic_meters_per_second / (pi * duct_1.diameter_meters^2 / 4 ) / c;  % 0.100 unitless
duct_2.Mach = -1.0 * flow_rate_cubic_meters_per_second / (pi * duct_2.diameter_meters^2 / 4 ) / c;  % 0.00628 unitless

% return

%% Part b

% No flow.

outlet_flanged = true;  % Flanged end.


TEST_FLAG = 1;  % 1: right-to-left.

frequency_set = 0:0.1:2.5e3;
    nFreq = length( frequency_set );
        TL = zeros( nFreq, 1 );


for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );


    T_total = [ 1 0; 0 1 ];


    if ( TEST_FLAG == 1 )
        
        % Right-to-left.
        %
        T_outlet = duct_segment_transfer_matrix( f, rho0, c, duct_2.length_meters, duct_2.area );
        T_contraction = [ 1  0; 0  duct_2.area/duct_1.area ];
        T_inlet = duct_segment_transfer_matrix( f, rho0, c, duct_1.length_meters, duct_1.area );

        Z = open_end_impedance( f, rho0, c, duct_2.length_meters, duct_2.area, outlet_flanged );

    else

        % Left-to-right.
        %
        T_outlet = duct_segment_transfer_matrix( f, rho0, c, duct_1.length_meters, duct_1.area );
        T_contraction = [ 1  0; 0  duct_1.area/duct_2.area ];
        T_inlet = duct_segment_transfer_matrix( f, rho0, c, duct_2.length_meters, duct_2.area );

        Z = open_end_impedance( f, rho0, c, duct_1.length_meters, duct_1.area, outlet_flanged );

    end


    T_net = T_inlet * T_contraction * T_outlet * T_total;


    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);

    if ( TEST_FLAG == 1 )
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  duct_1.area*T12/(rho0*c)  +  (rho0*c)*T21/duct_2.area  +  T22 ) / 2 )^2 );
    else
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  duct_2.area*T12/(rho0*c)  +  (rho0*c)*T21/duct_1.area  +  T22 ) / 2 )^2 );
    end


end  % End:  for f = frequency_set

TL_part_b = TL;

% return

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

    T_net = T_inlet * T_contraction * T_outlet * T_total;
    
    
    % Z = open_end_impedance( f, rho0, c, duct_2.length_meters, duct_2.area, outlet_flanged );
    Z = open_end_impedance( f, rho0, c, duct_1.length_meters, duct_1.area, outlet_flanged );
    

    T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
        TL( frequency_index ) = 10 * log10( abs( ( T11  +  duct_1.area*T12/(rho0*c)  +  (rho0*c)*T21/duct_2.area  +  T22 ) / 2 )^2 );
 

end  % End:  for f = frequency_set

TL_part_c = TL;

% return

%% Plot

Y_LIMITS = [ 0  50 ];

figure( ); ...
    plot( frequency_set, TL_part_b );  hold on;
    plot( frequency_set, TL_part_c, '--' );
    % plot( frequency_set, TL_partc, 'LineStyle', '--' );  grid on;
    %     legend( ...
    %         'Simple Expansion Chamber', ...
    %         'Double-tuned Expansion Chamber', ...
    %         'Cascaded Double-tuned Expansion Chamber', ...
    %         'Location', 'SouthOutside' );
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





























