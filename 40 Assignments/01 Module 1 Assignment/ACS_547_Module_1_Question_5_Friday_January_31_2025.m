

%% Synopsis

% Problem 5, Intake Duct Silencer

% From ACS 597 - Noise Control Applications, Lecture 3 (Wednesday, January 22, 2025)

% Material from lecture 4.

% Assumptions:  source is high-impedence;  load is Z_open;  no mean flow;

% 



%% To Do

% CHECK:  equation for volumetric flow rate;  

% Focus on interpretation;  "big picture" understanding.

% Do not repeat the prompt\question.

% Code should have "meat" of understanding;  avoid over-commenting.

% Code should use intuitive variable names.



%% Environment

close all; clear; clc;
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

% The -35 dB transmission loss dip\notch is at about 1,150 Hz.
f = 1150;

% Total attenuation required is +35 dB to make it 0 dB.

l = 0.0381;  % meters (1.5 inches)
h = 0.0127;  % meters (0.5 inches)

% The total attenuation of the lining from Figure 8.37 (Bies et. al., Fifth Edition),
%
% m = ( h_area( 0.1016 ) - h_area( 0.0254 ) ) / h_area( 0.0254 )
%
% ASSUMPTION:  Cross-section area of lined duct is the same as the inlet.
m = 1;

k = 2*pi*1.15e3;

length_of_expansion_chamber = 0.127;  % meters

% kL = k*length_of_expansion_chamber

% ASSUMPTION:
%
% Expansion ratio is m = 1;
% Assume peak difference is 10 dB.
% Total attunation of lining is 10 dB.
%
% The attenuation rate is about 10 dB / 0.127 meters or 78.7 dB per meter.

% return

%% Part b

l = 0.0381;  % meters (1.5 inches)

h = 0.011255;  % meters
    h_validate = 0.5*sqrt( ( pi*( duct_2.diameter_meters - 2*l  )^2 ) / 4 );  % Same value.



%% Part c

% The liner thickness ratio is,
l / h;  % 3.3852 unitless

% The normalized frequency is,
( 2 * h ) / ( 343 / f );  % 0.075471 unitless



%% Part d

% Assume the attenuation rate from Part a is 18 dB per meter.

% Use bottom, right subplot (16).

% The approximate resistivity parameter is 16.



%% Part e

% Calculate the flow resistivity.

R1 = 16 * rho0*c / l;  % 1.74e5 kg / m^3*s



%% Placeholder

% % Flow present (use Mach numbers).
% % Helmholtz resonator in place (between lefthand duct and expansion).
% 
% outlet_flanged = true;  % Flanged end.
% 
% 
% % Resonance
% w_o = 2*pi*136.6;  % Estimated from plot.
% 
% 
% % helmholtz_diameter_cavity = 
% % helmholtz_diameter_neck = 1e-3;
% %     helmholtz_L01 = 0.82 * ( 1 - 1.33*(helmholtz_diameter_neck/helmholtz_diameter_cavity ) );
% %     %
% %     epsilon = helmholtz_diameter_cavity / duct_1.length_meters;
% %     % helmholtz_L02 = 
% % 
% % helmholtz_volume = 1e-3;
% % 
% % % keyboard
% % 
% % helmholtz_L_neck = 1e-3;
% % Q = 2;
% % 
% % R_A = rho0*c / Q * sqrt( L_e / ( pi*helmholtz_diameter_neck^2/4 * helmholtz_volume ) );
% 
% 
% 
% frequency_set = 0:0.1:2.5e3;
%     nFreq = length( frequency_set );
%         TL = zeros( nFreq, 1 );
% 
% 
% for frequency_index = 1:1:nFreq
% 
%     f = frequency_set( frequency_index );
% 
% 
%     T_total = [ 1 0; 0 1 ];
% 
%     T_outlet = duct_segment_transfer_matrix_flow( f, rho0, c, duct_2.length_meters, duct_2.area, duct_2.Mach );
% 
%     T_expansion = duct_expansion_connection_transfer_matrix( rho0, c, duct_2.area, duct_1.area, duct_1.Mach );
% 
% 
%     % Z_A = 1j*rho0*2*pi*f(frequency_index) * L_e / ( pi*helmholtz_neck_diameter^2/4 )  -  1j*rho0*c^2/(helmholtz_volume*2*pi*f(frequency_index))  +  R_A;
%         % T_Helmholtz = [ 1  0;  1/Z_A  1 ];
% 
% 
%     T_inlet = duct_segment_transfer_matrix_flow( f, rho0, c, duct_1.length_meters, duct_1.area, duct_1.Mach );
% 
% 
%     % T_net = T_inlet * T_Helmholtz * T_expansion * T_outlet * T_total;
%     T_net = T_inlet * T_expansion * T_outlet * T_total;
% 
% 
%     Z = open_end_impedance( f, rho0, c, duct_2.length_meters, duct_2.area, outlet_flanged );
% 
%     T11 = T_net(1, 1);  T12 = T_net(1, 2);  T21 = T_net(2, 1);  T22 = T_net(2, 2);
%         TL( frequency_index ) = 10 * log10( abs( ( T11  +  duct_1.area*T12/(rho0*c)  +  (rho0*c)*T21/duct_2.area  +  T22 ) / 2 )^2 );
% 
% 
% end  % End:  for f = frequency_set
% 
% TL_part_d = TL;

% return

%% Plot

% Y_LIMITS = [ 0  50 ];
% 
% figure( ); ...
%     plot( frequency_set, TL_part_d );  hold on;
%     % plot( frequency_set, TL_part_c, '-' );
%     % plot( frequency_set, TL_part_d, '-.' );  grid on;
%     % legend( ...
%     %     'No Flow - Part b', ...
%     %     'Flow - Part c', ...
%     %     'Flow and Resonator - Part d', ...
%     %     'Location', 'SouthOutside' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Transmission Loss [dB]' );
%     title( 'Transmission Loss Profiles' );
%     %
%     Ax = gca;
%         Ax.XAxis.TickLabelInterpreter = 'latex';
%         Ax.YAxis.TickLabelInterpreter = 'latex';
%     %
%     % axis( [ -50  5e3+50  Y_LIMITS ] );
%     %
%     if ( PRINT_FIGURES == 1 )
%         exportgraphics( gcf, 'Figure TL All Profiles.pdf', 'Append', true );
%     end

% return

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





























