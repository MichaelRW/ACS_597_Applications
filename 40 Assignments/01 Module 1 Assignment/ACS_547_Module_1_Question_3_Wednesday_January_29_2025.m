


%% Synopsis

% Question 3 - Bugle Recorder



%% Note(s)

% For the lowest frequency, use 1 duct sgement with an open-ended
% impedance (see the example of the horn in class).



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



%% Constants and Anonymous Functions

rho0 = 1.21;  % Density of air (kg per cubic-meter).
c = 343;  % Speed of sound in air (meters per second).


h_R_A = @( rho0, c , S, k, delta_mu, D, w, gamma, h, epsilon, M ) ...
    ( rho0*c/S )  * ( ( (k*delta_mu*D*w) / (2*S) )*( 1 + (gamma - 1)*sqrt(5/(3*gamma)) )  +  0.288*k*delta_mu*log10((4*S)/(pi*h^2))  +  (epsilon*S*k^2)/(2*pi)  +  0.7*M );
%
% See Equation 8.34 on page 479 of Bies et al (2024).



%% Define Shape

L_mouth_piece = 0.09;  % Meters

pipe.inner_diameter = 0.009;  % Meters
pipe.thickness = 0.004;  % Meters

% The recorder is unflanged.

hole_diameter = 0.006;  % Meters



%% Part a

% The estimated total length of the recorder is 0.325 meters.

% The estimated length of the pipe extension is 0.235 meters.

% Estimation was done by trial-and-error.

pipe_net_length = 0.325;
pipe_area = pi*0.009^2/4;

flanged = false;


frequency_set = 1:1:2e3;

nFreq = length( frequency_set );
    A = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );

    T_total = [ 1 0; 0 1 ];

    T_segment = duct_segment_transfer_matrix( f, rho0, c, pipe_net_length, pipe_area );

    T_total = T_segment * T_total;

    Z = open_end_impedance( f, rho0, c, 0, pipe_area, flanged );

    T11 = T_total(1, 1);  T12 = T_total(1, 2);    
        A( frequency_index ) = -10*log10( abs( T11 + T12 / Z )^2 );

end

A_parta = A;  clear A;

% figure( ); ...
%     plot( frequency_set, A_parta, 'LineWidth', 0.8 );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
%         legend( 'C5' );

% return

%% Part b

epsilon = 0.006 / 0.009;  % 0.67

a = 0.006 / 2;


L_o = a * ( 0.9326 - 0.6196*epsilon );  % Lecture 3, Slide 11
L_e = 0.004 + 2*L_o;
    Z_A = 1j * rho0 * (2 * pi * f) * L_e / ( pi*0.006^2/4 );

    
k = 2*pi*f/c;  % The wave number for the respective frequency.
S_hole = pi/4*(0.006)^2;  % squared-meters
mu = 1.83e-5;  % kg per meter-second;  online reference.
delta_mu = sqrt( (2 * mu ) / ( 2*pi*f * rho0 ) );
D = pi * 0.006;
w = 2*pi*f;
gamma = 1.4;
h = 0.003;  % Larger of the edge radius or delta_mu.
Mach_number = 0;
R_A = h_R_A( rho0, c, S_hole, k, delta_mu, D, w, gamma, h, epsilon, Mach_number );
%
% Z_A = Z_A + R_A
T_Branch = [ 1  0;  1/Z_A  1 ];


duct_lengths = 0.235/4 * ones( 4, 1 );


frequency_set = 0:0.1:2e3;

nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );

    T_total = [ 1 0; 0 1 ];

    T1 = duct_segment_transfer_matrix( f, rho0, c, duct_lengths(4), pipe_area );  % Duct - Outlet
    T2 = [ 1 0;  0 1 ];  % Hole
    T3 = duct_segment_transfer_matrix( f, rho0, c, duct_lengths(4), pipe_area );  % Duct
    T4 = [ 1 0;  0 1 ];  % Hole
    T5 = duct_segment_transfer_matrix( f, rho0, c, duct_lengths(4), pipe_area );  % Duct
    T6 = [ 1 0;  0 1 ];  % Hole
    T7 = duct_segment_transfer_matrix( f, rho0, c, duct_lengths(4), pipe_area );  % Duct
    T8 = duct_segment_transfer_matrix( f, rho0, c, 0.09, pipe_area );  % Duct - Inlet

    
    T_total =  T8 * T7 * T6 * T5 * T4 * T3 * T2 * T1 * T_total;

    Z = open_end_impedance( f, rho0, c, 0, pipe_area, flanged );
    
    T11 = T_total(1, 1);  T12 = T_total(1, 2);
        A( frequency_index ) = -10*log10( abs( T11 + T12 / Z )^2 );
            
end


A_partb = A;  clear A;

figure( ); ...
    plot( frequency_set, A_partb );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Amplification Versus Recorder Length' );

return

%% Part b Verification

% a = 0.009 / 2;  % Meters
%     L_o = 0.61*a;  % Slide 18 of Lecture 2 slide set.
% 
% f = 698;  % Hz
%     k = 2*pi*f/c;  % The wave number for the respective frequency.    
% 
% S = pi/4*(0.009)^2;  % squared-meters
% 
% 
% f = 0:1:5e3;
% 
% nFreq = length( f );
%     A = zeros( nFreq, 1 );
% 
% 
% L1 = 0.2
%     L2 = 0.325 - L1;
% 
% 
% for iFreq = 1:1:nFreq
% 
%     k = 2*pi*f(iFreq)/c;
% 
%     T_total = [ 1 0; 0 1 ];
% 
%     % End duct.
%     T_1 = [ ...
%     cos(k*L1),                           1j*rho0*c/S*sin(k*L1); ...
%     1j*S/(rho0*c)*sin(k*L1),      cos(k*L1) ...
%     ];
% 
% 
%     % Orifice side branch.
%     epsilon = 0.006 / 0.009;  % 0.67
%         a = 0.006 / 2;
% 
%     L_o = a * ( 0.9326 - 0.6196*epsilon );  % Lecture 3, Slide 11
%         L_e = 0.004 + 2*L_o;
%     %
%     Z_A = 1j * rho0 * 2 * pi * f(iFreq) * L_e / ( pi*0.006^2/4 );
%         T_Branch = [ 1  0;  1/Z_A  1 ];
%     %
%     % R_A is neglected (energy loss).
% 
% 
%     % Front duct.
%     T_2 = [ ...
%     cos(k*L2),                           1j*rho0*c/S*sin(k*L2); ...
%     1j*S/(rho0*c)*sin(k*L2),      cos(k*L2) ...
%     ];
% 
% 
%     T_total = T_2 * T_Branch * T_1 * T_total;    
% 
%     T11 = T_total(1, 1);  T12 = T_total(1, 2);
%         A( iFreq ) = -10*log10( abs( T11 + T12 / Z )^2 );
% 
% end
% 
% 
% figure( ); ...
%     plot( f, A );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
%     title( 'Amplification Versus Frequency' );



%% Plot Amplification Profiles

figure( ); ...
    plot( frequency_set, A );  hold on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
        legend( ...
                'Simple Expansion Chamber', ...
                'Double-tuned Expansion Chamber', ...
                'Cascaded Double-tuned Expansion Chamber', ...
                'Location', 'SouthOutside' );



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


