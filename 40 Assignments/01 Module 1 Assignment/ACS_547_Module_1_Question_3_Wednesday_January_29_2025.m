


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
set( 0, 'DefaultFigureWindowStyle', 'docked' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Constants and Anonymous Functions

rho0 = 1.21;  % Density of air (kg per cubic-meter).
c = 343;  % Speed of sound in air (meters per second).


h_R_A = @( rho0, c , S, k, delta_mu, D, w, h ) ...
    ( rho0*c/S )  *  (  ( (k*delta_mu*D*w) / (2*S) )*1.4364   +   0.288*k*delta_mu*log10((4*S)/(pi*h^2))   +   (0.5*S*k^2)/(2*pi)  );

% h_RA_term_1 = @( rho0, c , S, k, delta_mu, D, w )  ( rho0*c/S )  *  ( (k*delta_mu*D*w) / (2*S) *1.4364 );
h_RA_term_1 = @( rho0, c , S, k, delta_mu, D, w )  ( rho0*c/S )  *  ( (k * sqrt( (2*3.178e-5) / (rho0*w) ) * D * 0.004 ) / (2*S) *1.4364 );

h_RA_term_2 = @( rho0, c , S, k, delta_mu, D, w, h ) ...
    ( rho0*c/S )  *  0.288*k*delta_mu*log10((4*S)/(pi*h^2));

h_RA_term_3 = @( rho0, c , S, k, delta_mu, D, w, h ) ...
    ( rho0*c/S )  *  (0.5*S*k^2)/(2*pi);


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

epsilon = 0.006 / 0.004;  % Diameter of the hold divided by diameter of pipe section (1.5).

switch ( 3 )

    case 0  % Original Value

        L_o = 0.00001;  % Estimate
            L_e = 0.004 + 2*L_o;

        fprintf( 1, '\nZero - Percentage change in pipe thickness:  %3.1f%%.\n\n', ( L_e - 0.004 ) / 0.004 * 100 );

    case 1  % Ingard (2010)

        a = 0.006 / 2;

        L_o = ( 0.6*a + 0.85*a ) / 2;
            L_e = 0.004 + 2*L_o;

        fprintf( 1, '\nIngard (2010) - Percentage change in pipe thickness:  %3.1f%%.\n\n', ( L_e - 0.004 ) / 0.004 * 100 );

    case 2  % Kurze and Riedel (2013)

        e = epsilon^2;

        a = 0.006 / 2;

        L_o = pi*a*( 1 - 1.47*e^0.5 + 0.47*e^1.5 );
            L_e = 0.004 + 2*L_o;

        fprintf( 1, '\nKurze and Riedel (2013) - Percentage change in pipe thickness:  %3.1f%%.\n\n', ( L_e - 0.004 ) / 0.004 * 100 );

    case 3  % Ji (2005)

        a = 0.006 / 2;

        L_o = a*( 0.9326 - 0.6196*epsilon);
            L_e = 0.004 + 2*L_o;

        fprintf( 1, '\nJi (2005) - Percentage change in pipe thickness:  %3.1f%%.\n\n', ( L_e - 0.004 ) / 0.004 * 100 );

    otherwise
        error ( '*** Invalid SWITCH Index ***' );

end


% All holes covered.
% 523 Hz - Length of pipe is 0.235 meters.

duct_lengths = 0.235/4 * ones( 4, 1 );  % 523 Hz


% frequency_set = 0:1:2e3;
frequency_set = 0:1:1e3;

nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );
    ZA_real = zeros( nFreq, 1 );  ZA_imaginary = zeros( nFreq, 1 );
    term_1_v = zeros( nFreq, 1 );  term_2_v = zeros( nFreq, 1 );  term_3_v = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );

    T_total = [ 1 0; 0 1 ];


    Z_A = 1j * rho0 * (2 * pi * f) * L_e / ( pi*0.006^2/4 );
        ZA_imaginary( frequency_index ) = Z_A;
    %
    % R_A = h_R_A( rho0, c, pi*(0.006)^2/4, 2*pi*f/c, sqrt( (2 * 1.83e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.006, 2*pi*f, 0.3 );
    
     term_1 = h_RA_term_1( rho0, c, pi*(0.006)^2/4, 2*pi*f/c, sqrt( (2 * 1.83e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.006, 2*pi*f );
        term_1_v( frequency_index ) = term_1;


     term_2 = h_RA_term_2( rho0, c, pi*(0.006)^2/4, 2*pi*f/c, sqrt( (2 * 1.83e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.006, 2*pi*f, 0.3 );
        term_2_v( frequency_index ) = term_2;
     term_3 = h_RA_term_3( rho0, c, pi*(0.006)^2/4, 2*pi*f/c, sqrt( (2 * 1.83e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.006, 2*pi*f, 0.3 );
        term_3_v( frequency_index ) = term_3;
    %
    R_A = term_1 + term_2 + term_3;

    % R_A =  R_A * 1e-5;
        ZA_real( frequency_index ) = R_A;
    %
    Z_A = Z_A + R_A;
        T_Hole = [ 1  0;  1/Z_A  1 ];



    switch ( 3 )
        
        case 0  % All holes covered.
            T1 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct - Outlet
            T2 = [ 1 0;  0 1 ];  % Hole
            T3 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct
            T4 = [ 1 0;  0 1 ];  % Hole
            T5 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct
            T6 = [ 1 0;  0 1 ];  % Hole
            T7 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct

        case 1  % Hole 1
            T1 = duct_segment_transfer_matrix( f, rho0, c, 0.05875 + 0.0290, pipe_area );  % Duct - Outlet
            T2 = T_Hole;
            T3 = duct_segment_transfer_matrix( f, rho0, c, 0.05875 - 0.0290, pipe_area );  % Duct
            T4 = [ 1 0;  0 1 ];  % Hole
            T5 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct
            T6 = [ 1 0;  0 1 ];  % Hole
            T7 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct

        case 2  % Hole 2
            T1 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct - Outlet
            T2 = [ 1 0;  0 1 ];  % Hole
            T3 = duct_segment_transfer_matrix( f, rho0, c, 0.05875 + 0.0230, pipe_area );  % Duct
            T4 = T_Hole;
            T5 = duct_segment_transfer_matrix( f, rho0, c, 0.05875 - 0.0230, pipe_area );  % Duct
            T6 = [ 1 0;  0 1 ];  % Hole
            T7 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct

        case 3  % Hole 3

            % OFFSET = 0.005;  % 882
            % OFFSET = 0.01;  % 862
            % OFFSET = 0.02;  % 823

            T1 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct - Outlet
            T2 = [ 1 0;  0 1 ];  % Hole
            T3 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct
            T4 = [ 1 0;  0 1 ];  % Hole
            T5 = duct_segment_transfer_matrix( f, rho0, c, 0.05875 + OFFSET, pipe_area );  % Duct
            T6 = T_Hole;
            T7 = duct_segment_transfer_matrix( f, rho0, c, 0.05875 - OFFSET, pipe_area );  % Duct

    end
    

    T8 = duct_segment_transfer_matrix( f, rho0, c, 0.09, pipe_area );  % Duct - Inlet
    
    T_total =  T8 * T7 * T6 * T5 * T4 * T3 * T2 * T1 * T_total;

    Z = open_end_impedance( f, rho0, c, 0, pipe_area, flanged );
    
    T11 = T_total(1, 1);  T12 = T_total(1, 2);
        A( frequency_index ) = -10*log10( abs( T11 + T12 / Z )^2 );
            
end


[ max_value, max_index ] = max( A );
    frequency_set( max_index )

figure( ); ...
    plot( frequency_set, A );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Amplification Versus Recorder Length' );


% ind_term_2 = ~ismissing( term_2_v );
% ind_term_3 = ~ismissing( term_3_v );
% 
% figure( ); ...
%     subplot( 2, 3, 1 ); ...
%         plot( frequency_set, ZA_real  );  title( 'Real-part of ZA' );
%     subplot( 2, 3, 2 ); ...
%         plot( frequency_set, term_1_v, 'LineStyle', '--' );  title( 'Term 1' );
%     subplot( 2, 3, 3 ); ...
%         plot( frequency_set( ind_term_2 ), term_2_v( ind_term_2 ), 'LineStyle', ':' );  title( 'Term 2' );
%     subplot( 2, 3, 4 ); ...
%         plot( frequency_set( ind_term_3 ), term_3_v( ind_term_3 ) );  title( 'Term 3' );
%     subplot( 2, 3, 5 ); ...
%         plot( frequency_set, imag( ZA_imaginary )  );  grid on;
%
    % legend( 'Real Part', 'Term 1', 'Term 2', 'Term 3', 'Imaginary Par' );
    % xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    % set( gca, 'YScale', 'Log' );

return

%% Plot Note Set

load( 'A_C5_Data.mat' );  % Variable(s):  A_C5


figure( ); ...
    plot( frequency_set, A_C5 );  hold on;
        text( 523, 25, 'C5' );
    plot( frequency_set, A_C5 );
        text( 523, 25, 'C5' );
    plot( frequency_set, A_C5 );
        text( 523, 25, 'C5' );
    plot( frequency_set, A_C5 );  grid on;
        text( 523, 25, 'C5' );
        %
        legend( 'C5', 'F5', 'A5', 'C6', 'Location', 'West' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Amplification Versus Recorder Length' );



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


