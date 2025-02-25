


%% Synopsis

% Question 3 - Bugle Recorder



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Define Constants and Anonymous Functions

rho0 = 1.21;  % Density of air (kg per cubic-meter).
c = 343;  % Speed of sound in air (meters per second).

h_RA_term_1 = @( rho0, c , S, k, delta_mu, D, w )  ( rho0*c/S )  *  ( (k * sqrt( (2*3.178e-5) / (rho0*w) ) * D * 0.004 ) / (2*S) *1.4364 );
h_RA_term_2 = @( rho0, c , S, k, delta_mu, D, w, h )  ( rho0*c/S )  *  0.288*k*3.178e-5*log10((4*S)/(pi*h^2));
h_RA_term_3 = @( rho0, c , S, k, delta_mu, D, w, h )  ( rho0*c/S )  *  (0.5*S*k^2)/(2*pi);
%
% See Equation 8.34 on page 479 of Bies et al (2024).



%% Problem 3a

% The estimated total length of the recorder is 0.325 meters.

% Estimation of hole location was done by trial-and-error.

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

%% Problem 3b

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


duct_lengths = 0.235/4 * ones( 4, 1 );  % 523 Hz, all holes covered.


frequency_set = 0:1:2e3;

nFreq = length( frequency_set );
    A = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );

    T_total = [ 1 0; 0 1 ];


    Z_A = 1j * rho0 * (2 * pi * f) * L_e / ( pi*0.006^2/4 );
    %
    term_1 = h_RA_term_1( rho0, c, pi*(0.006)^2/4, 2*pi*f/c, sqrt( (2 * 3.178e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.006, 2*pi*f );
    term_2 = h_RA_term_2( rho0, c, pi*(0.006)^2/4, 2*pi*f/c, sqrt( (2 * 3.178e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.006, 2*pi*f, 0.3 );
    term_3 = h_RA_term_3( rho0, c, pi*(0.006)^2/4, 2*pi*f/c, sqrt( (2 * 3.178e-5 ) / ( 2*pi*f * rho0 ) ), pi * 0.006, 2*pi*f, 0.3 );
        R_A = term_1 + term_2 + term_3;
    %
    Z_A = Z_A + R_A;
        T_Hole = [ 1  0;  1/Z_A  1 ];


    switch ( 0 )

        case 0  % All holes covered.  523 Hz
            T1 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct - Outlet
            T2 = [ 1 0;  0 1 ];  % Hole
            T3 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct
            T4 = [ 1 0;  0 1 ];  % Hole
            T5 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct
            T6 = [ 1 0;  0 1 ];  % Hole
            T7 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct

            % Check
            % 4 * 0.05875 = 0.235;

        case 1  % Hole 1 uncovered.  698 Hz
            T1 = duct_segment_transfer_matrix( f, rho0, c, 0.08775, pipe_area );  % Duct - Outlet
            T2 = T_Hole;
            T3 = duct_segment_transfer_matrix( f, rho0, c, 0.02975, pipe_area );  % Duct
            T4 = [ 1 0;  0 1 ];  % Hole
            T5 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct
            T6 = [ 1 0;  0 1 ];  % Hole
            T7 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct

            % Check
            % 0.08775 + 0.02975 + 2*0.05875 = 0.235;

        case 2  % Hole 2 uncovered.  880 Hz

            offset = 0.00625;  % 874

            T1 = duct_segment_transfer_matrix( f, rho0, c, 0.08775, pipe_area );  % Duct - Outlet
            T2 = T_Hole;
            T3 = duct_segment_transfer_matrix( f, rho0, c, 0.0505, pipe_area );  % Duct
            T4 = T_Hole;
            T5 = duct_segment_transfer_matrix( f, rho0, c, 0.038, pipe_area );  % Duct
            T6 = [ 1 0;  0 1 ];  % Hole
            T7 = duct_segment_transfer_matrix( f, rho0, c, 0.05875, pipe_area );  % Duct

            % Check
            % 0.08775 + 0.0505 + 0.038 + 0.05875 = 0.235

        case 3  % Hole 3 uncovered.  1,046 Hz

            T1 = duct_segment_transfer_matrix( f, rho0, c, 0.08775, pipe_area );  % Duct - Outlet
            T2 = T_Hole;
            T3 = duct_segment_transfer_matrix( f, rho0, c, 0.0505, pipe_area );  % Duct
            T4 = T_Hole;
            T5 = duct_segment_transfer_matrix( f, rho0, c, 0.02975, pipe_area );  % Duct
            T6 = T_Hole;
            T7 = duct_segment_transfer_matrix( f, rho0, c, 0.067, pipe_area );  % Duct

            % Check
            % 0.08775 + 0.0505 + 0.02975 + 0.067 = 0.235 m

    end

        T8 = duct_segment_transfer_matrix( f, rho0, c, 0.09, pipe_area );  % Duct - Inlet

        T_total =  T8 * T7 * T6 * T5 * T4 * T3 * T2 * T1 * T_total;

        Z = open_end_impedance( f, rho0, c, 0, pipe_area, flanged );

        T11 = T_total(1, 1);  T12 = T_total(1, 2);
            A( frequency_index ) = -10*log10( abs( T11 + T12 / Z )^2 );
            
end


% figure( ); ...
%     plot( frequency_set, A );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );



%% Plot Note Set

load( 'A_C5_Data.mat' );  % Variable(s):  A_C5
load( 'A_F5_Data.mat' );  % Variable(s):  A_F5
load( 'A_A5_Data.mat' );  % Variable(s):  A_A5
load( 'A_C6_Data.mat' );  % Variable(s):  A_C6


h_figure_1 = figure( ); ...
    plot( frequency_set, A_C5 );  hold on;
        text( 523, 25, 'C5' );
    plot( frequency_set, A_F5 );
        text( 698, -1, 'F5' );
    plot( frequency_set, A_A5 );
        text( 880, -14, 'A5' );
    plot( frequency_set, A_C6 );  grid on;
        text( 1042, -23, 'C6' );
        %
        legend( 'C5, 523 Hz', 'F5, 698 Hz', 'A5, 880 Hz', 'C6, 1046 Hz', 'Location', 'NorthWest' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Notes Spectrums for a Bugle Recorder' );



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
        exportgraphics( h_figure_1, 'Assignment 1 - Question 3 Bugle Recorder Note Spectrums.pdf', 'Append', true );
end

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


