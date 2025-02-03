


%% Synopsis

% Question 3 - Bugle Recorder



%% Environment

close all; clear; clc;
% restoredefaultpath;

set( 0, 'DefaultFigureWindowStyle', 'docked' );

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


h_RA = @( rho0, c , S, k, delta_mu, D, w, gamma, h, epsilon, M ) ...
    ( rho0*c/S )  * ( ( (k*delta_mu*D*w) / (2*S) )*( 1 + (gamma - 1)*sqrt(5/(3*gamma)) )  +  0.288*k*delta_mu*log10((4*S)/(pi*h^2))  +  (epsilon*S*k^2)/(2*pi)  +  0.7*M );
%
% See Equation 8.34 on page 479 of Bies et al (2024).



%% Define Shape

mouth_piece_length = 0.09;  % Meters

pipe_diameter = 0.009;  % Meters
pipe_thickness = 0.004;  % Meters

% The recorder is unflanged.

hole_diameter = 0.006;  % Meters



%% Part a


% The total length of the recorder, including the 0.09 meter long mouthpiece, is L.

a = 0.009 / 2;  % Meters
    L_o = 0.61*a;  % Slide 18 of Lecture 2 slide set.

f = 523;  % Hz
    k = 2*pi*f/c;  % The wave number for the respective frequency.    

S = pi/4*(0.009)^2;  % squared-meters


duct_lengths = 0.11262 * ones( 4, 1 );  % 523 Hz.
    k = 2*pi*f/c;  % The wave number for the respective frequency.    


frequency_set = 0:0.1:1e3;

nFreq = length( frequency_set );
    TL = zeros( nFreq, 1 );

for frequency_index = 1:1:nFreq

    f = frequency_set( frequency_index );

    T_total = [ 1 0; 0 1 ];


    % Duct segments.
    T_segment_1 = [ ...
        cos(k* duct_lengths( 1 ) ),  1j*rho0*c/S*sin(k* duct_lengths( 1 ) ); ...
        1j*S/(rho0*c)*sin(k* duct_lengths( 1 ) ),  cos(k* duct_lengths( 1 ) )  ];
    %
    T_segment_2 = [ ...
        cos(k* duct_lengths( 2 ) ),  1j*rho0*c/S*sin(k* duct_lengths( 2 ) ); ...
        1j*S/(rho0*c)*sin(k* duct_lengths( 2 ) ),  cos(k* duct_lengths( 2 ) )  ];
    %
    T_segment_3 = [ ...
        cos(k* duct_lengths( 3 ) ),  1j*rho0*c/S*sin(k* duct_lengths( 3 ) ); ...
        1j*S/(rho0*c)*sin(k* duct_lengths( 3 ) ),  cos(k* duct_lengths( 3 ) )  ];
    %
    T_segment_4 = [ ...
        cos(k* ( duct_lengths( 4 ) + 90e-3 ) ),  1j*rho0*c/S*sin(k* ( duct_lengths( 4 ) + 90e-3 ) ); ...
        1j*S/(rho0*c)*sin(k* ( duct_lengths( 4 ) + 90e-3 ) ),  cos(k* ( duct_lengths( 4 ) + 90e-3 ) )  ];

        
    % Holes
    T_hole_1 = [ 1 0;  0 1 ];
    T_hole_2 = [ 1 0;  0 1 ];
    T_hole_3 = [ 1 0;  0 1 ];

    % L_e = L + L_o;


    % Net transfer function.
    T_total = T_segment_4 * T_hole_3 * T_segment_3 * T_hole_2 * T_segment_2 * T_hole_1 * T_segment_1 * T_total;
        T11 = T_total(1, 1);  T12 = T_total(1, 2);

    Z = open_end_impedance( f, rho0, c, duct_lengths(4), S, 0 );  % Unflanged.
        A( frequency_index ) = -10*log10( abs( T11 + T12 / Z )^2 );
            
end


figure( ); ...
    plot( frequency_set, A );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Amplification Versus Recorder Length' );

[ max_value, max_index ] = max( A );
    frequency_set( max_index )

return

%% Part b

% L_net = 0.325;  % Meters - First Peak
L_net = 0.653;  % Meters - Second Peak

a = 0.009 / 2;  % Meters
    L_o = 0.61*a;  % Slide 18 of Lecture 2 slide set.

f = 698;  % Hz
    k = 2*pi*f/c;  % The wave number for the respective frequency.    

S = pi/4*(0.009)^2;  % squared-meters


test_lengths = 0:0.001:0.5;
    test_lengths = L_net - test_lengths;

nLengths = length( test_lengths );
    A = zeros( nLengths, 1 );

for iLength = 1:1:nLengths

    L = test_lengths(iLength);
        L_duct_2 = L;
        L_duct_1 = L_net - L_duct_2;

    T_total = [ 1 0; 0 1 ];


    % End section.
    T_1 = [ ...
    cos(k*L_duct_1),                           1j*rho0*c/S*sin(k*L_duct_1); ...
    1j*S/(rho0*c)*sin(k*L_duct_1),      cos(k*L_duct_1) ...
    ];


    % Orifice side branch.
    epsilon = 0.006 / 0.009;  % 0.67
        a = 0.006 / 2;
    %
    L_o = a * ( 0.9326 - 0.6196*epsilon );  % Lecture 3, Slide 11
        L_e = 0.004 + 2*L_o;
    %
    Z_A = 1j * rho0 * (2 * pi * f) * L_e / ( pi*0.006^2/4 );
    %
    k = 2*pi*f/c;  % The wave number for the respective frequency.
    S_hole = pi/4*(0.006)^2;  % squared-meters
    mu = 1.83e-5;  % kg per meter-second;  online reference.
        delta_mu = sqrt( (2 * mu ) / ( 2*pi*f * rho0 ) );
    D = pi * 0.006;
    w = 2*pi*f;
    gamma = 1.4;
    h = 0.003;  % Larger of the edge radius or delta_mu.
    Mach_number = 0;
        R_A = h_RA( rho0, c, S_hole, k, delta_mu, D, w, gamma, h, epsilon, Mach_number );
        %
        % Z_A = Z_A + R_A
            T_Branch = [ 1  0;  1/Z_A  1 ];


    % Section next to mouthpiece.
    T_2 = [ ...
    cos(k*L_duct_2),                           1j*rho0*c/S*sin(k*L_duct_2); ...
    1j*S/(rho0*c)*sin(k*L_duct_2),      cos(k*L_duct_2) ...
    ];


    T_total = T_2 * T_Branch * T_1 * T_total;
        T11 = T_total(1, 1);  T12 = T_total(1, 2);

    A( iLength ) = -10*log10( abs( T11 + T12 / Z )^2 );

end


figure( ); ...
    plot( test_lengths * 1e3, A );  grid on;
        set( gca, 'XDir', 'reverse' );
    xlabel( 'Offset from End of 162.6 mm Length Recorder [mm]' );  ylabel( 'Amplitude [dB]' );
    title( 'Amplification Versus Offset from End of Recorder' );

return

%% Part b Verification

a = 0.009 / 2;  % Meters
    L_o = 0.61*a;  % Slide 18 of Lecture 2 slide set.

f = 698;  % Hz
    k = 2*pi*f/c;  % The wave number for the respective frequency.    

S = pi/4*(0.009)^2;  % squared-meters


f = 0:1:5e3;

nFreq = length( f );
    A = zeros( nFreq, 1 );


L1 = 0.2
    L2 = 0.325 - L1;


for iFreq = 1:1:nFreq

    k = 2*pi*f(iFreq)/c;

    T_total = [ 1 0; 0 1 ];

    % End duct.
    T_1 = [ ...
    cos(k*L1),                           1j*rho0*c/S*sin(k*L1); ...
    1j*S/(rho0*c)*sin(k*L1),      cos(k*L1) ...
    ];


    % Orifice side branch.
    epsilon = 0.006 / 0.009;  % 0.67
        a = 0.006 / 2;

    L_o = a * ( 0.9326 - 0.6196*epsilon );  % Lecture 3, Slide 11
        L_e = 0.004 + 2*L_o;
    %
    Z_A = 1j * rho0 * 2 * pi * f(iFreq) * L_e / ( pi*0.006^2/4 );
        T_Branch = [ 1  0;  1/Z_A  1 ];
    %
    % R_A is neglected (energy loss).


    % Front duct.
    T_2 = [ ...
    cos(k*L2),                           1j*rho0*c/S*sin(k*L2); ...
    1j*S/(rho0*c)*sin(k*L2),      cos(k*L2) ...
    ];


    T_total = T_2 * T_Branch * T_1 * T_total;    

    T11 = T_total(1, 1);  T12 = T_total(1, 2);
        A( iFreq ) = -10*log10( abs( T11 + T12 / Z )^2 );

end


figure( ); ...
    plot( f, A );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Amplitude [dB]' );
    title( 'Amplification Versus Frequency' );



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


