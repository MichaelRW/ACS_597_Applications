


%% Synopsis

% Bugle Recorder

% From ACS 597 - Noise Control Applications, Lecture 3 (Wednesday, January 22, 2025)

% Material from lectures 2 and 3.

% Assumptions:  source is high-impedence;  load is Z_open;  no mean flow;



%% To Do

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



%% Define Shape

L_mouth_piece = 0.09;  % Meters

pipe.inner_diameter = 0.009;  % Meters
pipe.thickness = 0.004;  % Meters

% The recorder is unflanged.

hole_diameter = 0.006;  % Meters



%% Part a

% % Determine the length of the recorder to produce 523 Hz.
% 
% % The total length of the recorder, including the 0.09 meter long mouthpiece, is L.
% 
% a = 0.009 / 2;  % Meters
%     L_o = 0.61*a;  % Slide 18 of Lecture 2 slide set.
% 
% f = 523;  % Hz
%     k = 2*pi*f/c;  % The wave number for the respective frequency.    
% 
% S = pi/4*(0.009)^2;  % squared-meters
% 
% 
% test_lengths = 0:0.0001:0.1;
%     test_lengths = test_lengths + 0.09;
% 
% 
% nLengths = length( test_lengths );
%     A = zeros( nLengths, 1 );
% 
% 
% for iLength = 1:1:nLengths
% 
%     L = test_lengths(iLength);
% 
%     T_total = [ 1 0; 0 1 ];
% 
%     L_e = L + L_o;
%         Z = 1j * rho0 * c / S * tan( k* L_e );
% 
%     T = [ ...
%     cos(k*L),                           1j*rho0*c/S*sin(k*L); ...
%     1j*S/(rho0*c)*sin(k*L),      cos(k*L) ...
%     ];
% 
% 
%     T_total = T * T_total;
%         T11 = T_total(1, 1);  T12 = T_total(1, 2);
% 
%     A( iLength ) = -10*log10( abs( T11 + T12 / Z )^2 );
% 
% end
% 
% 
% figure( ); ...
%     plot( test_lengths * 1e3, A );  grid on;
%     xlabel( 'Total Recorder Length [mm]' );  ylabel( 'Amplitude [dB]' );
%     title( 'Amplification Versus Recorder Length' );



%% Part b

% L_net = 0.1626;  % Meters
% 
% a = 0.009 / 2;  % Meters
%     L_o = 0.61*a;  % Slide 18 of Lecture 2 slide set.
% 
% f = 698;  % Hz
%     k = 2*pi*f/c;  % The wave number for the respective frequency.    
% 
% S = pi/4*(0.009)^2;  % squared-meters
% 
% 
% test_lengths = 0:0.00005:0.16;
%     test_lengths = L_net - test_lengths;
%         min( test_lengths )
% 
% 
% nLengths = length( test_lengths );
%     A = zeros( nLengths, 1 );
% 
% 
% for iLength = 1:1:nLengths
% 
%     L = test_lengths(iLength);
%         L_duct_2 = L;
%         L_duct_1 = L_net - L_duct_2;
%     %
%     % Assume duct lengths are not offset by radius of hole\orifice side branch.
% 
%     % fprintf( 1, '', );
% 
% 
%     T_total = [ 1 0; 0 1 ];
% 
% 
%     % End duct.
%     T_1 = [ ...
%     cos(k*L_duct_1),                           1j*rho0*c/S*sin(k*L_duct_1); ...
%     1j*S/(rho0*c)*sin(k*L_duct_1),      cos(k*L_duct_1) ...
%     ];
% 
% 
% 
%     % Orifice side branch.
%     epsilon = 0.006 / 0.009;  % 0.67
%         a = 0.006 / 2;
% 
%     L_o = a * ( 0.9326 - 0.6196*epsilon );  % Lecture 3, Slide 11
%         L_e = 0.004 + 2*L_o;
%     %
%     Z_A = 1j * rho0 * 2 * pi * f * L_e / ( pi*0.006^2/4 );
%         T_Branch = [ 1  0;  1/Z_A  1 ];
%     %
%     % R_A is neglected (energy loss).
% 
% 
% 
%     % Front duct.
%     T_2 = [ ...
%     cos(k*L_duct_2),                           1j*rho0*c/S*sin(k*L_duct_2); ...
%     1j*S/(rho0*c)*sin(k*L_duct_2),      cos(k*L_duct_2) ...
%     ];
% 
% 
% 
%     % End termination.
%     a = 0.009/2;  L_o = 0.61*a;
%     L_e = L_duct_1 + L_o;  Z = 1j * rho0 * c / S * tan( k* L_e );
% 
% 
%     T_total = T_2 * T_Branch * T_1 * T_total;
%         T11 = T_total(1, 1);  T12 = T_total(1, 2);
% 
%     A( iLength ) = -10*log10( abs( T11 + T12 / Z )^2 );
% 
% end
% 
% 
% figure( ); ...
%     plot( test_lengths * 1e3, A );  grid on;
%         set( gca, 'XDir', 'reverse' );
%     xlabel( 'Offset from End of 162.6 mm Length Recorder [mm]' );  ylabel( 'Amplitude [dB]' );
%     title( 'Amplification Versus Offset from End of Recorder' );



%% Part b Verification

a = 0.009 / 2;  % Meters
    L_o = 0.61*a;  % Slide 18 of Lecture 2 slide set.

f = 698;  % Hz
    k = 2*pi*f/c;  % The wave number for the respective frequency.    

S = pi/4*(0.009)^2;  % squared-meters


f = 0:1:5e3;

nFreq = length( f );
    A = zeros( nFreq, 1 );


% L_net = 0.1626;  % Too short;  Revisit Part a.
%     L1 = 0.04725;        % 1,259 Hz;  1,650 Hz
%     L1 = 0.04;             % 1,228 Hz;  1,827 Hz 
%     L1 = 0.01;             % 1,058 Hz;  2,108 Hz
%     L1 = 0.005;           % 1,041 Hz;  2,081 Hz

L_net = 0.325;
    % L1 = 0.01;         % 529 Hz;  1,057 Hz
    L1 = 0.006;            % 525 Hz;
    % L1 = 0.005;       % 524 Hz;  1,049 Hz
    
L2 = L_net - L1;


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


    % End termination.
    a = 0.009/2;  L_o = 0.61*a;
    L_e = L1 + L_o;  Z = 1j * rho0 * c / S * tan( k* L_e );


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


