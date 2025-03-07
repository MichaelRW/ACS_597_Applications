


%% Synopsis

% From the recorded video on Dynamic Vibration Absorber (DVA).



%% Note(s)

% Dynamic vibration absorption is useful if:
%
%   a.)  Can not get away from the resonance frequency;  can not tune the
%   system to the resonance is away from the frequency of the forcing
%   function.
%
%   b.)  



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( './00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Example 1

m =  50;  % kg
k = 100e3;  % N\m

wo = sqrt( k / m);  % 44.7 radians\s
    fo = wo/(2*pi);  % 7.1 Hz


% Add a Dynamic Vibration Absorber (DVA)
m_dva = 1;  % kg
k_dva = 1800;  % N\m

wo_dva = sqrt( k_dva / m_dva );  % 42.4 radians\s
    fo_dva = wo_dva/(2*pi);  % 6.8 Hz


% Solve
m1 = m_dva;
k1 = 0;  % Spring above DVA mass;  does not exist;  set to zero.

m2 = m;
k2 = k;  % 100e3 N/m
    k2 = k*(1+0.1*1j);

k3 = k_dva;  % 1,800 N/m
    k3 = k_dva*(1+0.1*1j);



%% Blocked Frequencies (Euqation 9.1 of the Notes)

w1 = sqrt( ( k1 + k3 ) / m1 );  % 42.4 radians\s;  smaller than the system resonant frequency.
    f1 = w1/(2*pi);  % 6.8 Hz


w2 = sqrt( ( k2 + k3 ) / m2 );  % 45.1 radians\s;  greater than the system resonant frequency.
    f2 = w2/(2*pi);  % 7.2 Hz



%% Coupled Frequencies (Equation 9.15 of the Notes)

mu_4 = ( k3^2 / (m1*m2) );  % 64.8e3 (radians\s)^4
    mu = ( mu_4 )^0.25;

w(1) = 0.5*( ( w1^2 + w2^2 )  +  sqrt( ( w1^2 - w2^2 )^2 + 4*mu_4 ) );
w(2) = 0.5*( ( w1^2 + w2^2 )  -  sqrt( ( w1^2 - w2^2 )^2 + 4*mu_4 ) );
    w = sqrt( w );


w_minus = w(2);              % 40.5 radians\s  (lower than 42.4 and 45.1 radians\s)
    f_minus = w(2)/(2*pi);  % 6.4 Hz
%
w_plus = w(1);                 % 46.9 radians\s (higher than 42.4 and 45.1 radians\s)
    f_plus = w(1)/(2*pi);    % 7.4 Hz
%
% Note(s):
%
%   The blocked frequencies are always between these two frequencies.


clear w mu_4;



%% Mode Shapes (Equations 9.18 and 9.19 of the Notes)

phi_plus = [ ...
    1; ...
    -(1/mu^2)*sqrt(m1/m2)*(w_plus^2 - w1^2), ...
    ];

phi_minus = [ ...
    1; ...
    (1/mu^2)*sqrt(m1/m2)*(w1^2 - w_minus^2), ...
    ];



%% Plot Admittance - Magnitude of Displacement over Force

F0 = 1;

w = 0:0.01:40*2*pi;
    f = w/(2*pi);


m = [ 1 50 ];
k = [  0  100e3*(1+0.1*1j)  1.8e3*(1+0.1*1j)  ];

X2 = F0/(m(1)*m(2))*k(3)./((w.^2-w_plus^2).*(w.^2-w_minus^2));
X2_original = F0./m(2)*1./(w2^2-w.^2);


figure( ); ...
    semilogy( f, abs(X2_original) );  hold on;
    loglog( f, abs(X2) );  grid on;
        legend( 'Original', 'With DVA' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );




%% ph


% m = [ 1 50 ];
% k = [ 0 100e3 1.8e3 ];
% 
% FRF = nDOF_direct_solution( m, k, [ ], 1:1:1e3, 'admittance' )


%       masses:  An nDOF vector of masses.
%       stiffnesses:  An nDOF+1 vector of stiffnesses, with the first and last entries being
%           connections to ground.
%       dampings:  An nDOF+1 vector of stiffnesses, with the first and last entries being
%       connections to ground.
%       freq:  An nFreq vector of frequencies.
%       FRF_type:  Is one of:  'admittance', 'mobility', 'accelerance', 'stiffness', 'impedance', or 'mass'.



%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


