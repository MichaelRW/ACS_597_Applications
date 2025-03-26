


%% Synopsis

% Problem 3 - Washing Machine Dynamic Vibration Absorber (DVA) Design



%% Note(s)

% At 300 RPM, the load frequency, this system will isolated vibration more
% than the 2 DOF and 1 DOF systems.

% The DVA behaves like a Helmholtz resonator, working as a mechnical notch
% filter at a particular frequency.

% DVAs work well in systems that have one dominate mode.  A good
% example of this is building sway and its compensation.  Based on the
% design and construction of a building, it will typically have a single,
% dominate mode.

% DVA is useful if you can not get away from the resonance frequency;  can
% not tune a 2 DOF system when the resonance frequency is away from the
% forcing frequency.



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



%% Dynamic Vibration Absorber (DVA)

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


if ( 0 )

    m = [ 1 50 ];
    k = [  0  100e3  1.8e3  ];
    % k = [  0  100e3*(1+0.1*1j)  1.8e3*(1+0.1*1j)  ];
    dampings = [ 0 0 0 ];

    X2 = F0/(m(1)*m(2))*k(3)./((w.^2-w_plus^2).*(w.^2-w_minus^2));
    X2_original = F0./m(2)*1./(w2^2-w.^2);

else

    m = [ 50 1 ];
    k = [ 100e3 1.8e3 0 ];
    % dampings = [ 0 0 0 ];

    X2 = F0/(m(2)*m(1))*k(1)./((w.^2-w_plus^2).*(w.^2-w_minus^2));
    X2_original = F0./m(1)*1./(w2^2-w.^2);

end


figure( 'Name', 'Displacement with DVA' ); ...
    loglog( f, abs(X2_original) );  hold on;
    loglog( f, abs(X2) );  grid on;
        legend( 'Original Displacment', 'Displacement With DVA', 'Location', 'NorthWest' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Displacment [m]' );




%% Admittance

m = [ 50 1 ];
k = [ 100e3 1.8e3 0 ];
dampings = [ 0 0 0 ];
%
% Damping can be added by using a dampings vector or appending a complex
% value to the respective spring stiffness.  THESE ARE NOT INTERCHANGABLE.


f = 0:0.1:100;  % No matrix singularity warning issued (set 0 to 0.1 if it occurs).

FRF = nDOF_direct_solution( m, k, dampings, f, 'admittance' );  % Symmetric matrix.
%
squeeze( FRF(100, :, :) );
%
%   -1.1099e-05   9.6547e-06
%   9.6547e-06  -0.00049166
%
% (1, 1) - Force on M1 and its associated displacement.
% (2, 2) - Force on M2 and its associated displacement.
%
% (1, 2) - Force on M1 and the displacement of M2.
% (2, 1) - Force on M2 and the displacment on M1.
%
%   These are interchangeable;  the same result.


FRF_org = nDOF_direct_solution_org( m, k, dampings, f, 'admittance' );  % Symmetric matrix.
%
squeeze( FRF_org(100, :, :) );
%
%   -1.1099e-05   9.6547e-06
%   9.6547e-06  -0.00049166


figure( 'Name', 'Admittance of DVA' ); ...
    loglog( f, abs( FRF( :, 1, 1 ) ) );  hold on;
    loglog( f, abs( FRF_org( :, 1, 1 ) ), 'LineStyle', '--' );  grid on;
        legend( 'Updated Function', 'Original Function', 'Location', 'NorthWest' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );



%% Clean-up

if ( ~isempty( findobj( 'Type', 'figure' ) ) )
    monitors = get( 0, 'MonitorPositions' );
        if ( size( monitors, 1 ) == 1 )
            autoArrangeFigures( 3, 4, 1 );
        elseif ( 1 < size( monitors, 1 ) )
            autoArrangeFigures( 2, 2, 2 );
        end
end

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)


