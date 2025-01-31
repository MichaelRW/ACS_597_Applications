

%% Synopsis

% This function calculates the transfer matrix for a rectangular duct segment.

% See slide 4 of Lecture 3 slide set.


function [ T ] = duct_segment_transfer_matrix_flow( f, rho0, c, L, S, mach_number )
%
%  f:  The working frequency (Hertz).
%  rho0:  The ratio of specific heats (1.2;  unitless).
%  c:  The speed of sound in air (meters per second).
%  L:  The length of the duct section (meters).
%  S:  The cross-sectional area of the duct (squared-meters).
%  mach_number:  The Mach number for the respective duct segment.


k = 2*pi*f/c;  % The wave number for the respective frequency.

k_c = k / ( 1 - mach_number );


phase_change = exp( -1j*mach_number*k_c*L );


T = [ ...
    cos( k_c * L ),  1j*rho0*c/S*sin( k_c * L ); ...
    1j*S/(rho0*c)*sin( k_c * L ),  cos( k_c * L ) ...
    ];


T = phase_change * T;


