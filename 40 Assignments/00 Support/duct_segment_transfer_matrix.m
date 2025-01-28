

%% Synopsis

% This function calculates the transfer matrix for a rectangular duct segment.

function [ T ] = duct_segment_transfer_matrix( f, rho0, c, L, S )
%
%  f:  The working frequency (Hertz).
%  rho0:  The ratio of specific heats (1.2;  unitless).
%  c:  The speed of sound in air (meters per second).
%  L:  The length of the duct section (meters).
%  S:  The cross-sectional area of the duct (squared-meters).


k = 2*pi*f/c;  % The wave number for the respective frequency.

T = [ ...
    cos(k*L),  1j*rho0*c/S*sin(k*L); ...
    1j*S/(rho0*c)*sin(k*L),  cos(k*L) ...
    ];
%
% See slide 4 of Lecture 3 slide set.


