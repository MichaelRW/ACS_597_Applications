

%% Synopsis

% This function calculates the complex impedence for a circular duct segment.

% ONLY FOR CIRCULAR CROSS SECTIONS.

function [ Z ] = open_end_impedance( f, rho0, c, L, S, flanged )
%
%   f:  The working frequency (Hertz).
%   rho0:  The ratio of specific heats (1.2;  unitless).
%   c:  The speed of sound in air (meters per second).
%   L:  The length of the duct section (meters).
%   S:  The cross-sectional area of the duct (squared-meters).
%   flanged:  The end state of the circular duct (boolean).


a = sqrt( S / pi );  % Diameter of circular section.

k = 2*pi*f/c;  % The wave number for the respective frequency.

if ( flanged )
    L_0 = 8*a/(3*pi);  % Slide 18 of Lecture 2 slide set.
else
    L_0 = 0.61*a;  % Slide 18 of Lecture 2 slide set.
end

L_e = L + L_0;
%
%  L_e:  Effective length.
%  L:  Length of the segment.
%  L_0:  Offset length.

Z = 1j*rho0*c/S*tan(k*L_e);  % End impedance.


