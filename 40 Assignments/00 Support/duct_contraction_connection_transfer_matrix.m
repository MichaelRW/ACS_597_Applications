

%% Synopsis

% This function calculates the transfer matrix for a rectangular duct segment.

% See slide 4 of Lecture 3 slide set.


function [ T ] = duct_contraction_connection_transfer_matrix( rho0, c, area_upstream, area_downstream, mach_number )


K_d = 0.5 * ( 1 - area_downstream/area_upstream );


T = [ ...
    1  rho0*c/area_downstream*K_d*mach_number; ...
    0  1 ...
    ];


