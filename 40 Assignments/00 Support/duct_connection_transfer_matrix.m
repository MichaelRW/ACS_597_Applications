

%% Synopsis

% This function calculates the transfer matrix for a rectangular duct segment.

function [ T ] = duct_connection_transfer_matrix( S_A, S_B )
%
%  S_A:  Area on the right.
%  S_B:  Area on the left.


T = [ ...
    1  0; ...
    0  S_A / S_B; ...
    ];
%
% See slide 9 of Lecture 2 slide set.


