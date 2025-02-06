

close all; clear; clc;

pause( 1 );


p = [ ...
    0 5 6; ...
    1 9 .1; ...
    0.2 10, 0.01; ...
    ];

q = [ ...
    -0.01 6 7; ...
    9e-3 10 20; ...
    0 0 0; ...
    ];

p = cat( 3, p, q )

number_of_values = 8;

[ sortedValues, sortedIndices ] = sort( p(:) );
    smallestValues = sortedValues(1:number_of_values);
    smallestIndices = sortedIndices(1:number_of_values);


[ x, y, z ] = ind2sub( size(p), smallestIndices );
    [ x y z ]

smallestValues


