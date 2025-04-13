

close all;  clear;  clc;



x = 1:1:10;
    x = x(:);

rng( 0 );
    y = randn( numel(x), 3 );

figure( ); ...
    plot( x, y );  grid on;



color_maps = { 'Blues', 'Greens', 'Oranges' };

for frequency_index = 1:1:3

    OFFSET = 3;
        cmap = slanCM( color_maps{ frequency_index }, numel( 4 ) + OFFSET );
            cmap(1:OFFSET, :) = [ ];
                cmap_vectors( frequency_index, :, :, : ) = flipud( cmap );

end


figure( ); ...
    plot( x, y );
        colororder( cmap_vectors )