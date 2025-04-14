

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



% https://www.mathworks.com/matlabcentral/answers/507355-colormap-for-multline-plot

% https://www.mathworks.com/matlabcentral/answers/405484-how-to-get-handles-of-the-line-object

% https://www.mathworks.com/matlabcentral/answers/2091796-saving-multiple-plots-into-one-handle

% https://www.mathworks.com/matlabcentral/answers/296360-how-can-i-change-the-line-color-e-g-color-g-doesn-t-work

% https://www.mathworks.com/matlabcentral/answers/438317-how-to-color-individual-lines-in-a-plot


