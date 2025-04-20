
function [ ] = plot_directivity( aName, )


figure( 'Name', aName ); ...
    temp2 = temp * fractions( 1 );
    polarplot( angle(temp2), abs(temp2) );
    hold on;
    for index = 2:1:numel( fractions )
        temp2 = temp * fractions( index );
        polarplot( angle(temp2), abs(temp2), 'Color', [ 0.00, 0.45, 0.74 ] );
    end
    grid on;


    