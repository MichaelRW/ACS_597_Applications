
function [ h, h2 ] = labelColorbar( laybel )

% Trap on.
if ( ~ischar( laybel ) ),  error( '*** Input must be a string. *** ');  end

h = colorbar( );  % Returns the colorbar object.
    h2 = get( h, 'YLabel' );  % Returns the text object belonging to the y-label.

h2.String = laybel;  % Set the label string equal to the input.
h2.Interpreter = 'Latex';

if ( nargout == 0 ),  clear;  end


