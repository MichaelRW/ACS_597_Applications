

close all;

figure( ); ...

cvars = round( abs( p ) );

cmap = flipud( jet( 11 ) );
    cvars_map = ( cvars - min(cvars) ) / ( max(cvars) - min(cvars) ) * (n-1) + 1;


temp2 = temp * fractions( 1 );
polarplot( angle(temp2), abs(temp2), 'Color', cmap(1, :) );
hold on;
for index = 2:1:numel( fractions )
    temp2 = temp * fractions( index );
    polarplot( angle(temp2), abs(temp2), 'Color', cmap(index, :) );
end
grid on;

colormap( cmap );

h = colorbar;
    set( h, 'YDir', 'reverse', 'TickLabels', [ linspace( 1, 0, 8 ) ] );

    clim([1.2 2.6])


% colorbar;

return

%%

close all;  clc;

keyboard;

% Define our data.
theta = linspace(0,(2*pi),30);
x=cos(theta);
y=sin(theta);
energy = linspace(0,1000,30);
numPoints = length(x);
% Get a colormap, a unique color for every energy level
cmap = jet(numPoints); % Initialize jet colormap.
% Get energy in the range 1 to numPoints so we can use that to get a row from the colormap.
qEnergy = imquantize(energy, numPoints);
for k = 1 : numPoints
	% Get the color for this energy level:
	thisEnergy = qEnergy(k);
	thisColor = cmap(thisEnergy);
	fprintf('Plotting point #%d at (%.3f, %.3f) with color (%.3f, %.3f, %.3f)\n',...
		k, x(k), y(k), cmap(k, 1), cmap(k, 2), cmap(k, 3));
	plot(x(k), y(k), '.', 'Color', cmap(k, :), 'MarkerSize', 40);
	hold on;
end
grid on;
axis square;
fprintf('Done running %s.m ...\n', mfilename);

return

%%

if ( 0 )

    % [ t, r ] = meshgrid( linspace(0, 2*pi, 361), linspace( -4, 4, 101) );
    [ t, r ] = meshgrid( linspace(0, 2*pi, 10), linspace( -4, 4, 10) );
        size( t ), size( r )
        [ x, y ] = pol2cart( t, r );
    %
    P = peaks(x,y);

else

    [ t, r ] = meshgrid( (2*pi)/N*n, fractions );  % t:  2-by-6;  r: 2-by-6
        [ x, y] = pol2cart( t, r );  % x:  2-by-6;  y: 2-by-6
    %
    temp = abs( p );  % 2-by-1

    P = [ repmat( temp(1), 1, numel( fractions ) );  repmat( temp(2), 1, numel( temp2 ) ) ]

end

% return

anglesDisplayed = [ 30 270 ]*pi/180;  % Angles displayed.
radiiDisplayed = [ 4 0.8 ];


axprop = {'DataAspectRatio',[1 1 8],'View', [-12 38],...
          'Xlim', [-4.5 4.5],       'Ylim', [-4.5 4.5],...
          'XTick',[-4 -2 0 2 4],    'YTick',[-4 -2 0 2 4]};

figure( 'color', 'white' ); ...
    polarplot3d( P, 'plottype', 'meshc', 'angularrange', anglesDisplayed, 'radialrange', radiiDisplayed, 'meshscale', 2, 'polargrid',{8 8} );
    set( gca, axprop{:} );