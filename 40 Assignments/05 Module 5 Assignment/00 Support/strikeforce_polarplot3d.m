


%% Environment

close all;  clear;  clc;
% restoredefaultpath

set( 0, 'DefaultFigureWindowStyle', 'docked' );



%%

[ t, r ] = meshgrid( linspace(0, 2*pi, 361), linspace( -4, 4, 101) );
    [ x, y ] = pol2cart( t, r );

P = peaks(x,y);


anglesDisplayed = [ 30 270 ]*pi/180;  % Angles displayed.
radiiDisplayed = [ 4 0.8 ];


axprop = {'DataAspectRatio',[1 1 8],'View', [-12 38],...
          'Xlim', [-4.5 4.5],       'Ylim', [-4.5 4.5],...
          'XTick',[-4 -2 0 2 4],    'YTick',[-4 -2 0 2 4]};


figure( 'color', 'white' ); ...
    polarplot3d( P, 'plottype', 'meshc', 'angularrange', anglesDisplayed, 'radialrange', radiiDisplayed, 'meshscale', 2, 'polargrid',{8 8} );
    set( gca, axprop{:} );


