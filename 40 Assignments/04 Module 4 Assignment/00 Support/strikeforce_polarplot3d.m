


%% Environment

close all;  clear;  clc;
% restoredefaultpath



%%

[t,r] = meshgrid(linspace(0,2*pi,361),linspace(-4,4,101));

[x,y] = pol2cart(t,r);

P = peaks(x,y);


t2 = [30 270]*pi/180;

r2 = [.8 4];
    r3 = fliplr(r2);


axprop = {'DataAspectRatio',[1 1 8],'View', [-12 38],...
          'Xlim', [-4.5 4.5],       'Ylim', [-4.5 4.5],...
          'XTick',[-4 -2 0 2 4],    'YTick',[-4 -2 0 2 4]};



%% Mesh plot with contours, overlay 8 by 8 polar grid

figure( 'color', 'white' ); ...
    polarplot3d( P, 'plottype', 'meshc', 'angularrange', t2, 'radialrange', r3, 'meshscale', 2, 'polargrid',{8 8} );
    set( gca, axprop{:} );


