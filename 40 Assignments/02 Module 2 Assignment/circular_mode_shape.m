function psi = circular_mode_shape(n_r, n_theta, a, plot_shape)

psi = bessel_extrema_table(n_theta,n_r);

if plot_shape
    r = linspace(0, a, 30);
    theta = linspace(0, 2*pi, 100);
    
    [R, TH] = meshgrid(r,theta);
    
    P = besselj(n_theta,psi*pi*R/a).*cos(TH*n_theta);
    Y = R.*cos(TH);
    Z = R.*sin(TH);
    
    figure
    pcolor(Y,Z,P)
    map = rwb_colormap;
    colormap(map)
    
    p_max = max(abs(P),[],'all');
    clim([-1 1]*p_max)
    shading interp
    axis equal
    xlabel('x')
    ylabel('y')
    colorbar
end

end

function psi_mn = bessel_extrema_table(m,n)

T = [0	1.2197	2.2331	3.2383	4.2411;
0.5861	1.697	2.714	3.7261	4.7312;
0.9722	2.1346	3.1734	4.1923	5.2036;
1.3373	2.5513	3.6115	4.6428	5.6624;
1.6926	2.9547	4.0368	5.0815	6.1103;
2.0421	3.3486	4.4523	5.5108	6.5494;
2.3877	3.7353	4.86	5.9325	6.9811;
2.7034	4.1165	5.2615	6.3477	7.4065];

psi_mn = T(m+1,n+1);

end

function map = rwb_colormap

nSteps = 50;
b = [ones(nSteps,1); linspace(1,0,nSteps)'];
r = [linspace(0,1,nSteps)'; ones(nSteps,1)];
g = [linspace(0,1,nSteps)'; linspace(1,0,nSteps)'];

map = [r g b];
end
