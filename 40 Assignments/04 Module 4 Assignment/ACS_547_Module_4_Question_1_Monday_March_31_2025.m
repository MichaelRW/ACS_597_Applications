


%% Synopsis

% Problem 1 - Simulation of Sources from Monopoles



%% To Do

% Confirm relationship of pressure dependence on distance for a given frequency.

% 



%% Note(s)

% Distance (source and receivers) are in units of meters.

% Source can be real valued or complex valued.  Complex valued sources have magnitude and phase information.

% Complex pressures are in units of Pascals.



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

set( groot, 'DefaultFigurePosition', [ 100  450    750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Parameters

c = 343;  % m\s
rho0 = 1.21;  % kg/m^3

f = [ 1e2  1e3  1e4 ].';  f_set = { '100 Hz', '1 kHz', '10 kHz' };
    lambda = c ./ f;  % 3.43 m, 0.343 m, and 0.0343 m

k = (2*pi*f) ./ c;      

radius_fractions = 0.1:1:4;



%% Problem 2a - Pressure Versus Distance for Monopole

% Monopole distance dependence.  Decay is 6 dB per double of distance.

xyz_sources = [ ...
    0, 0, 0; ...
    ];
%
Q_sources = 1;


x = 0.01:0.01:10;
    y = zeros( length( x ), 1 );  z = zeros( length( x ), 1 );
        xyz_receivers = [ x(:) y(:) z(:) ];


for index = 1:1:numel( f )
    p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, f( index ), rho0, c );
        L( index, : ) = 10*log10( abs(p).^2 );
end


% figure( 'Name', 'Monople - Pressure Magnitude Versus Distance' ); ...
%     semilogx( x, L );  grid on;
%         legend( f_set, 'Location', 'NorthEast' );
%     xlabel( 'Distance [m]' );  ylabel( 'Pressure Magnitude [$dB$]' );

% return

%% Problem 2b - Monopole Directivity Pattern

xyz_sources = [ ... 
    0, 0, 0; ...
    ];
%
Q_sources = 1;


r = radius_fractions;  % m

theta = 0:0.01:2*pi;

xyz_receivers = nan( numel( r ), numel( theta), 3 );

for radius_index = 1:1:numel( r )

    x = r( radius_index )*cos( theta );  y = r( radius_index )*sin( theta );  z = zeros( size( x ) );
        xyz_receivers( radius_index, :, :, : ) = [ x(:) y(:) z(:) ];

end


% color_maps = { 'Gray', 'Bone', 'Copper' };
color_maps = { 'Blues', 'Greens', 'Oranges' };
        



    for frequency_index = 1:1:numel( f )
    % for frequency_index = 1

        h = figure( 'Name', 'Monople - Spherical Spreading' ); ...

        OFFSET = 3;
            cmap = slanCM( color_maps{ frequency_index }, numel( r ) + OFFSET );
                cmap(1:OFFSET, :) = [];

        tick_marks = [ ];

        for radius_index = 1:1:numel( r )
    
            p = sum_of_monopoles( xyz_sources, Q_sources, squeeze( xyz_receivers( radius_index, :, : ) ), f( frequency_index ), rho0, c );
                L = 10*log10( abs(p).^2 );
                    polarplot( theta, L, 'Color', cmap( radius_index, : ) ); hold on;
                        tick_marks = [ L(1) tick_marks ];

        end


        colormap( gca( h ), cmap );
            h_colorbar = colorbar( );
                h_colorbar.Ticks = linspace( 0, 1, 4 );
                h_colorbar.TickLabels = num2cell( round( tick_marks ) );

                keyboard
            %
            h_colorbar.Label.String = 'Sound Pressure [dB]';

        rlim( [ 0, 100 ] );

    end

return
                    





figure
polarplot(a, r1)
hold on
polarplot(a, r2)
hold off

    h1 = subplot( 1, 2, 1 ); ...
        plot( nan, nan );  hold on;
        xlabel( 'Angle [radians]' );  ylabel( 'Pressure Magnitude [$dB$]' );
    subplot( 1, 2, 2 ); ...
        h2 = polarplot( nan, nan );  hold on;


for index = 1:1:numel( r )

    x = r( index )*cos( theta );
    y = r( index )*sin( theta );
    z = zeros( size( x ) );
        xyz_receivers = [ x(:) y(:) z(:) ];

    p = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, f(1), rho0, c );
    p2 = abs(p).^2;
        L = 10*log10( p2 );

    subplot( h1 ); ...
        plot( theta, L.' );  grid on;

    % polarplot( h2 ); ...
        polarplot( theta, L );  grid on;

end

return

%% ph

return

N = 128;
    n = 0:1:N;
        temp = exp( 1i.*(2*pi)/N*n);

% Set of unit magnitude points around the origin.
x = real( temp ).';  y = imag( temp ).';
    xy_unit_set = [ x, y ];
        xyz_unit_set = [ xy_unit_set  zeros( size( xy_unit_set, 1 ), 1 ) ];

% Scale by fractions.


source_1 = [ 0, 0, 0 ];
source_2 = [ 5*lambda, 0, 0 ];
% source_2 = [ 0.6*lambda, 0, 0 ];
% source_2 = [ 0.1*lambda, 0, 0 ];
    dipole_xyz_sources = [ source_1;  source_2 ];


k = (2*pi) / lambda;
    d = 1 / k

    aValue = 6e-1;

    dipole_xyz_sources = [ ...
        aValue  aValue  0; ...
        aValue  -aValue  0; ...
        -aValue  aValue  0; ...
        -aValue  -aValue  0 ];  % Lateral Quadrupole

% dipole_Q_sources = repmat( Q_sources, 2, 1 );
% dipole_Q_sources = [ Q_sources;  0.5.*Q_sources ];

     dipole_Q_sources = [ ...
         1 + 1i; ...
         1 - 1i; ...
         1 + 1i; ...
         1 - 1i ];


figure; ...
    polarplot( nan, nan );  hold on;


for fraction_index = 1:1:numel( fractions )

    for xy_index = 1:1:size( xy_unit_set, 1 )

        % p_monopole = sum_of_monopoles( xyz_sources, Q_sources, xyz_unit_set.*fractions( fraction_index ), f, rho0, c );
        % p_monopole = sum_of_monopoles( dipole_xyz_sources, dipole_Q_sources, xyz_unit_set.*fractions( fraction_index ), f, rho0, c );
        p_monopole = sum_of_monopoles( dipole_xyz_sources, dipole_Q_sources, xyz_unit_set, f, rho0, c );
            temp = xy_unit_set(:, 1) + 1i*xy_unit_set(:, 2 );
                polarplot( angle(temp), abs(p_monopole)./1e3 );

        % keyboard

    end

    % keyboard;

end

figure


return






% Dipole - Equal Strength, In-phase Sources
% source_1 = [ 0, 0, 0 ];
% source_2 = [ lambda, 0, 0 ];
%     dipole_xyz_sources = [ source_1;  source_2 ];
% 
% dipole_Q_sources = repmat( Q_sources, 2, 1 );
% 
% p_monopole = sum_of_monopoles( dipole_xyz_sources, Q_sources, xyz_receivers, f, rho0, c );


    
%% Clean-up

if ( ~isempty( findobj( 'Type', 'figure' ) ) )
    monitors = get( 0, 'MonitorPositions' );
        if ( size( monitors, 1 ) == 1 )
            autoArrangeFigures( 3, 4, 1 );
        elseif ( 1 < size( monitors, 1 ) )
            autoArrangeFigures( 2, 2, 2 );
        end
end

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)

p = 10 + 1i*5;  % Pascals;  sinusoid

p_mag = abs( p )  % 11.18 Pascals

p_rms = p_mag / sqrt(2)  % 7.91 Pascals RMS

p_dB_SPL = 20*log10( p_rms / 20e-6 )  % 111.94 dB SPL Z


p_dB_SPL_verify = convert_complex_pressure_to_dB_SPL( p )


