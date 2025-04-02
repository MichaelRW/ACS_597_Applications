


%% Synopsis

% Problem 2 - Simulation of Sources from Monopoles



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

set( groot, 'DefaultFigurePosition', [ 100  750    750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Anonymous Function Definitions



%% Strikeforce

c = 343;  % m\s
rho0 = 1.21;  % kg/m^3

f = 1e3;
    lambda = c / f;  % 0.343 m

fractions = [ 0.5  0.75  1 1.25  1.5  2  2.5  3  3.5  5  10 ];

xyz_sources = [ 0, 0, 0 ];
Q_sources = 1 + 1i*1;



%% Problem 2a - Pressure Versus Distance for Monopole

x = lambda .* fractions;
    xyz_receivers = [ x.'  zeros( numel(x) , 1 )  zeros( numel(x) , 1 ) ];

p_monopole = sum_of_monopoles( xyz_sources, Q_sources, xyz_receivers, f, rho0, c );

% figure( 'Name', 'Monopole - Pressure Versus Distance' ); ...
%     loglog( x, abs(p_monopole) );  grid on;
%     xlabel( 'Distance [m]' );  ylabel( 'Pressure [dB]' );

% return

%% Problem 2b - Directivity Patterns

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


% return



% Monopole
% figure( 'Name', 'Monopole - Directivity' ); ...
%     temp2 = temp * fractions( 1 );
%     polarplot( angle(temp2), abs(temp2) );
%     hold on;
%     for index = 2:1:numel( fractions )
%         temp2 = temp * fractions( index );
%         polarplot( angle(temp2), abs(temp2), 'Color', [ 0.00, 0.45, 0.74 ] );
%     end
%     grid on;


% Dipole - Equal Strength, In-phase Sources
% source_1 = [ 0, 0, 0 ];
% source_2 = [ lambda, 0, 0 ];
%     dipole_xyz_sources = [ source_1;  source_2 ];
% 
% dipole_Q_sources = repmat( Q_sources, 2, 1 );
% 
% p_monopole = sum_of_monopoles( dipole_xyz_sources, Q_sources, xyz_receivers, f, rho0, c );


    
%%

% return

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


