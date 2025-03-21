


%% Synopsis

% This script-file exercises the nDOF_direct_solution function script-file.



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( './00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% 1 DOF Example

masses = 50;
stiffnesses = [ 100e3  0 ];
dampings = [ 0  0 ];
    % dampings = [ 0  (1 + 0.1*1j) ];

wo = sqrt( 100e3 / 50 );
    fo = wo/(2*pi);  % 7.1 Hz

frequencies = 0:0.1:100;


admittance = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );  % 4001-by-2-by-2
%
figure( ); ...
    loglog( frequencies, abs( admittance ) );  grid on;
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );


impedance_original_function = nDOF_direct_solution_org( masses, stiffnesses, dampings, frequencies, 'impedance' );  % 4001-by-2-by-2
impedance_updated_function = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'impedance' );  % 4001-by-2-by-2
%
figure( ); ...
    loglog( frequencies, abs( impedance_original_function ) );  hold on;
    loglog( frequencies, abs( impedance_updated_function ) );  grid on;
        legend( 'Original Function', 'Updated Function', 'Location', 'North' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Impedance [$\frac{N\cdots}{m}$]' );

return

%% 2 DOF Example

masses = [ 1  50 ];
stiffnesses = [ 0  1800 100e3 ];

dampings = [ 0  0  0 ];
dampings = [ 0  (1 + 0.1*1j)  (1 + 0.1*1j) ];

frequencies = 0:0.01:40;

[ FRF ] = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );  % 4001-by-2-by-2


admittance = zeros( numel( frequencies ), 1 );

for index = 1:1:numel( frequencies )
    temp = diag( squeeze( FRF( index, :, : ) ) );
        temp2 = temp(1) + temp(2)*1j;
        admittance( index ) = abs( temp2 );
end

clear temp1 temp2;


% figure( ); ...
%     semilogy( frequencies, admittance );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );


figure( 'Name', 'Admittance' ); ...
    semilogy( frequencies, abs( FRF( :, 1, 1 ) ) );  hold on;
    semilogy( frequencies, abs( FRF( :, 1, 2 ) ) );
    semilogy( frequencies, abs( FRF( :, 2, 1 ) ) );
    semilogy( frequencies, abs( FRF( :, 2, 2 ) ) );  grid on;
        legend( 'M1 Movement with F1', 'M1 Movement with F2', 'M2 Movement with F1', 'M2 Movement with F2', 'Location', 'North' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
    
figure( 'Name', 'Dynamic Stiffness' ); ...
    semilogy( frequencies, 1 ./ abs( FRF( :, 1, 1 ) ) );  hold on;
    semilogy( frequencies, 1 ./ abs( FRF( :, 1, 2 ) ) );
    semilogy( frequencies, 1 ./ abs( FRF( :, 2, 1 ) ) );
    semilogy( frequencies, 1 ./ abs( FRF( :, 2, 2 ) ) );  grid on;
        legend( 'M1 Movement with F1', 'M1 Movement with F2', 'M2 Movement with F1', 'M2 Movement with F2', 'Location', 'North' );
    xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );



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








































