


%% Synopsis

% Problem 4 - Traffic Noise Modeling

% See Lecture 22 on Monday, April 7, 2025
%   ""D:\15 Downloads\00 GitHub\ACS_547\35 Lectures\22 Monday, April 7, 2025\Lecture 22 - Coherence effects - filled.pptx""



%% Note(s)

% Incoherent line source.

% Distance (source and receivers) are in units of meters.

% Source can be real-valued or complex-valued.  Complex-valued sources have magnitude and phase information.

% Sources have units of volume velocity, m^3/s.

% Complex pressures are in units of Pascals.



%% Environment

close all; clear; clc;
% restoredefaultpath;

addpath( genpath( './00 Support' ), '-begin' );

% set( groot, 'DefaultFigurePosition', [ 100  750    750  500 ] );  % x, y, width, height

set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );



%% Contants and Parameters

source_separation_distance = 8;  % m

average_A_weighted_source_pressure = 86;  % dB SPL A at 1 m

distance_to_house_from_line_source = 200;  % m



%% Part 4a

L1 = average_A_weighted_source_pressure;
    L2 = L1;
        L3 = L1;

% sound_pressure_net = 10*log10( 10^(L1/10) + 10^(L2/10 ) + 10^(L3/10 ) );
sound_pressure_net = 10*log10( 10^(L1/10) + 10^(L2/10 ) )


data = [ L1, L2, L3 ];

average_sound_level = 10*log10( mean( 10.^( data./10 ) ) )
LA_eq = 10*log10( ( 1/numel( data ) ) * sum ( 1.*10.^( data ./ 10 )) )

return

%% Part 4b



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


