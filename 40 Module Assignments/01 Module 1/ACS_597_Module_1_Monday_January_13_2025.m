

%% Synopsis

% Homework Set 1 - Cut-on Frequencies in Ducts and Pipes

% Note:  Send draft of report before submission for comments.
%
% Dimensions are annontated in the class notes.



%% Environment

close all; clear; clc;
% restoredefaultpath;

set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.0 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Define Values and Functions

c = 343;  % The speech of sound in meters per second.


h_f_cut_on_rectangular_duct = @( c, L )  0.5 .* c ./ L;
%
% c - The speed of sound.
% L - The largest cross-section dimension of the rectangular duct.


h_f_cut_on_circular_duct = @( c, d )  0.568 .* c ./ d;
%
% c - The speec of sound.
% L - The diameter of the circular duct.



%% Problem 1a

% The cross-sectional dimensions for the rectangular duct are:  Lx = 12 cm and Ly = 20 cm.

% The largest dimension is Ly = 20 cm or 0.2 m.

% The cut-on frequency is,
h_f_cut_on_rectangular_duct( c, 0.2 );  % 857.5 Hz

% Fco = 858;  Hz






%% Problem 1b






% b.)  A = pir^2 = LxLy
% 
% 
% c.)  fco increases with water.
% 
% 
% d.)  
% 
% 
% Latex document.  Units not italized.  space between number and unit
% 
% 
% Focus on interpretation.
% 
% 
% Plot ideas:
%     1.)  Example:  Temperation, C;  Temperature [C] - y-axis and x-axis
%     2.)  
    


%% Clean-up

if ( ~isempty( findobj( 'Type', 'figure' ) ) )
    monitors = get( 0, 'MonitorPositions' );
        if ( size( monitors, 1 ) == 1 )
            autoArrangeFigures( 2, 2, 1 );
        elseif ( 1 < size( monitors, 1 ) )
            autoArrangeFigures( 2, 2, 1 );
        end
end

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)



%% Delete


% %% Fireworks Time Series
%
% [ s, fs ] = audioread( 'Test_shoot_A.wav' );
% %
% assert( fs == 48e3, '*** Sample rate is not 48 kHz. ***' );
% 
% sample_period_seconds = 1 / fs;
%     time_indices = ( 0:1:( size( s, 1 ) - 1 ) ) .* sample_period_seconds;
% 
% figure( ); ...
%     plot( time_indices, s( :, 1 ) );  hold on;  % Channel 1
%     plot( time_indices, s( :, 2 ) );  % Channel 2
%     plot( time_indices, s( :, 3 ) );  % Channel 3
%     plot( time_indices, s( :, 4 ) );  grid on;  % Channel 4
%     legend( 'Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Location', 'NorthWest' );
%     xlabel( 'Time [seconds]' );  ylabel( 'Amplitude [WU]' );  title( 'Fireworks Recording' );
%     xlim( [ time_indices( [ 1 end ] ) ] );  ylim( [ -0.15 +0.34 ] );
%     %
%     if ( PRINT_FIGURES == 1 ), print( fullfile( '.', 'HS9_Q1_Time_Series.png' ), '-dpng' );  end


