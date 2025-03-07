


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



%% Calculate Admittances

masses = [ 1  50 ];
stiffnesses = [ 0  1800 100e3 ];
dampings = [ 0  0  0 ];

frequencies = 0:0.01:40;

[ FRF ] = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );  % 4001-by-2-by-2
%
    % 'admittance'
    % 'mobility'
    % 'accelerance'
    % 'stiffness'
    % 'impedance'
    % 'mass'




%% Clean-up

fprintf( 1, '\n\n\n*** Processing Complete ***\n\n\n' );



%% Reference(s)








































