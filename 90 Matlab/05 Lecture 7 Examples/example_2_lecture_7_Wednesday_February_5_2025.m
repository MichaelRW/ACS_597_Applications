


%% Synopsis

% ACS 547, Lecture 7 - Example 2, Sabine Rooms

% See slide 15 on "Lecture 07 - Sabine rooms - Filled.pptx".



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 1.5 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Information

room.width = 10;  room.length = 15;  room.height = 6;  % meters

target_dB_level = 85;  % dB

% A room is 10 m × 15 m × 6 m high.  The room currently contains several machines 
% and workers.  Four new machines are to be installed.  You are tasked with making sure 
% the new machines do not bring the overall sound levels above 85 dB.

% The machines are available in two models:
%
%   “Loud”, which produce sound power of 94 dB re: 1 pW
%   “Quiet”, which produce sound power of 84 dB re: 1 pW, but cost $1,600 more each


% What measurements must be done?

% What treatment options are available if the levels are too high?

% What is the most cost-effective solution?



%% Step 1 - Baseline Measurements

% Measure the existing sound pressure levels with the machines TURNED ON.
Lp = 75;  % dB re: 20e-6 Pascals (references varies with person and location).


% Measure the reverberation time of the room with the machines TURNED OFF.
T60 = 2.4;  % seconds


% Calculate the volume and area of the room.
room_volume = room.width * room.length * room.height;  % 900 m^3
room_area = 2*(room.width * room.height)  +  2*(room.length * room.height)  +  2*(room.width * room.length );  % 600 m^2


%% Step 2 - Calculate Baseline Conditions

h_average_aborption_coefficient = @( volume, area, c, T60 )  ( 55.25 .* volume) ./ ( area .* c .* T60 );

h_room_constant = @( area, average_absorption_coefficient )  ( area .* average_absorption_coefficient ) ./ ( 1 - average_absorption_coefficient  );


alpha_average = h_average_aborption_coefficient( room_volume, room_area, 343, T60 );  % 0.1 unitless

room_constant = h_room_constant( room_area, alpha_average );  % 66.2 m^2 or Sabines



%% Step 3 - Initial Predicted Level

Lp_original = 75;  % dB re: 20e-6 Pascals - Level with the existing machines.


% Calculate the sound pressure level for one new machine.
Lw_quiet = 84;  % dB
    Lp_new_quiet_one_machine = Lw_quiet + 10*log10( 4 / room_constant );  % 71.7 dB

Lw_loud = 94;  % dB
    Lp_new_loud_one_machine = Lw_loud + 10*log10( 4 / room_constant );  % 81.7 dB


% Calculate the new total sound pressure levels for both types of machines.
Lp_combined_quiet = 10*log10( 10^(75/10)  +  4*10^(Lp_new_quiet_one_machine/10) );  % 79.6 dB

Lp_combined_load = 10*log10( 10^(75/10)  +  4*10^(Lp_new_loud_one_machine/10) );  % 88.0 dB



%% Treatment Options

% Option:  Bought loud machines and spent money on absorbing material?

room_constant_old = room_constant;

% Target level is 85 dB.  With the 4 loud machines the new level is 88 dB.

% A pressure level difference of -3 dB is required to meet the target level.

room_constant_new = room_constant_old / 10^( -3 / 10 );  % 134 m^2 or Sabines

% If only the celing was covering in tile.
alpha_new = ( ( room.width * room.length ) * 0.6  +  ( 2*(room.width * room.height) + 2*(room.length * room.height) + room.width*room.length ) * 0.1 ) / room_area;  % 0.225 unitless

% The new room constant with the tiling added to the celing.
room_constant_new_with_treatment = room_area * alpha_new / ( 1 - alpha_new );  % 174.2 m^2 or Sabines
%
% The ceiling treatment (174.2 m^2) exceeds the required level (134 m^2),
% which should be more than adequate.



%% Cost Analysis

loud_machine_cost = 0;

% 4 quiet machines, which will meet the target pressure level of 85 dB.
additional_cost_quiet_machines = 4 * ( 1600 + loud_machine_cost )


% 4 loud machines with ceiling treatment that will meet the target pressure level of 85 dB.
additional_cost_loud_machines = 4 * loud_machine_cost + 450.0;


cost_difference = additional_cost_loud_machines - additional_cost_quiet_machines;  % -$5,950.00
%
% Buying the loud machines and use a ceiling treatment saves $5,950.00.



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

%% Overall, A-weighted Level

% Note(s):
%
%   The above analysis was done using the unweighted sound level for the 500 Hz octave band.
%
%   If an overall level is to be calculated (i.e., across a set of octave bands), then this analysis
%   must be done for all octave band center frequencies.  Once the unweighted sound pressure
%   levels at the location are determned, the respective octave band A-weighting offsets are
%   applied.
%
%   The overall A-weighted sound pressure levels is then calculated logarithmically add using the
%   expression,
%
%   10*log10( sum( 10^(Lp_a/10) )


