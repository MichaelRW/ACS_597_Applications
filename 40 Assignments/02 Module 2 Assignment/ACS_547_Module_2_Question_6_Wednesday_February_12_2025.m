


%% Synopsis

% Lecture 11, Wednesday, February 19, 2025

% The compressor elevated above the ground.



%% Environment

close all; clear; clc;
% restoredefaultpath;

% addpath( genpath( '' ), '-begin' );
addpath( genpath( '../40 Assignments/00 Support' ), '-begin' );

% set( 0, 'DefaultFigurePosition', [  400  400  900  400  ] );  % [ left bottom width height ]
set( 0, 'DefaultFigurePaperPositionMode', 'manual' );
set( 0, 'DefaultFigureWindowStyle', 'normal' );
set( 0, 'DefaultLineLineWidth', 0.8 );
set( 0, 'DefaultTextInterpreter', 'Latex' );

format ShortG;

pause( 1 );

PRINT_FIGURES = 0;



%% Define Compressor

compressor.width = 1;  % m
compressor.depth = 1;  % m
compressor.height = 2;  % m
    compressor.area = 2*(compressor.width * compressor.depth) + 2*(compressor.width * compressor.height) + 2*(compressor.depth * compressor.height);  % 3 m^2
    compressor.volume = compressor.width * compressor.depth * compressor.height;  % m^3

compressor.power_level = 105;  % dB re: 1e-12 Watts
compressor.frequency = 50;  % Hz

c = 343;  % m/s

rho0 = 1.2;  % kg/m^3  CHECK



%% Sound Level Target

sound_level_target = 82;  % dB re: 20e-6 Pascals



%% Define Workshop

R = 40;  % m^2 or Sabins



%% Define Close-fitting Enclosure

helmholtz_factor = (2 * pi * compressor.frequency) / c;  % 0.92 m
%
%  For a small enclosure, k*d << 1.  Therefore d << 1.1.
d = 0.75;
    % helmholtz_factor * d;  % 0.69

enclosure.width = compressor.width + d;  % m
enclosure.depth = compressor.depth + d;  % m
enclosure.height = 3;  % m
    enclosure.area = 2*(enclosure.width * enclosure.depth) + 2*(enclosure.width * enclosure.height) + 2*(enclosure.depth * enclosure.height);
    enclosure.volume = enclosure.width * enclosure.depth * enclosure.height;

enclosure.E = 3.6e12;  % Pascals
enclosure.thickness = 3.81e-2;  % m
enclosure.density = 800;  % kg/m^3
enclossure.poisson_ratio = 0.25;  % Unitless

% Clamped boundary conditions.



%% Define Air Intage

air_intake_radius = 10e-2;  % m



%% Insertion Loss

% For the insertion loss to be high, we need:
%
%   1.)  Compliance of the air to be high;  volume of enclosure must be large.
%   2.)  Compliance of each enclosure wall to be low.


% The correction factor for clamped walls.  See Figure 12.4 on slide 9 of the Lecture 11 notes.
aspect_ratio = enclosure.height / enclosure.width;  % 1.7
    correction_factor = 2;  % Approximate value read from the Figure 12.4.

 h_wall_compliance = @( wall_area, correction_factor, bending_stiffness )


Ca = enclosure.volume / ( rho0 * c^2 )

return

%%

Lp_10_meters = Lw  +  10*log10( machine.D /( 4 * pi * machine.distance^2 ) );  % dB re: 20e-6 Pa
    % [ octave_band_frequencies  Lp_10_meters ]
%
%  The value of R is infinite.  The machine is outside in open air.

octave_band_IL = Lp_10_meters - 30;
    % [ octave_band_frequencies  octave_band_IL ]



%% Anonymous Function for Insertion Loss

h_IL_large = @( Sw, alpha_w, Si, alpha_i, TL )  10*log10(  1  +  (Sw*alpha_w  +  Si*alpha_i)./(Sw + Si)*10^(TL/10)  );



%% Find Values of TL and Aborption that will Meet the Target Insertion Loss - Ground Reflecting

% Assumption(s):
%
%   1.)  The ground is a hard reflecting survice
%   2.)  The enclosure is a cube.
%   3.)  There is no noise transmission through the ground.

enclosure.dimension = 2;  % m
    enclosure.area = 6 * enclosure.dimension^2;  % 20 m^2


% https://www.controlnoise.com/wp-content/uploads/2022/02/Acoustic-Enclosures-Datasheet.pdf
% https://www.cecoenviro.com/wp-content/uploads/2023/12/Acoustic-Enclosures-8pp-A4-web.pdf
% https://www.controlnoise.com/product/acoustic-enclosures/


switch ( 3 )

    case 1

        % 250 Hz - QBV-2
        alpha_w = 0.27;  % From specification sheet.
        TL = 20;  % From specification sheet.
            IL_estimates(1) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 500 Hz - QBV-2
        alpha_w = 0.96;  % From specification sheet.
        TL = 29;  % From specification sheet.
            IL_estimates(2) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 1 kHz - QBV-2
        % alpha_w = 1.13;  % From specification sheet.
        alpha_w = 0.99;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 40;  % From specification sheet.
            IL_estimates(3) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 2 kHz - QBV-2
        % alpha_w = 1.08;  % From specification sheet.
        alpha_w = 0.99;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 50;  % From specification sheet.
            IL_estimates(4) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 4 kHz - QBV-2
        alpha_w = 0.99;  % From specification sheet.
        TL = 55;  % From specification sheet.
            IL_estimates(5) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );


    case 2

        % Note:  Abosrption values are carried over from QBV-2.

        % 250 Hz - QBV-3
        alpha_w = 0.27;  % From specification sheet.
        TL = 25;  % From specification sheet.
            IL_estimates(1) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 500 Hz - QBV-2
        alpha_w = 0.96;  % From specification sheet.
        TL = 33;  % From specification sheet.
            IL_estimates(2) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 1 kHz - QBV-2
        % alpha_w = 1.13;  % From specification sheet.
        alpha_w = 0.99;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 46;  % From specification sheet.
            IL_estimates(3) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 2 kHz - QBV-2
        % alpha_w = 1.08;  % From specification sheet.
        alpha_w = 0.99;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 53;  % From specification sheet.
            IL_estimates(4) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 4 kHz - QBV-2
        alpha_w = 0.99;  % From specification sheet.
        TL = 58;  % From specification sheet.
            IL_estimates(5) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );


    case 3

        % Note:  Abosrption values are carried over from QBV-2.

        % 250 Hz - QBV-3
        alpha_w = 0.99;  % From specification sheet.
        TL = 39;  % From specification sheet.
            IL_estimates(1) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 500 Hz - QBV-2
        alpha_w = 0.96;  % From specification sheet.
        TL = 59;  % From specification sheet.
            IL_estimates(2) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 1 kHz - QBV-2
        % alpha_w = 1.13;  % From specification sheet.
        alpha_w = 0.99;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 68;  % From specification sheet.
            IL_estimates(3) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 2 kHz - QBV-2
        % alpha_w = 1.08;  % From specification sheet.
        alpha_w = 0.99;  % From specification sheet.  See comment on slide 28 of Lecture 10.
        TL = 67;  % From specification sheet.
            IL_estimates(4) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );
        
        % 4 kHz - QBV-2
        alpha_w = 0.91;  % From specification sheet.
        TL = 72;  % From specification sheet.
            IL_estimates(5) = h_IL_large( enclosure.area, alpha_w, machine.area, machine.absorption, TL );

end


[ octave_band_IL    IL_estimates.'    (octave_band_IL - IL_estimates.') ]


% What is the most restrictive case?



%% Find Values of TL and Aborption that will Meet the Target Insertion Loss - Ground with Cover

% Assumption(s):
%
%   1.)  The ground is covered with the absorption material.
%   2.)  The enclosure is a cube.



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


