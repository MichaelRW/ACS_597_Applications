

%% Environment

close all;  clear;  clc;
% restoredefaultpath;

format LongG;

pause( 1 );



%% Bending Wave Speed

c = 343;

f = 0.1:1:20e3;
    w = 2*pi*f;

h_bending_wave_speech = @( D, w, ms )  sqrt( ( D.*w ) / ms );


figure( ); ...
    plot( f, ones( size(f) )*c, 'Marker', 'none' );  hold on;
    plot( f, h_bending_wave_speech( 31.4, w, 9.36 ) );  grid on;
