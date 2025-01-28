%% Radiation from a Quadrupole Longitudinal Orientation
% This script will plot the theoretical sound pressure level radiated by a
% longitudinal quadrupole versus angle.
%
% Reference 
%
% * David A. Bies and Colin H. Hansen (2009), 
%   "Engineering Noise Control: Theory and Practice", 
%   Fourth edition, Spon Press.
%
% Carl Howard 10 January 2013

%% Theory
% The sound pressure from a dipole source is given by
% 
% $$
%   p_{\textrm{dipole}}  = -j \frac{ \rho_0 c_0 k^2 S \, u \, d}{
%                                   4 \pi r} \, \cos(\theta) \,
%                                  e^{j(\omega t - kr)}
% $$
%    
% where 
%
% * $S$ is the surface area of the sphere of radius $a$ given by
% 
% $$    S = 4 \pi a^2   $$
%
% * $u$ is the velocity of the expanding and contracting 
% surface of the sphere, 
% * $\rho_0$ is the density of the acoustic fluid, 
% * $c_0$ is the speed of sound of the acoustic fluid,
% * $\omega=2 \pi f$ is the angular frequency, 
% * $k=\omega/c_0$ is the wavenumber, 
% * $r$ is the distance from the centre of the dipole to the 
% measurement location, 
% * $d$ is the distance between the two out of phase monopole sources.

%% Define the parameters of the dipole

%%% 
% Radius of the monopole sphere $a$
a_radius        = 0.01;         %[m]

%%%
% Velocity of the surface of the sphere $u$
u_velocity      = 1.0;        % [m/s]

%%%
% Separation distance between the out-of-phase monopoles $d$
h_x_separation    = 0.1;          % [m]
L_y_separation    = 0.1;        % [m]

%%%
% Surface area of each monopole sphere
S_area          = 4*pi*a_radius^2;

%%%
% Volume velocity of each monopole in the dipole $Q = S \, u$
Q_volume_velocity = S_area * u_velocity;


%% Properties of the fluid
%

%%% 
% Density of the fluid $\rho$
rho_density     =   1.2041;     % [kg/m^3]

%%% 
% Speed of sound of the fluid $c$
c_speed_sound   =   343.024;    % [m/s]

%%% 
% Reference Sound Pressure $p_{\rm{ref}}$
p_ref           =   20e-6;      % [Pa]

%% Analysis configuration
% Harmonic frequency of analysis
freq            = 500;          % [Hz]
omega           = 2*pi*freq;    % [radians/sec]
k_wavenum       = omega/c_speed_sound; 

%%% 
% Angle over which the SPL is to be calculated $\theta$.
% For the |Dirplot| matlab function, the angle |theta| 
% has to be in degrees, and for a full circle 
% is between -180 to 180 degrees
theta           = [-180:1:180]; % [degrees]

%%% 
% Radius where SPL is to be evaluated
r_radius        = 4;            % [m]


%% Calculate the sound pressure of a quadrupole
% The sound pressure from a dipole source is given by
% 
% $$
%   p_{\textrm{dipole}}  = -j \frac{ \rho_0 c_0 k^2 S \, u \, d}{
%                                   4 \pi r} \, \cos(\theta) \,
%                                  e^{j(\omega t - kr)}
% $$
%
p_quadrupole      = 0;


% From Eq. 5.65 Bies and Hansen edition 5
Q_rms = Q_volume_velocity / sqrt(2);

% Eq 5.68
W_long_quad = (rho_density*c_speed_sound) * ...
        (2*k_wavenum.^3 * h_x_separation * L_y_separation * Q_rms)^2 /    ...
        (5*pi * (1+k_wavenum^2 * a_radius^2) );


% From Eq. 5.69 Bies and Hansen edition 5
% Define theta=0, so that it is in the plane of the lateral quadrupoles
p_squared=  (5* rho_density*c_speed_sound ) *     ...
            W_long_quad* (cos(theta*pi/180)).^4  / ...
            (4*pi*r_radius^2);




                    
%%%
% The sound pressure level in decibels is given by
%
% $$ L_p = 20 \log_{10} 
%               \left[
%                   \frac{p_{\rm{dipole}}}{\sqrt{2} \, p_{\rm{ref}}}
%               \right]
% $$
%
SPL_long_quad      = 10*log10(p_squared/p_ref.^2);



%% Rescale the SPL plot
% For the book we will rescale the results so that it looks like a
% directivty plot, by normalising the SPL to the maximum SPL. This is not
% strictly correct as directivity should be normalised by the 
% average intensity.
max_SPL = max(SPL_long_quad);

SPL_long_quad_norm=SPL_long_quad-max_SPL;

%Alter the results where the SPL_dipole_norm < -100 to SPL=0
indexes=find(SPL_long_quad_norm<-60);
SPL_long_quad_norm(indexes)=-60;

%% Plot the SPL versus angle

p1h=polarplot( (theta*pi/180),SPL_long_quad_norm,'k-');

% Set the page size for the figure
set(1, 'units', 'centimeters', 'pos', [0 0 10 10])

set(p1h(1),'linewidth',3);

% % Set some of the font properties
hA=gca;
fontname ='Times New Roman';
set(hA,'FontName',fontname);
set(hA,'FontSize',9);
set(hA,'defaultaxesfontname',fontname);
set(hA,'defaulttextfontname',fontname);

% Set orientation of 0 degrees. Default is on the horizontal = 'right'
% If you want 0 deg on vertical top, then change to 'top'
set(hA,'ThetaZeroLocation','top');    

% Change the direction of the angles to clockwise is increasing.
set(hA,'ThetaDir','clockwise');

% Set the limits for the radial axis
set(hA,'RLim',[-60 0]);    

% Set angle for Radial Axis numbers to 15 deg
set(hA,'RAxisLocation',15);  

% Try to move position of tick labels
rTL=get(hA,'RTickLabels');
% for jj=1:length(rTL)
%     rTL_v2{jj,1}=strjoin({'     ',rTL{jj}});
% end;
% Delete the first label at the origin
rTL_v2=rTL;
rTL_v2{1}='';
set(hA,'RTickLabels',rTL_v2);


% make gridlines more visible
set(hA,'GridAlpha',0.8);       

%% Plot another version
%
% From 
% Fundamentals of General Linear Acoustics
% By Finn Jacobsen, Peter Moller Juhl
% Eq 9.15
%
% Equation 9.15
p_long_quad_v2 = -1i * rho_density * c_speed_sound * (h_x_separation)^2 * ...
        k_wavenum^3 * Q_volume_velocity / (pi*r_radius) *   ...
        exp(-1i*k_wavenum*r_radius) *   ...
        (   ...
        cos(theta*pi/180).^2 * ...
        ( 1 + 3/(1i*k_wavenum*r_radius) - 3/(k_wavenum*r_radius)^2 ) ...
        - ( 1 / (j*k_wavenum*r_radius) )    +   ...
         ( 1 / (k_wavenum * r_radius)^2 )   ...
        );
    
SPL_long_quad_v2 = 20*log10( abs(p_long_quad_v2)/sqrt(2)/p_ref);

figure(2)
p2=plot(theta,SPL_long_quad,theta,SPL_long_quad_v2,'.');
l2h=legend('Bies and Hansen','Jacobssen');
axis([-180 180 40 90]);
xlabel('Angle Theta [degrees]');
ylabel('SPL [dB re 20 micro-Pa]');
