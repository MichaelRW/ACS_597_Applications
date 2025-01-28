%% Sound Pressure Level Along A Duct With A Temperature Gradient Using The 4-Pole Method
% This Matlab script can be used to plot the pressure
% distribution along a cylindrical duct, where there is a 
% driving velocity at one end (caused by the motion of a piston)
% and the other end is rigid (closed),
% and a _linear temperature gradient_ of the gas along the duct.
%
%% Introduction
% The 4 pole (or transmission line) analysis method is used to calculate
% the pressure distribution along a duct.
%
% The theory and equations are adapted from: 
%
% * R. I. Sujith (1996),
%   "Transfer matrix of a uniform duct with an axial mean
%   temperature gradient",
%   The Journal of the Acoustical Society of America,
%   vol 100, no 4, p2540-2542.
% * C.Q. Howard (2013),
%   "Comments on: "Transfer matrix of a uniform duct with an axial mean
%   temperature gradient", [J. Acoust. Soc. Am., 100, 2540-2542 (1996)].
%
% Note there are several errors in Eqs. (13)-(16) in the Sujith paper. 
% These equations have been re-derived and corrected in the book.
%
% Other relevant references are:
%   
% * Leo L. Beranek and Istvan L. Ver (1992),
%   "Noise and Vibration Control Engineering: Principles and Applications",
%   John Wiley and Sons, New York, USA, ISBN: |0471617512|, 
%   page 377.
% * M. L. Munjal and M. G. Prasad (1986),
%   "On plane-wave propagation in a uniform pipe in the presence of a mean
%   flow and a temperature gradient",
%   The Journal of the Acoustical Society of America,
%   vol 80, no 5, p1501-1506. 
% * K.S. Peat (1988),
%   "The transfer matrix of a uniform duct with a linear temperature
%   graident",
%   Journal of Sound and Vibration
%   vol 123, no 1, p45-53.
%
% It is mentioned in the Sujith paper that 
% the methods by Munjal and Prasad, and the Peat papers are only valid
% for small temperature gradients. 
%
% Keywords: 
% pipe ; duct ; tube ; 4 pole ;  four pole ; transmission line ;
% pressure distribution
%
% Carl Howard 12 March 2012

%close all; clear all

%% Properties of the fluid
c_speed_sound=343.34;   % [m/s] Speed of sound in air at 20 degree c 
rho_density=1.2021;     % [kg/m^3] Density of air 

gamma = 1.4;            % ratio of specific heats
R_universal = 8.3144621;% Universal gas constant
M_molar_mass = 0.029;   % Molar mass of air

% Specific gas constant
R_specific = R_universal / M_molar_mass;

P_atmospheric =  101325; % [Pa] Atmospheric pressure

%% Analysis frequency
% The analysis frequency is $f$, and 
% the circular frequency is $\omega=2 \pi f$
%
freq=200;               % [Hz] Frequency of analysis
omega = 2*pi*freq;      % [radians /s] circular frequency

%%% 
% Wavelength $\lambda$
lambda = c_speed_sound / freq;

%%%
% Wave number $k$
k_wave_no=(2*pi*freq)/c_speed_sound;

%% Properties of the duct
a_radius=0.1/2;                    % [m] duct radius

%%%
% Cross sectional area of circular duct $S$
S_duct_area=pi*a_radius^2;     % [m]

%%%
% Length of duct $L$
L_duct_len=3.0;                   % [m] Length of the duct

%%%
% Define a vector containing the positions of where pressures will be
% calculated (microphone location) in increments of $\Delta z$
number_divs = 1000;                 %[no units] number of duct segments
delta_z=L_duct_len/number_divs;   %[m] length of each segment
mic_loc=[0:delta_z:L_duct_len];   %[m] microphone locations


%% Define the acoustic excitation
% Define the particle velocity at source and rigid ends of the duct
u_particle_vel_1=0;     %[m/s] closed rigid end at z=L_duct_len [m]
u_particle_vel_2=1;     %[m/s] driving end at z=0 [m]

%%%
% Acoustic particle velocity
%
% In Beranek and Ver, they used mass volume velocity $\rho S u$. 
% Munjal and Prasad, and Peat, they use $S u$, because the density $\rho$
% varies with location and this term is included in the 4-pole transmission
% matrix. In Sujith, he uses just particle velocity $u$.
V_1=u_particle_vel_1;     %[kg /s] at z=L_duct_len [m]
V_2=u_particle_vel_2;     %[kg /s] at z=0 [m]


%% Temperature at each end of the duct
% Note that the formulation used by Sujith cannot handle
% $T_1=T_2$. If the user sets $T_1=T_2$ then $m=0$ and this will cause 
% equations to attempt to evaluate 0/0 which will generate |NaN| errors.
% In addition the term $\nu$ will evaluate to 0, and cause
% issues with the terms $1/ \nu$.

% This is the temperature range used in the Ansys Workbench and 
% Ansys Mechanical APDL models.
Temp_1 = 400 + 273;      % [Kelvin] at z=L_duct_len [m]
Temp_2 = 20 + 273;      % [Kelvin] at z=0 [m]


% Temp_1 = 400.1 + 273;      % [Kelvin] at z=L_duct_len [m]
% Temp_2 = 400 + 273;      % [Kelvin] at z=0 [m]

% This is for a constant temperature duct
% Temp_1 = 22.1 + 273;        % [Kelvin]  at z=L_duct_len [m]
% Temp_2 = 22 + 273;          % [Kelvin] at z=0 [m]

if (Temp_1-Temp_2==0),
    error('This script cannot handle T1=T2.')
end;


%%%
% The gas has a linear temperature gradient between the ends of the duct of
% the form
%
% $$ T(z) = T_2 + m z $$
%
% where
%
% $$ m = \frac{T_1 - T_2}{L} $$
%
% It is assumed that at $z_1=L$ and $z_2=0$.
%
m_temp_grad = (Temp_1 - Temp_2) / L_duct_len;


%%%
% Define the densities
% From Bies and Hansen (2009) Eq. 1.8, p17-18
% 
% $$ \rho = \frac{M P_{\textrm{static}}}{R T} $$
%
% where
%
% * $P_{\textrm{static}} = 101325 \, \textrm{Pa}$ atmospheric pressure 
%   (assuming the gas in the tube is not pressurised.)
% * $M=0.029 \, \textrm{kg.mol}^{-1}$ molecular weight of air
% * $R=8.314 \, \textrm{J.mol}^{-1}.\,{}^{\circ} \textrm{K}^{-1}$ 
%   universal gas constant
% * $T$ temperature in Kelvin.
%
rho_1 = P_atmospheric * M_molar_mass / ( R_universal * Temp_1);
rho_2 = P_atmospheric * M_molar_mass / ( R_universal * Temp_2);


%% Calculate pressure at the closed end of the duct
% Define the 4-pole transmission matrix
%
% First define the constant $\nu$ [ Sujith, Eq. (5)]
nu = abs(m_temp_grad)/2 * sqrt(gamma*R_specific);

%%%
% Calculate two constant terms used numerous times in the evaluation of the
% elements of the transmission matrix
omega_nu_sqrt_Temp_1 = omega/nu*sqrt(Temp_1);
omega_nu_sqrt_Temp_2 = omega/nu*sqrt(Temp_2);

%%%
% For a duct with a temperature gradient the terms of the transmission
% matrix are given by [see the book or Howard (2013), the equations in
% Sujith (1996) have several errors].
%
T_11 = pi/(2*nu) * omega * sqrt(Temp_1) *        ...
        (           ...
            besselj(1,omega_nu_sqrt_Temp_1)     ...
            * bessely(0,omega_nu_sqrt_Temp_2)     ...
            - besselj(0,omega_nu_sqrt_Temp_2)     ...
            * bessely(1,omega_nu_sqrt_Temp_1)     ...
        );
    
T_12 =  1i* pi/(2*nu) * omega * sqrt(Temp_1)      ...
        * abs(m_temp_grad)/m_temp_grad                                      ...
        * rho_1 * sqrt(gamma*R_specific*Temp_1)             ...
        * (                                             ...
           besselj(0,omega_nu_sqrt_Temp_2)          ...
            * bessely(0,omega_nu_sqrt_Temp_1)     ...
            - besselj(0,omega_nu_sqrt_Temp_1)     ...
            * bessely(0,omega_nu_sqrt_Temp_2)     ...
        );
    
T_21 =  1i * pi/(2*nu) * omega * sqrt(Temp_1)           ...
        * m_temp_grad / abs(m_temp_grad)                      ...
        /( rho_2*sqrt(gamma*R_specific*Temp_2) )        ...
        * (                                             ...
           besselj(1,omega_nu_sqrt_Temp_2)          ...
            * bessely(1,omega_nu_sqrt_Temp_1)     ...
            - besselj(1,omega_nu_sqrt_Temp_1)     ...
            * bessely(1,omega_nu_sqrt_Temp_2)     ...
        );

T_22 = pi/(2*nu) * omega * sqrt(Temp_1)           ...
        * rho_1 * sqrt(gamma*R_specific*Temp_1)         ...
        /( rho_2* sqrt(gamma*R_specific*Temp_2) )        ...
        * (                                             ...
           besselj(1,omega_nu_sqrt_Temp_2)          ...
            * bessely(0,omega_nu_sqrt_Temp_1)     ...
            - besselj(0,omega_nu_sqrt_Temp_1)     ...
            * bessely(1,omega_nu_sqrt_Temp_2)     ...
        );



T = [ T_11    T_12    ;   T_21    T_22];


%%%
% Calculate pressure at the closed end of the duct $P_1$
P_1 = (V_2 - T_22 * V_1) / T_21;
P_2 = T_11 * P_1 + T_12 * V_1;

P_closed_end=P_1;
P_driving_end=P_2;


%% Calculate the pressure distribution along the length of the duct
% We know the pressure and velocity at the driving end of the duct 
% $P_2, u_2$ and at the closed end $P_1, u_1$. We can define a 4-pole or
% transmission matrix for a small length of duct of length $\Delta z$ and
% incrementally calculate the pressure and velocity at the start of the
% duct.
%
% *Note* that the difference between this analysis with a temperature
% gradient along the duct, and with a constant temperature profile, is that
% the transmission matrix will vary along the length of the duct, as the
% properties of the gas changes at each location.
%
% We will be effectively traversing a "microphone" 
% from the closed end of the duct to the driving end and calculating
% the pressure and velocity at each point.
%
% There are |number_divs| small duct segments 
% (or divisions). 
% The location of the "microphone" is 
%
% $$ z_n = L - n \times \Delta z$$
%
% where $n = 1 \cdots$ |number_divs|. 
% Note that we do not bother to calculate the response again 
% at $n=0$, $z=L$ because we already know that it is $[P_1 ; u_1]$.
%
% The sound pressure and velocity at microphone location $z_n$ is
%
% $$ \left[ \matrix{ P_{\textrm{mic}} \cr u_{\textrm{mic}} } \right] =
%           \left[\mathbf{T}_{\Delta z}\right]^{n} 
%       \left[ \matrix{ P_1 \cr u_1 } \right] $$
% 
% where $n = 1 \cdots$ |number_divs|. 
%
% Define an empty matrix that contains the:
%
% * pressure (top row) and
% * velocity (bottom row) 
% along the length of the duct.
p_v=zeros(2,length(mic_loc));

%%%
% Insert the first column of the matrix with the pressure and velocity at
% the driving end.
p_v(:,end) = [ P_1 ; V_1];


%%%
% Calculate the response along the length of the duct.
% Note that there are |number_divs+1| microphone positions.
% The constants will have to re-calculated at each location, and we'll 
% name them with the suffix |_dash|.
for index=length(mic_loc)-1:-1:1

    Temp_1_dash = (Temp_1 - Temp_2)/L_duct_len*mic_loc(index+1)  ...
                    + Temp_2;
    Temp_2_dash =(Temp_1 - Temp_2)/L_duct_len*mic_loc(index)  ...
                    + Temp_2;

    rho_1_dash = P_atmospheric * M_molar_mass / ( R_universal * Temp_1_dash);
    rho_2_dash = P_atmospheric * M_molar_mass / ( R_universal * Temp_2_dash);

    m_temp_grad_2 = (Temp_1_dash - Temp_2_dash) / delta_z;

    nu_dash = abs(m_temp_grad_2)/2 * sqrt(gamma*R_specific);

    omega_nu_sqrt_Temp_1_dash = omega/nu_dash*sqrt(Temp_1_dash);
    omega_nu_sqrt_Temp_2_dash = omega/nu_dash*sqrt(Temp_2_dash);


    T_11 = pi/(2*nu_dash) * omega * sqrt(Temp_1_dash) *     ...
            (                                               ...
                besselj(1,omega_nu_sqrt_Temp_1_dash)        ...
                * bessely(0,omega_nu_sqrt_Temp_2_dash)      ...
                - besselj(0,omega_nu_sqrt_Temp_2_dash)      ...
                * bessely(1,omega_nu_sqrt_Temp_1_dash)      ...
            );

    T_12 =  1i* pi/(2*nu_dash) * omega * sqrt(Temp_1_dash)      ...
            * abs(m_temp_grad)/m_temp_grad                      ...
            * rho_1_dash * sqrt(gamma*R_specific*Temp_1_dash)   ...
            * (                                                 ...
               besselj(0,omega_nu_sqrt_Temp_2_dash)             ...
                * bessely(0,omega_nu_sqrt_Temp_1_dash)          ...
                - besselj(0,omega_nu_sqrt_Temp_1_dash)          ...
                * bessely(0,omega_nu_sqrt_Temp_2_dash)          ...
            );

    T_21 =  1i * pi/(2*nu_dash) * omega * sqrt(Temp_1_dash)     ...
            * m_temp_grad / abs(m_temp_grad)                    ...
            /( rho_2_dash*sqrt(gamma*R_specific*Temp_2_dash) )  ...
            * (                                                 ...
               besselj(1,omega_nu_sqrt_Temp_2_dash)             ...
                * bessely(1,omega_nu_sqrt_Temp_1_dash)          ...
                - besselj(1,omega_nu_sqrt_Temp_1_dash)          ...
                * bessely(1,omega_nu_sqrt_Temp_2_dash)          ...
            );

    T_22 = pi/(2*nu_dash) * omega * sqrt(Temp_1_dash)           ...
            * rho_1_dash * sqrt(gamma*R_specific*Temp_1_dash)   ...
            /( rho_2_dash* sqrt(gamma*R_specific*Temp_2_dash) ) ...
            * (                                                 ...
               besselj(1,omega_nu_sqrt_Temp_2_dash)             ...
                * bessely(0,omega_nu_sqrt_Temp_1_dash)          ...
                - besselj(0,omega_nu_sqrt_Temp_1_dash)          ...
                * bessely(1,omega_nu_sqrt_Temp_2_dash)          ...
            );

   
    T_dash = [ T_11    T_12    ;   T_21    T_22];

    p_v(:,index)=T_dash*p_v(:,index+1);
end;

%% Calculate Sound Pressure Levels
% The sound pressure level at the microphone location is calculated as
%
% $$ \textrm{SPL}_{\textrm{mic}} = 
%     20 \log_{10} 
%           \left( 
%               \frac{|p_{\textrm{RMS}}|}{20 \times 10^{-6} } 
%           \right)     $$
%
% $$    \quad
%   = 20 \log_{10} 
%           \left( 
%               \frac{|p|}{\sqrt{2}} \times \frac{1}{20 \times 10^{-6} } 
%           \right)     
% $$
SPL_mic=20*log10(abs(p_v(1,:)/sqrt(2))/20e-6);

%%%
% For comparison, the SPL at the driving and closed end of the duct are
SPL_driving_end=20*log10(abs(P_2)/20e-6/sqrt(2));
SPL_closed_end=20*log10(abs(P_1)/20e-6/sqrt(2));

%% Plot the SPL along the length of the duct
figure
set(gca,'fontsize',16);
plot1h=plot(    mic_loc,SPL_mic,    ...
            0,SPL_driving_end,'x',    ...
            L_duct_len,SPL_closed_end,'o');
set(plot1h,'linewidth',2);
set(plot1h(2),'markersize',10);
xlabel('Mic Position Along Duct [m]'); 
ylabel('Sound Pressure Level [dB re 20uPa]');
title('Sound Pressure Level Along A Piston-Rigid Duct');
leg1h=legend( 'SPL along duct',       ...
        'SPL at driving end',   ...
        'SPL at rigid end',    ...
        'location','southwest');
axis([0 L_duct_len 100 150]);
box on
grid;

%% Plot the pressure along the length of the duct
figure
plot2h=plot(mic_loc,real(p_v(1,:)),     ...
            mic_loc,imag(p_v(1,:)),'--',     ...
            0,imag(P_2),'x',            ...
            L_duct_len,imag(P_1),'o');


set(gca,'fontsize',16);
set(plot2h,'linewidth',2);
set(plot2h(3),'markersize',10);
set(plot2h(4),'markersize',10);
xlabel('Mic Position Along Duct [m]'); 
ylabel('Pressure [Pa]');
title('Pressure Along A Piston-Rigid Duct');

leg2h=legend( 'Theory: Real pressure',       ...
        'Theory: Imag pressure',       ...
        'Im(P_2) at driving end',   ...
        'Im(P_1) at rigid end',    ...
        'location','EastOutside');
    
axis([0 L_duct_len -600 600]);
box on
grid;


%% Plot the particle velocity along the duct
u_mic = p_v(2,:);

figure
plot3h=plot(    mic_loc,real(u_mic),    ...
            mic_loc,imag(u_mic),'--',    ...
            0,u_particle_vel_2,'x',    ...
            L_duct_len,u_particle_vel_1,'o');
set(gca,'fontsize',16);
set(plot3h,'linewidth',2);
set(plot3h(3),'markersize',10);
set(plot3h(4),'markersize',10);
xlabel('Mic Position Along Duct [m]'); 
ylabel('Particle Velocity [m/s]');
title('Particle Velocity Along A Piston-Rigid Duct');
leg3h=legend( 'Theory: Real particle vel.',       ...
            'Theory: Imag particle vel.',       ...
        'u_2 at driving end',   ...
        'u_1 at rigid end',    ...
        'location','EastOutside');
box on
grid on;
