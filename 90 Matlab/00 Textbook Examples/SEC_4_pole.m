%% Transmission Loss of an Expansion Chamber Silencer
% This script is used to calculate the transmission loss of a 
% single expansion chamber silencer
%%
% The goal of this script is to generate Fig 10.12, p386 from
% the book:
%
% * Beranek and Ver (1992),
%   "Noise and Vibration Control Engineering:
%   Principles and Applications",
%   John Wiley and Sons, New York, USA.
%   ISBN: 0-471-61751-2
%
% Keywords: 
% pipe ; duct ; 4 pole ;  four pole ; transmission line ;
% silencer; muffler; expansion chamber
%
% Carl Howard 24 April 2013

%% Define properties of the fluid
c = 343.24;          %[m/s] speed of sound
rho = 1.2041;        %[kg/m^3] density of air

%% Define the properties of the main duct
a = 0.05/2;        %[m] radius of main duct
S_duct = pi*(a^2);       %[m^2] area of main duct

%% Define the properties of the quarter wave tube
% Beranek and Ver (1992) Fig 10.12 shows the TL evaluated at 
% area ratios of 
% |N_ratio_areas| = 2, 4, 8, 16, 32, 64
% where 
%
% $$ 
%   N = \frac{S_2}{S_1} 
%   = \frac{S_{\textrm{expansion}}}{S_{\textrm{duct}}} 
% $$
%
% where 
%
% * $S_{\textrm{expansion}}$ is the cross sectional area of the 
% expansion chamber and
% * $S_{\textrm{duct}}$ is the cross sectional area of the main
% duct.
%
N_ratio_areas = 16;                  %[] ratio of areas
S_expansion = N_ratio_areas * S_duct;%[m^2] area of expansion chamber
L_expansion = 0.5;                   %[m] length of the expansion chamber

% S_expansion=pi*(0.2/2)^2;
% N_ratio_areas=S_expansion/S_duct;

%%%
% This is what would normally be used to define the area of the 
% expansion chamber.
%
%   a_expansion = 0.1;               %[m] radius of the expansion chamber
%   S_expansion = pi*a_expansion^2;    %[m^2] area of expansion chamber

%% Define the analysis frequency range
% Beranek and Ver (1992) Fig 10.12 has the x-axis range from 
% $(kL / \pi) = 0 \cdots 4$.
% Hence this determines the analysis frequency range.
k_max = 4 *pi / L_expansion;      %[] maximum wavenumber of analysis
f_max = k_max * c / (2*pi);     %[Hz] maximum analysis frequency
freq_spacing=1;               %[Hz] frequency increment
f = [1:freq_spacing:f_max];     %[Hz] Excitation frequency 
omega= 2*pi*f;                  %[rad/s] circular frequency
k = omega / c;                  % vector containing wavenumbers
kL = k*L_expansion;               % a parameter used numerous times

%% Calculate the transmission loss
% The transfer matrix for a resonator is given by
% [see Beranek and Ver (1992)], 
%
% $$ \mathbf{T} = 
%       \left[ 
%           \matrix{ T_{11} & T_{12}   \cr 
%                  T_{21} & T_{22} }  
%       \right]
%   = \left[ 
%       \matrix{
%       \cos (k L)  &   \frac{ j c_0 }{ S_{\textrm{expansion}} } \sin (k L)      \cr
%       \frac{j S_{\textrm{expansion}} }{c_0} \sin (k L)  &   \cos(k L)
%       }
%       \right] 
%   $$
%
T_11 = cos(kL);
T_12 = 1i*c/S_expansion * sin(kL);
T_21 = 1i*S_expansion/c * sin(kL);
T_22 = cos(kL);


%% Calculate the transmission loss
% The transmission loss is calculated from 
% Beranek and Ver (1992), Eq. (10.10) p374 as
%
% $$\textrm{TL} = 20 \log_{10} 
%   \left| 
%       \frac{T_{11}+ \frac{S}{c} T_{12} + \frac{c}{S} T_{21} + T_{22}}{2}
%   \right|     $$
%
TL_4pole = 20*log10( abs( T_11+S_duct/c*T_12 +c/S_duct*T_21 + T_22)/2);


%%%
% Jacobsen (2011) defines a single line equation for the 
% transmission loss of an expansion chamber as 
% [ Finn Jacobsen, (2011), September,
% "Propagation of sound waves in ducts", 
% Denmark Technical University, Note no 31260,
% <http://web-files.ait.dtu.dk/fjac/p_home_page/notes/Duct_acoustics.pdf> ]
%
% [see also Bies and Hansen (2009),
% "Engineering Noise Control", Eq. (9.99), p464]
%
% $$ \textrm{TL} = 10 \times \log_{10} 
%   \left[ 1 + \frac{1}{4} 
%       \left( \frac{S_1}{S_2} - \frac{S_2}{S_1} \right)^2
%               \sin^2 kL 
%   \right]
% $$
%
% where 
%
% * $S_1=$ |S_duct| is the cross sectional area of the 
%   upstream and downstream ducts.
% * $S_2=$ |S_expansion| is the cross sectional area of the 
%   expansion chamber, and is always $S_2 > S_1$
% * $L$ is the length of the expansion chamber.
%
TL_B_and_H = 10*log10( 1 + 0.25*(1/N_ratio_areas - N_ratio_areas)^2    ...
                    *( sin(kL) ).^2   );


%% Plot the transmission loss
figno=figure(1);
x_axis=k*L_expansion/pi;          % normalised frequency
%p1h=plot(x_axis,TL_4pole,'-',x_axis(1:20:end),TL_B_and_H(1:20:end),'r.');
p1h=plot(x_axis,TL_4pole,'k:');
grid on

% Set the page size for the figure
%set(figno, 'units', 'centimeters', 'pos', [0 0 10 10])
set(figno, 'units', 'centimeters', 'pos', [0 0 10 8])

% Set some of the font properties
hA=gca;

fontname ='Times New Roman';
set(hA,'defaultaxesfontname',fontname);
set(hA,'defaulttextfontname',fontname);

% make gridlines more visible
set(hA,'GridAlpha',0.8);       

set(hA,'fontsize',9);
set(p1h,'linewidth',1);
% l1h=legend(  ['Beranek: N=',num2str(N_ratio_areas)],   ...
%             ['B&H: N=',num2str(N_ratio_areas)]);uistack(p1h(2),'top');
l1h=legend(  ['SEC: $N=',num2str(N_ratio_areas),'$']);

l1h=legend(  ['SEC: $N=2$'],    ...
             ['SEC: $N=4$'],    ...
             ['SEC: $N=8$'],    ...
             ['SEC: $N=16$']    );


set(l1h,'Interpreter','Latex','FontSize',9,'FontName','Times New Roman');        

axis([0 4 0 30]);

xlabel('Normalised Frequency $kL/\pi$','Interpreter','latex','FontName','Times New Roman');
ylabel('Transmission Loss [dB]','Interpreter','latex','FontName','Times New Roman');


set(hA,'FontUnits','points','FontWeight','normal', ...
    'FontSize',9,'FontName','Times New Roman');
