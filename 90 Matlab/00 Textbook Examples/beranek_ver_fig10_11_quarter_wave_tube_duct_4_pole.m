%% Transmission Loss of a Side Branch Resonator
% This script is used to calculate the transmission loss of a 
% quarter wave tube (also known as a transverse tube) resonator, 
% attached to a main duct.
%
%%
% The goal of this script is to generate Fig 10.11, p384 from
% the book:
%
% *   Beranek and Ver (1982),"Noise and Vibration Control Engineering:
%       Principles and Applications".
%
% Keywords: 
% pipe ; duct ; 4 pole ;  four pole ; transmission line ;
% quarter wave tube; transverse tube
%
% Carl Howard 18 December 2012

%% Define properties of the fluid
c = 343.24;          %[m/s] speed of sound
rho = 1.2041;        %[kg/m^3] density of air

%% Define the properties of the main duct
L = 4;              %[m] length of main duct
a = 0.100/2;        %[m] radius of main duct
S = pi*(a^2);       %[m^2] area of main duct

%% Define the properties of the quarter wave tube
% Beranek and Ver (1992) Fig 10.11 shows the TL evaluated at 
% area ratios of 
% |N_ratio_areas| = 0.25, 0.5, 1, 2, 4, 8,
% where 
%
% $$N = \frac{S_2}{S_1} = \frac{S_{\textrm{QWT}}}{S_{\textrm{duct}}} $$
%
% where $S_{\textrm{QWT}}$ is the cross sectional area of the quarter wave
% tube, and $S_{\textrm{duct}}$ is the cross sectional area of the main
% duct.
%
N_ratio_areas = 0.5;                %[] ratio of areas
S_quarter = N_ratio_areas * S;      %[m^2] area of quarter wave tube
L_quarter = 1.45 + 0.61*0.025;                   %[m] length of the quarter wave tube

%%%
% This is what would normally be used to define the area of the quarter
% wave tube.
%a_QWT = 0.1;               %[m] radius of the quarter wave tube
%S_quarter = pi*a_QWT^2;    %[m^2] area of quarter wave tube

%% Define the analysis frequency range
% Beranek and Ver (1992) Fig 10.11 has the x-axis range from 
% $(kL / \pi) = 0 \cdots 4$.
% Hence this determines the analysis frequency range.
k_max = 4 *pi / L_quarter;      %[] maximum wavenumber of analysis
f_max = k_max * c / (2*pi);     %[Hz] maximum analysis frequency
%freq_spacing=0.1;               %[Hz] frequency increment
freq_spacing=0.01;               %[Hz] frequency increment

f = [1:freq_spacing:f_max];     %[Hz] Excitation frequency 
omega= 2*pi*f;                  %[rad/s] circular frequency
k = omega / c;                  % vector containing wavenumbers

%% Calculate the transmission loss
% The transfer matrix for a resonator is given by
% (see Beranek and Ver (1992), Eq. (10.20), p379)
%
% $$ \mathbf{T} = 
%       \left[ 
%           \matrix{ T_{11} & T_{12}   \cr 
%                  T_{21} & T_{22} }  
%       \right]
%   = \left[ 
%           \matrix{ 1 & 0   \cr 
%                   \frac{1}{Z_r} & 1   }  
%       \right] $$
%
% where (see p380 Beranek and Ver (1992))
%
% * $Z_r=Z_t+Z_c$ is the combined impedance of the resonator.
% * $Z_{\textrm{t}}$ is the impedance of the throat connecting the pipe 
%   to the cavity, which is not considered in the analysis for Fig 10.11.
% * $Z_{\textrm{c}}$ is the impedance of the cavity, which is zero in 
%           the case of a quarter wave tube, 
% * $Z_{\textrm{tt}}$ is the impedance of the transverse tube defined as
%       (see Eq. (10.21) Beranek and Ver (1992) )
%
% $$ Z_{\textrm{tt}}=\frac{-jc}{S_{\textrm{QWT}}}
%                       \cot(k L_{\textrm{QWT}}) $$
%
% The transmission loss is calculated from 
% Beranek and Ver (1992), Eq. (10.10) p374 as
%
% $$\textrm{TL} = 20 \log_{10} 
%   \left| 
%       \frac{T_{11}+ \frac{S}{c} T_{12} + \frac{c}{S} T_{21} + T_{22}}{2}
%   \right|     $$
%
% (Note that a more efficient way of calculating this result 
% is to remove the constant values out of the |for| loop, 
% but will be done this way so that it is
% easier to follow the sequence of calculations.)
%

for index=1:length(k);
    Z_tt = -1i*c/S_quarter*cot(k(index)*L_quarter);
    Z_c = Z_tt;
    Z_t = 0;
    Z_r = Z_t + Z_c;

    T_11 = 1;
    T_12 = 0;
    T_21 = 1/Z_r;
    T_22 = 1;

    TL(index) = 20*log10( abs( T_11+S/c*T_12 +c/S*T_21 + T_22)/2);
end;

%% Plot the transmission loss
x_axis=k*L_quarter/pi;          % normalised frequency
figure(1);
p1=plot(x_axis,TL,'k:');
grid on
figno=gcf;
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
set(p1,'linewidth',1);
axis([0 4 0 40]);
xlabel('Normalised Frequency $kL/\pi$','Interpreter','latex','FontName','Times New Roman');
ylabel('Transmission Loss [dB]','Interpreter','latex','FontName','Times New Roman');
%title('Predicted TL of a Transverse Tube Resonator');
l1h=legend(['N=',num2str(N_ratio_areas)]);
set(l1h,'Interpreter','Latex','FontSize',9,'FontName','Times New Roman');        

set(hA,'FontUnits','points','FontWeight','normal', ...
    'FontSize',9,'FontName','Times New Roman');

