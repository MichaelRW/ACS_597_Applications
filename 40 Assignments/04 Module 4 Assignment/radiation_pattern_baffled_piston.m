%% Radiation pattern of a baffled circular piston in an infinite baffle
% This script will plot the beam pattern of a baffled circular piston 
%
% References for this work are:
%
% * Lawrence E. Kinsler, Austin R. Frey, Alan B. Coppens, James V. Sanders,
%   "Fundamentals of Acoustics", 4th edition, Wiley, 1999.
%
% The goal is to generate Figure 7.4.5, p183 in 
% Kinsler et al. (1999)
%
% Carl Howard, 27 November 2012
%
% Filename: |radiation_pattern_baffled_piston.m|

%% Define properties of air
c_speed_air=343;        % [m/s] speed of sound
rho_dens_air=1.21;      % [kg/m^3] density of air

%% Define dimensions of piston
% The piston radius is arbitrarily defined $a=0.1$m
%a_piston_radius=0.1;    % [m] radius of the piston

a_piston_radius=0.01;    % [m] radius of the piston


%% Define the analysis frequency range
% Figure 7.4.5, p183 Kinsler et al. has the frequency
% defined as $ka=10$. Hence $k = 10 / a$
k_wavenum = 10 / (a_piston_radius);     % wavenumber



%%
% Calculating alternative representations of frequency
%
% $$\omega = k c $$
%
omega = k_wavenum*c_speed_air;  % [rad/s] frequency

%%
% and also
%
% $$ f = \omega / (2 \pi) $$
%
freq = omega / (2*pi);          % [Hz] frequency

%%
% The maximum velocity of the piston is arbitrarily defined here as
%U_max_vel = 1e-3*2*pi*freq;     % [m/s] maximum velocity of the piston

U_max_vel = 1;

 
%% Calculate the pressure distribution
% The pressure is given by Eq.(7.4.17) in Kinsler et al.
%
% $$ p(r,\theta,t) = \frac{j}{2} \rho_0 c U_0 \frac{a}{r} k a 
%                      \left[
%                           \frac{2 J_1 (k a \sin \theta)}{ka \sin \theta}
%                       \right]
%                       e^{j(\omega t - kr)}    $$
%
% where $J_1$ is the Bessel function.
%
% Note that the value of $[2 J_1(x)/x]$ is unity when $x=0$.
%
% However, $J_1(0)=0$ and Matlab does not resolve $0/0 = 1$.
%
% The workaround is to make sure that $x \neq 0$.

%%
% Define the range of angles $\theta$ for the beam pattern
theta1 = [-pi/2:0.001:0.001];  % [radians]
theta2 = [0.001:0.001:pi/2];   % [radians]

%%
% Note that if $\theta=0$, will cause problems
% when trying to evaluate the beam pattern as the denominator contains
% the term $\sin(\theta)=\sin(0)=0$, hence $1/\sin(0) =\infty=$|NaN|.

%%
% Define an arbitrary radius
r_radius = 1;               % [m] distance from front of piston

%%
% Define a shorthand variable for $ka$
ka = k_wavenum*a_piston_radius;

%%
% The pressure variation with angle is given by
pressure1 = (1i/2)*rho_dens_air * c_speed_air * U_max_vel    ...
            *(a_piston_radius/r_radius)*(ka)    ...
            *(                                  ...
                2*besselj(1,ka*sin(theta1)) ./  ...
                (ka * sin(theta1))              ...
            );
pressure2 = (1i/2)*rho_dens_air * c_speed_air * U_max_vel    ...
            *(a_piston_radius/r_radius)*(ka)    ...
            *(                                  ...
                2*besselj(1,ka*sin(theta2)) ./  ...
                (ka * sin(theta2))              ...
            );
        
%%        
% Calculate the pressure at $\theta=0$ 
% noting that $[2 J_1(0)/0]=1$,
p_0 = (1i/2)*rho_dens_air * c_speed_air * U_max_vel    ...
            *(a_piston_radius/r_radius)*(ka)    ...
            *(                                  ...
               1                                ...
            );

%%
% Paste all the results together
theta       = [theta1, 0, theta2];
pressure    = [pressure1, p_0, pressure2];
        
%%
% Calculate the beam radiation pattern relative 
% to the pressure at $\theta=0$
b_theta =   20*log10(abs(pressure/20e-6)) -   ...
            20*log10(abs(p_0/20e-6)) ;

%% Calculate the angles where pressure is zero
%
% There are pressure nodes at angles $\theta_m$ given by
% the solution to 
%
% $$ J_1(ka \sin \theta_{m}) = 0 $$
%
% Define $v_{m}=ka \sin \theta_m$
%
% The zeros of the Bessel function can be solved numerically.
%
% We know that $\theta_m \leq \pi/2$, hence the maximum value of $v$ is
% $v_{\rm max} = ka$
%
% Our initial guesses are
%
% $v_{m}=1+\sqrt{2}+(m-1) \pi +1+1^{0.4}$
%
% This can be rearranged for $m$ as
%
% $$m = \frac{ka - 3 -\sqrt{2}}{\pi} +1 $$
%
% In this case
m_max = floor( (ka - 3 -sqrt(2))/(pi)+1 );


%%
% Find the |m_max| zeros
if m_max<1,
    warning('For the parameters of the piston, there are no angles where the pressure is zero.');
else
    for jj=1:m_max,
        % The initial guess x0 for the zeros are
        x0=1+sqrt(2)+(jj-1)*pi+1+1^0.4;
        % Use an inline equation and solve for the zeros
        [v(jj), FVAL, EXITFLAG] = fzero(inline('besselj(1,v)','v'), x0);
        % Check if the values were calculated correctly.
        if (EXITFLAG>0),
            % Found a zero.... good, continue to the next zero
        end;
        if (EXITFLAG<0),
            % Bad, could not solve for the zero... 
            %should end program here so you can figure what is wrong
            error('Unable to solve for the zeros of Bessel function.');
        end;    
     end;    
    %%
    % Calculate the equivalent angles for the pressure zeros
    theta_m=asin(v/ka);             %[radians] 

    % The polar plot in Matlab does not handle -ve values well
    % and will reflect through the origin, rotating by $\pi$.
    % Instead, the beam pattern will be shifted by +50dB,
    %
    result_offset=50;               %[dB] offset for graphs

    % Create vectors for the pressure null lines
    theta_zeros = [theta_m; theta_m];
    p_0_lines = repmat([0;result_offset],1,m_max);

end;


        
%% Plot the results

%Alter the results where the SPL_dipole_norm < -100 to SPL=0
indexes=find(b_theta<-60);
b_theta(indexes)=-60;

p1h=polarplot( theta,b_theta,'k-');

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
set(hA,'RAxisLocation',5);  

% Try to move position of tick labels
rTL=get(hA,'RTickLabels');
% for jj=1:length(rTL)
%     rTL_v2{jj,1}=strjoin({' ',rTL{jj}});
%  end;
% Delete the first label at the origin
rTL_v2=rTL;
rTL_v2{1}='';
set(hA,'RTickLabels',rTL_v2);


% make gridlines more visible
set(hA,'GridAlpha',1.0);

set(hA,'ThetaLim',[-90 90]);
% set(hA,'ThetaZeroLocation','top')
% set(hA,'ThetaDir','clockwise')

return;




%%
% The following code uses the inbuilt Matlab |polar| function to plot the
% directivity pattern, but it does not look very good.
%
% The polar plot in Matlab does not handle -ve values well
% and will reflect through the origin, rotating by $\pi$.
% The workaround is to shift the beam pattern by +50dB.

% result_offset=50;               %[dB] offset for graphs

% Create vectors for the pressure null lines
% theta_zeros = [theta_m; theta_m];
% p_0_lines = repmat([0;result_offset],1,m_max);
% p1=polar(theta,(b_theta)+result_offset);
% set(gca,'Fontsize',16);
% ax=axis;
% ax(1)=0;
% axis(ax);
% % Over-plot the lines showing the angle of the pressure zeros
% hold on
% p2=polar(theta_zeros,p_0_lines,'r--');
% hold off

%%
% The following code will plot using the directivity pattern using a script
% called |dirplot| which can be obtained from the Matlab File Exchange web
% site:
% <http://www.mathworks.com/matlabcentral/fileexchange/1251-dirplot>
p1=Dirplot(theta*180/pi,(b_theta),'-',[0 -50,5]);
set(gca,'Fontsize',16);
% Set font size, line thickness and add graph labels
set(p1(:),'Linewidth',2);
title('Beam Pattern of a Baffled Circular Piston');

if exist('theta_zeros'),
    hold on
    p_0_lines = repmat([-50 ; 0],1,m_max);
    p2=Dirplot(theta_zeros*180/pi,p_0_lines,'r--',[0 -50 5]);
    hold off
    set(p2(:),'Linewidth',2);
    l1=legend('Directivity Pattern','Pressure Nulls');
else
	l1=legend('Directivity Pattern');
end;


grid
return;

