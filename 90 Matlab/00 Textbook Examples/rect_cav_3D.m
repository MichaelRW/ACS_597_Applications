function [pressure]=rect_cav_3D;
% This will calculate the pressure response inside a cavity for active
% control



% Put in the values for a rectangular cavity.

% Dimensions of the cavity
cavity.len_x = 2.5;    %m
cavity.len_y = 3.0;    %m
cavity.len_z = 5.0;    %m
cavity.volume = cavity.len_x .* cavity.len_y .* cavity.len_z;   %m^3

% Parameters of the fluid
fluid.dens      = 1.29;     % kg/m^3
fluid.speed     = 343;      % m/s 

%cavity.damping = 0.01;
cavity.damping = 0.00;


% Make sure that 0<source_z<cavity.len_z
% if (source_z>cavity.len_z)|(source_z<0),
%     error ('Source_z is outside dimensions of cavity');
% end;    

% Location of the source
source.loc_x = 0.5;   %m
source.loc_y = 0.5;   %m
source.loc_z = 0.5;   %m

% amplitude of the source
source.amplitude = 0.01;   % m^3/s

% 
% Location of the microphones in the cavity
num_mics = 1;   % The number of error microphones

mic.loc_x(1) = 2.0;       %m
mic.loc_y(1) = 1.5;       %m
mic.loc_z(1) = 2.5;       %m

% mic.loc_x(2) = 0.100;       %m
% mic.loc_y(2) = 0.10;       %m
% mic.loc_z(2) = 1.51;       %m


speedair  = fluid.speed;
densair  = fluid.dens;

% The number of acoustic modes to consider for the cavity.
ncavmodes   =   2000;



%============================================================
% natural frequencies of the cavity are:
%============================================================
% create a matrix containing all the mode index permutations
fprintf('===Start calculation of the resonance frequencies of the cavity\n');

% n_max^3 the number of modes to calculate the resonance frequencies
% which will be sorted and truncated to num_modes
n_max = 30;

mode_index=[];
for ii=0:n_max
    for jj = 0:n_max
        for kk = 0:n_max
            mode_index = [ mode_index ; ii jj kk ];
        end;
    end;
end;

% calculate the resonance frequencies of the room
res_freq = sqrt(    (mode_index(:,1) / cavity.len_x).^2 +   ...
                    (mode_index(:,2) / cavity.len_y).^2 +   ...    
                    (mode_index(:,3) / cavity.len_z).^2 );

res_freq = res_freq .* fluid.speed ./ 2;

% sort the resonance frequencies
temp=sortrows([res_freq mode_index],1);
res_freq = temp(:,1);
mode_index = temp(:,2:4);

% use only first 1:num_modes
res_freq = res_freq(1:ncavmodes);
mode_index = mode_index(1:ncavmodes,:);


% calculate the modal volume
% see Beranek and Ver book P162
modal_volume = zeros(size(mode_index));
index2change = find(mode_index(:,1)==0);
modal_volume(index2change,1) = 1;
index2change = find(mode_index(:,1)>0);
modal_volume(index2change,1) = 0.5;
index2change = find(mode_index(:,2)==0);
modal_volume(index2change,2) = 1;
index2change = find(mode_index(:,2)>0);
modal_volume(index2change,2) = 0.5;
index2change = find(mode_index(:,3)==0);
modal_volume(index2change,3) = 1;
index2change = find(mode_index(:,3)>0);
modal_volume(index2change,3) = 0.5;

modalvolume = cavity.volume .* modal_volume(:,1) .* modal_volume(:,2) .* modal_volume(:,3);
modalvolume = modalvolume.';

cavfreq = res_freq.';
cav_mode_index = mode_index;

% evaluate the mode shape at the measurement location
% phi_eval =  cos(cav_mode_index(:,1).*pi.*mic.loc_x./ cavity.len_x) .*   ...
%             cos(cav_mode_index(:,2).*pi.*mic.loc_y./ cavity.len_y) .*   ...
%             cos(cav_mode_index(:,3).*pi.*mic.loc_z./ cavity.len_z);
% phi_eval = phi_eval(2:ncavmodes,:); % have to drop the first mode        

phi_error_mic = zeros(ncavmodes,num_mics);
for ii=1:num_mics,
	phi_error_mic(:,ii) =  ...
                cos(cav_mode_index(:,1).*pi.*mic.loc_x(ii)./ cavity.len_x) .*   ...
                cos(cav_mode_index(:,2).*pi.*mic.loc_y(ii)./ cavity.len_y) .*   ...
                cos(cav_mode_index(:,3).*pi.*mic.loc_z(ii)./ cavity.len_z);
end; %end of ii
phi_error_mic=phi_error_mic.';
%check the size of the phi_error_mic
if size(phi_error_mic)~=[num_mics,ncavmodes],
    error('Error in the size of the phi_error_mic matrix');
end;


phi_source =    cos(cav_mode_index(:,1)*pi*source.loc_x/cavity.len_x) .*   ...
				cos(cav_mode_index(:,2)*pi*source.loc_y/cavity.len_y) .*   ...
				cos(cav_mode_index(:,3)*pi*source.loc_z/cavity.len_z);
phi_source = phi_source.';


%============================================================
fprintf('===End of calculation of the resonance frequencies of the cavity \n');
%============================================================



% Do the big frequency loop
freq_vect = [1:200];
end_freq = freq_vect(length(freq_vect));


% Calculate the primary acoustic loading
%primary_load = densair*speedair^2*(phi_source./modalvolume)*source.amplitude;
% Note that the equation for the RHS is
%   RHS = rho * c^2 / (modal_vol) * (vol_accel)
% where the modal_vol and the vol_accel are modal terms 
% we shall assume that the units are volume accel so that we don't have to
% multiply this term by the frequency when we get to the big frequency loop

% For some reason unresolved at the moment, the density has to be removed
primary_load = speedair^2*(phi_source./modalvolume)*source.amplitude;



% Define empty matrices to speed up the calculations
J_none = zeros(num_mics,length(freq_vect));

warning off

for ii=1:length(freq_vect),
    %disp(['Freq ',num2str(freq_vect(ii)),' of ',num2str(end_freq)]);
    omega = 2*pi*freq_vect(ii);

% For the comparison with FastBEM, we need to alter the source loading from
% a volume acceleration (which is used in Ansys) to a volume velocity
% source. To convert from vol accel to vol velocity, we just divide by the
% frequency vector
primary_load2 = primary_load ./ (1i*omega*ones(size(primary_load)));
    
    

    % This is the active control calculation equations
    % See Nelson and Elliot Appendix A.5
    % equations A5.21
    %
    %   The error is given by the equation:
    %   e = d + Cx
%     small_d = 1i*omega* phi_error_mic *   ...
%                 (primary_load ./ ( (2*pi*cavfreq).^2 +2*1i*cavity.damping*(2*pi*cavfreq)*omega - omega^2 )).';
    
    small_d =  1i*omega* phi_error_mic *   ...
                (primary_load2 ./ ( (2*pi*cavfreq).^2 +2*1i*cavity.damping*(2*pi*cavfreq)*omega - omega^2 )).';
   
    % The value of the error without active control is
    J_none(1:num_mics,ii) = small_d;
    
end; % of the big frequency loop

warning on

disp('End of calculations')

pressure=J_none;



%% Plot the results compared with FastBEM
% The matrix fastbem2 should have columns [ freq SPL ] from the FastBEM
% % analysis
% 
% theory=pressure;
% 
% p1h=plot(freq_vect,20*log10(abs(theory)/20e-6/sqrt(2)),'k-',fastbem2(:,1), fastbem2(:,2),'k.');
% 
% % Set the page size for the figure
set(1, 'units', 'centimeters', 'pos', [0 0 10 10])

% Set some of the font properties
hA=gca;

fontname ='Times New Roman';
set(hA,'defaultaxesfontname',fontname);
set(hA,'defaulttextfontname',fontname);

% make gridlines more visible
set(hA,'GridAlpha',0.8);       

% Turn on the minor Tick marks
set(hA,'XMinorTick','on','YMinorTick','on')
hA.YAxis.MinorTickValues = [0:10:100];
hA.XAxis.MinorTickValues = [0:10:200];
set(hA,'ticklength',[0.0200    0.0250])

set(p1h,'linewidth',1);
l1=legend('Matlab','BEM','Ansys');
set(l1,'Interpreter','Latex','FontSize',9,'FontName','Times New Roman');        
axis([0 200 0 100]);
xlabel('Frequency [Hz]','Interpreter','latex','FontName','Times New Roman');
ylabel('Sound Pressure Level [dB re 20 $\mu$Pa]','Interpreter','latex','FontName','Times New Roman');
set(hA,'FontUnits','points','FontWeight','normal', ...
    'FontSize',9,'FontName','Times New Roman');
%axis square


