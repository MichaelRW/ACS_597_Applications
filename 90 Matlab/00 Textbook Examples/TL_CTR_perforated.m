%% Transmission Loss of a Concentric Tube Resonator Silencer
% This script is used to calculate the transmission loss of a 
% Concentric Tube Resonator (CTR) silencer comprising a section of
% perforated tube in the middle of the expansion chamber section.

% Keywords: 
% pipe ; duct ; 4 pole ;  four pole ; transmission line ;
% silencer; muffler; expansion chamber
%
% Carl Howard 3 November 2016

%% Define the properties of air

% The NCEJ paper says c= 480 m/s at T=300C on p255
% c_0 = 480;       %[m/s] speed of sound
% rho_0 = 1.17;     %[kg/m^3] density of air

% These are the default values used by ANSYS
% c_0 = 343.24;       %[m/s] speed of sound
% rho_0 = 1.2041;     %[kg/m^3] density of air

% The are the values in teh NCEJ paper on p254 at T=28.8C c=348.34 m/s
c_0 = 348.34;       %[m/s] speed of sound
rho_0 = 1.17;     %[kg/m^3] density of air


%% Define properties of gas flow
U_1 = 0;            % [m/s] velocity of flow in Duct 1 (main duct)
U_2 = 0;            % [m/s] velocity of flow in Duct 2 (outside perf tube)
M = 0;              % [no units] Mach number, no flow for this example

%% Define the properties of the main duct
% for this example, it is assumed that the cross-sectional areas of
% the upstream and downstream ducts are the same.
diam_main_duct = 0.050;%[m] diam of main duct
radius_main_duct = diam_main_duct/2;         %[m] radius of main duct

S_duct = pi*(radius_main_duct^2);  %[m^2] area of main duct

L_upstream=0.2;         % [m] length of upstream duct before silencer
L_downstream=0.2;       % [m] length of downstream duct after silencer


%% Define the properties of the expansion tube
% Beranek and Ver (2006) Fig 9.10(c) shows the TL evaluated at 
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

% N_ratio_areas = 16;                  %[] ratio of areas
% S_expansion = N_ratio_areas * S_duct;%[m^2] area of expansion chamber
% radius_expansion_chamber = sqrt(S_expansion / pi);  % [m] radius of expansion chamber
% diam_expansion_chamber = 2*radius_expansion_chamber; % [m] diam of expansion chamber

diam_expansion_chamber =3.0*diam_main_duct; % [m] diam of expansion chamber
radius_expansion_chamber = diam_expansion_chamber/2;  % [m] radius of expansion chamber
S_expansion = pi*radius_expansion_chamber^2;%[m^2] area of expansion chamber

L_expansion = 0.401;                   %[m] length of the expansion chamber

% GOod values are 0.399, 0.0846, 0.185; that cause it to match with
% ansys_CTR_4.
% Another option is 0.401; 0.0849; 0.186

% The dimensions in the ANSYS model are:
% 0.401, 0.0825, 0.183, 
% Eq 11 and 13 result in 13.8 mm and 18 mm

% The dimensions in the ANSYS model for Fig 10 are:
% 0.401, 0.0891, 0.1894, 

N_ratio_areas = S_expansion / S_duct;

% The length of the tube in the downstream contraction section is
%L_ext_tube_contraction = L_expansion /2;
%L_ext_tube_contraction = L_expansion /4;
%L_ext_tube_contraction= 0.0899;
%L_ext_tube_contraction= 0.0825;

% excellent agreement with ANSYS
L_ext_tube_contraction= L_expansion/4 - 0.0130;

% The length of the contraction section for Fig 10 of the NCEJ paper for
% the 1-D TMMP (4-pole) model is L_b = 0.0891 m
%L_ext_tube_contraction= 0.0891;


% The length of the tube in the upstream expansino section is
%L_ext_tube_expansion = L_expansion /4;
%L_ext_tube_expansion = L_expansion /2;

%L_ext_tube_expansion=0.195;
%L_ext_tube_expansion=0.1825

% excellent agreement with ANSYS
L_ext_tube_expansion=L_expansion/2 - 0.0130; 

% The length of the expansion section for Fig 10 of the NCEJ paper for
% the 1-D TMMP (4-pole) model is L_a = 0.1894 m
%L_ext_tube_expansion=0.1894;


% The length of the straight section between the upstream expansion with
% the extension tube, and the downstream contraction that has an extension
% tube is: 
%
% $$ L_{\textrm{straight}} = L_{\textrm{expansion}} - L_{6} - L_{2} $$
%
% (Note: unfortunately the notation for $L_2$ is confusing because
% in the textbook $L_2$ corresponds to the length of the extension tube in
% the expansion / contraction segements, whereas in this Matlab code the
% lengths of the extension tubes are defined as $L_1$ for the downstream
% contraction section, and $L_3$ for the upstream expansion section.)
L_CTR_straight = L_expansion - L_ext_tube_contraction  ...
                    - L_ext_tube_expansion;
                
% Do a sanity check to make sure that the length of the straight section is
% greater than 0
if (~(L_CTR_straight>0)),
    error('The length of the pipe in expansion section is not positive. Check dimensions');
end;


%%%
% Normally one would define the diameter of the expansion chamber duct
% and then calculate the area as follows:
%
%   a_expansion = 0.1;               %[m] radius of the expansion chamber
%   S_expansion = pi*a_expansion^2;  %[m^2] area of expansion chamber

%% Define the analysis frequency range
% Beranek and Ver (2006) Fig 9.10(c) has the x-axis range from 
%
% $$(kL / \pi) = 0 \cdots 4$$ 
%
% Hence this determines the analysis frequency range.
k_max = 5 *pi / L_expansion;    %[] maximum wavenumber of analysis
f_max = k_max * c_0 / (2*pi);   %[Hz] maximum analysis frequency
freq_spacing=0.1;               %[Hz] frequency increment
f = [1:freq_spacing:f_max];     %[Hz] Excitation frequency 
omega= 2*pi*f;                  %[rad/s] circular frequency
k_0 = omega / c_0;              % vector containing wavenumbers
k_c = k_0 / (1-M^2);            % wave number adjusted for Mach number


%% Outlet straight section 
% The transmission matrix for the outlet duct straight segment is
T_1_11    = cos(k_c * L_downstream);
T_1_12    = 1i*(c_0/S_duct) * sin(k_c * L_downstream);
T_1_21    = 1i*(S_duct/c_0) * sin(k_c * L_downstream);
T_1_22    = cos(k_c * L_downstream);
%T_1     = [ T_1_11    T_1_12    ;   T_1_21    T_1_22];

%% Extended Tube Inside Downstream Contraction 
% The transmission matrix for the tube inside the downstream contraction is
T_2_11    = cos(k_c * L_ext_tube_contraction);
T_2_12    = 1i*(c_0/S_duct) * sin(k_c * L_ext_tube_contraction);
T_2_21    = 1i*(S_duct/c_0) * sin(k_c * L_ext_tube_contraction);
T_2_22    = cos(k_c * L_ext_tube_contraction);
%T_2     = [ T_2_11    T_2_12    ;   T_2_21    T_2_22];


%% Perforated tube section

% From Mechel Formula of Acoustics
% Chapter K: Muffler Acoustics, p 807

%sigma_porosity = 0.196;    % porosity

% The NCEJ paper Fig 10 has sigma = 0.27;
sigma_porosity = 0.27;    % [ratio] porosity of perforations

% The NCEJ paper Fig 10  has these values
d_h = 0.003;    % [m] hole diameter
t_h = 0.002;    % [m] thickness of the perforated tube

d_1 = diam_main_duct; % [m] diameter of the perforated tube
d_2 = diam_expansion_chamber; % [m] diameter of expansion chamber

% Length of annular duct from upstream end of expansion chamber to 
% start of perforated tube, which is also the same as the 
% length of the extended tube in the upstream end where 
% there is an expansion.
L_a = L_ext_tube_expansion;   % [m] 

% Length of annular duct from downstream end of end of perforated tube 
% to the downstream end of the expansion chamber, which is also the same 
% as the length of the extended tube in the downstream end where 
% there is a contraction.
L_b = L_ext_tube_contraction;

Mach = M;           % [no units] Mach number through peforated tube

% From Munjal 2014, p121
%k_0 = omega / c_0;
M_1 = U_1 / c_0;    % mach number in Duct 1 (main duct)
M_2 = U_2 / c_0;    % mach number in Duct 2 (in the expansion chamber M_2=0)

% For perforates in stationary media, 
% Eq. 6 Mechel 2008, Chapter K Munjal p808
zeta = (0.006 + 1i*k_0*(t_h+0.75*d_h) )/sigma_porosity;

% For perforates with grazing flow, Eq. 7
% Eq. 7 Mechel 2008, Chapter K Munjal p808
%zeta = (7.337e-3 * (1+72.23*Mach) + 1i*2.2245e-5*(1+51*t_h)*(1+204*d_h)*freq)/sigma;

% Define the wave numbers
% Mechel 2002 Formulas of Acoustics, p807 Eq. 4
% Munjal 2014 , p121-122
k_a = sqrt(k_0.^2 - 1i*4*k_0 ./ (d_1 * zeta));
k_b = sqrt(k_0.^2 - 1i*4*k_0*d_1 ./ ( (d_2^2-d_1^2)*zeta) );


% Eqs. 3a, 3b Mechel 2008, Chapter K Munjal p807
% Munjal 2014 , p121
alpha_1 = -1i*M_1 / (1-M_1^2) * (k_a.^2 + k_0.^2)./k_0;
alpha_2 = k_a.^2 / (1-M_1^2);
alpha_3 = 1i * M_1 / (1-M_1^2) * (k_a.^2 - k_0.^2)./k_0;
alpha_4 = -(k_a.^2 - k_0.^2) / (1-M_1^2);
alpha_5 = 1i*M_2 / (1-M_2^2) * (k_b.^2 - k_0.^2)./k_0;
alpha_6 = -(k_b.^2 - k_0.^2) / (1-M_2^2);
alpha_7 = -1i*M_2 / (1-M_2^2) * (k_b.^2 + k_0.^2)./k_0;
alpha_8 = k_b.^2 / (1-M_2^2);




%% Extended Tube Inside Upstream Expansion
% The transmission matrix for the tube inside the upstream expansion is
T_4_11    = cos(k_c * L_ext_tube_expansion);
T_4_12    = 1i*(c_0/S_duct) * sin(k_c * L_ext_tube_expansion);
T_4_21    = 1i*(S_duct/c_0) * sin(k_c * L_ext_tube_expansion);
T_4_22    = cos(k_c * L_ext_tube_expansion);
%T_4     = [ T_4_11    T_4_12    ;   T_4_21    T_4_22];

%% Upstream straight duct segement at inlet
% The transmission matrix for the upstream inlet straight segment  is
T_5_11    = cos(k_c * L_upstream);
T_5_12    = 1i*(c_0/S_duct) * sin(k_c * L_upstream);
T_5_21    = 1i*(S_duct/c_0) * sin(k_c * L_upstream);
T_5_22    = cos(k_c * L_upstream);

%T_5     = [ T_5_11    T_5_12    ;   T_5_21    T_5_22];


%% Calculate the transmission loss for the concentric tube resonator silencer (CTR)
% Combine all transmission matrices and then calculate the 
% transmission loss $\textrm{TL}$ using 
% Beranek and Ver (2006), Eq. (9.11) p286 as
% 
% $$ \textrm{TL} = 10 \log_{10} 
%   \left[ \left( \frac{1+M_{n}^{2}}{1+M_{1}^{2}} \right)^2 
%   \left( \frac{Y_1}{4 \, Y_n } \right)    
%   \left| T_{11} + \frac{T_{12}}{Y_{1}} 
%           + Y_{n} T_{21} + \frac{Y_n T_{22}}{Y_1}
%   \right|^{2} 
%   \right]     $$
%
%
% Define some constants outside of the |for| loop
M_n = M;        % Mach number at the upstream end of the muffler
M_1 = M;        % Mach number at the downstream end of the muffler

% We are assuming the inlet diameter is the same as the outlet diameter
Y_n = c_0 / S_duct;
Y_1 = c_0 / S_duct;

% Define an empty matrix for the transmission loss of the CTR
% to speed up the |for| loop
TL_CTR = zeros(size(k_0));

%%%
% Now execute the |for| loop. The elements of each transmission matrix were
% calculated as vectors, and now need to be combined to evaluate the total
% transmission matrix as:
%
% $$ \mathbf{T}_{\textrm{total}} = 
%               \mathbf{T}_3 \mathbf{T}_2 \mathbf{T}_1
% $$
%


for jj=1:length(k_0),

    % The following has to be put into a FOR loop
    alpha_matrix = [ -alpha_1(jj)   -alpha_3(jj)    -alpha_2(jj)    -alpha_4(jj)    ;
                     -alpha_5(jj)   -alpha_7(jj)    -alpha_6(jj)    -alpha_8(jj)    ;
                     1          0           0           0       ;
                     0          1           0           0       ];



    % Calculate the eignvectors and eigenvalues of the alpha matrix.
    % Note that the eigenvectors are in columns, and the eigenvalues are
    % actualy returned in a diagonal matrix.
    [psi_eigenvectors,beta_eigenvalues] = eig(alpha_matrix);

    % Extract the entries along the diagonal of the eigenvalue matrix
    beta_eigenvalues = diag(beta_eigenvalues);

    % We will rescale the eigenvectors, which is done in Munjal's textbook
    psi_eigenvectors_norm=psi_eigenvectors./repmat(psi_eigenvectors(1,:),4,1);
    
    % From Munjal 2014, page 124 and Mechel p807
    % The description in Munjal's textbook for psi can cause confusion as
    % typically one would consider $\psi_{1,i}$ as the 1st eigenvector, and
    % each of the i'th elements.... but no! 
    % $\psi_{1,i}$ is meant to be interpretted as the first of all 
    % eigenvectors!
    % For interested readers, this can be confirmed by using the
    % expressions on Munjal (2014) p122 for the relationships between 
    % $\psi_{1,i}$, $\psi_{2,i}$, $\psi_{3,i}$, $\psi_{4,i}$ and confirm
    % that the rescaled eigenvectors |psi_eigenvectors_norm| conforms to
    % these expressions. For example, the interested reader could put a
    % dbstop into this script at this point HERE and confirm that
    % |psi_eigenvectors_norm(4,:)| 
    % is the same as |psi_eigenvectors_norm(2,:).*psi_eigenvectors_norm(3,:)|
    % and also that |psi_eigenvectors_norm(4,:)| 
    % is the same as |psi_eigenvectors_norm(2,:)./beta_eigenvalues.'|
    
    % First calculate A_start_perf
    z=0;
    A_1_row = psi_eigenvectors_norm(3,:).*exp(beta_eigenvalues.'*z);
    A_2_row = psi_eigenvectors_norm(4,:).*exp(beta_eigenvalues.'*z);

    % There are likely errors in the formulas for A_3,i and A_4,i in the
    % publications by Munjal in the 2014 book and the Mechel chapter as the
    % formulas are inconsistent!!
    % There is possibly a missing psi in the Mechel / Munjal 2008 for A_3_row
    % and it unclear if the denominator in Munjal 2nd ed 2008 for
    % A_3_row should be M_1 instead of as written as M_2.
    % For the example M=0 so it is irrelevant.
    A_3_row = -1*psi_eigenvectors_norm(1,:).*exp(beta_eigenvalues.'*z) ./ (1i*k_0(jj) + M_2 * beta_eigenvalues.');

    %A_3_row = -1*exp(beta_eigenvalues.'*z) ./ (1i*k_0(jj) + M_1 * beta_eigenvalues.');

    A_4_row = -1*psi_eigenvectors_norm(2,:).*exp(beta_eigenvalues.'*z)  ./ (1i*k_0(jj) + M_2 * beta_eigenvalues.');

    A_start_perf = [    A_1_row;    ...
                        A_2_row;    ...
                        A_3_row;    ...
                        A_4_row     ];


    % From Munjal 2014, page 124 and Mechel p807
    % First calculate A_start_perf
    z=L_CTR_straight;
    A_1_row = psi_eigenvectors_norm(3,:).*exp(beta_eigenvalues.'*z);
    A_2_row = psi_eigenvectors_norm(4,:).*exp(beta_eigenvalues.'*z);

    % There are likely errors in the formulas for A_3,i and A_4,i in the
    % publications by Munjal in the 2014 book and the Mechel chapter as the
    % formulas are inconsistent!!
    A_3_row = -1*psi_eigenvectors_norm(1,:).*exp(beta_eigenvalues.'*z) ./ (1i*k_0(jj) + M_2 * beta_eigenvalues.');

    %A_3_row = -1*exp(beta_eigenvalues.'*z) ./ (1i*k_0(jj) + M_1 * beta_eigenvalues.');

    A_4_row = -1*psi_eigenvectors_norm(2,:).*exp(beta_eigenvalues.'*z)  ./ (1i*k_0(jj) + M_2 * beta_eigenvalues.');

    A_end_perf  = [     A_1_row;    ...
                        A_2_row;    ...
                        A_3_row;    ...
                        A_4_row     ];
                    

    % Munjal 2014 , p123, Eq 3.117
    T_CTR = A_start_perf * inv(A_end_perf);
                    
    T_11 = T_CTR(1,1);
    T_12 = T_CTR(1,2);
    T_13 = T_CTR(1,3);
    T_14 = T_CTR(1,4);

    T_21 = T_CTR(2,1);
    T_22 = T_CTR(2,2);
    T_23 = T_CTR(2,3);
    T_24 = T_CTR(2,4);

    T_31 = T_CTR(3,1);
    T_32 = T_CTR(3,2);
    T_33 = T_CTR(3,3);
    T_34 = T_CTR(3,4);

    T_41 = T_CTR(4,1);
    T_42 = T_CTR(4,2);
    T_43 = T_CTR(4,3);
    T_44 = T_CTR(4,4);

    X_1 = -1i * tan(k_0(jj) * L_a);
    X_2 = 1i * tan(k_0(jj) * L_b);

    F_1 = T_42 + X_2 * T_44 - X_1 * (T_22 + X_2 * T_24);
    % The eqn in Munjal 1987 textbook p136 is different,
    % but I think it is wrong... the first term is T_32
    % whereas other publications have T_42
    %F_1 = T_32 + X_2 * T_44 - X_1 * (T_22 + X_2 * T_24);

    A_1 = (X_1 * T_21 - T_41) / F_1;
    A_2 = T_12 + X_2 * T_14;
    B_1 = (X_1 * T_23 - T_43) / F_1;
    B_2 = T_32 + X_2 * T_34;


    T_3_11 = T_11 + A_1 * A_2;
    T_3_12 = Y_1*(T_13 + B_1 * A_2);
    T_3_21 = (T_31 + A_1 * B_2) / Y_n;
    T_3_22 = Y_1 / Y_n * (T_33 + B_1 * B_2);



    T_1 = [ T_1_11(jj)      T_1_12(jj);     ...
            T_1_21(jj)      T_1_22(jj)  ];

    T_2 = [ T_2_11(jj)      T_2_12(jj);     ...
            T_2_21(jj)      T_2_22(jj)  ];

    T_3     = [ T_3_11    T_3_12    ;       ...
                T_3_21    T_3_22];

    T_3 = [ 1           M_1 *Y_n; ...
            M_1/Y_n     1       ]   *   ...
            T_3 *                       ...
            inv([1      M_2 * Y_1       ;   ...
                M_2/Y_1     1           ]);
            
    T_4 = [ T_4_11(jj)      T_4_12(jj);
            T_4_21(jj)      T_4_22(jj)  ];

    T_5 = [ T_5_11(jj)      T_5_12(jj);
            T_5_21(jj)      T_5_22(jj)  ];

        
    T_total = T_5 * T_4 * T_3 * T_2 * T_1;
    
% Extract the elements from the matrix
    T_11 = T_total(1,1);
    T_12 = T_total(1,2);
    T_21 = T_total(2,1);
    T_22 = T_total(2,2);

% Calculate the transmission loss
    TL_CTR(jj) = 10*log10(    ...
               ( (1+M_n^2)/(1+M_1^2) )^2  * (Y_1 / (4 * Y_n) ) * ...
               abs( T_11 + (T_12/Y_1) + (Y_n*T_21) + (Y_n*T_22/Y_1) )^2 ...
               );
    
    
end;



%% Calculate Transmission Loss for Simple Expansion Chamber silencer (SEC)
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
% [ see also Beranek and Ver (2006),
% "Noise and Vibration Control", Eq. (9.25), p293]
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
TL_SEC = 10*log10( 1 + 0.25*(1/N_ratio_areas - N_ratio_areas)^2    ...
                    *( sin(k_0*L_expansion) ).^2   );


%% Plot the transmission loss for the CTR and SEC
% x_axis=k_0*L_expansion/pi;          % normalised frequency
% p1h=plot(x_axis,TL_CTR,'k-',x_axis(1:20:end),TL_SEC(1:20:end),'k-.');
% grid
% 
% % Set the page size for the figure
% set(1, 'units', 'centimeters', 'pos', [0 0 10 10])
% 
% % Set some of the font properties
% hA=gca;
% 
% fontname ='Times New Roman';
% set(hA,'defaultaxesfontname',fontname);
% set(hA,'defaulttextfontname',fontname);
% 
% % make gridlines more visible
% set(hA,'GridAlpha',0.8);       
% 
% % Turn on the minor Tick marks
% set(hA,'XMinorTick','on','YMinorTick','on')
% hA.YAxis.MinorTickValues = [0:2:60];
% hA.XAxis.MinorTickValues = [0:0.1:4];
% set(hA,'ticklength',[0.0200    0.0250])
% 
% set(p1h,'linewidth',2);
% l1=legend(  ['CTR: $N=',num2str(N_ratio_areas),'$'],   ...
%             ['SEC: $N=',num2str(N_ratio_areas),'$']);
% set(l1,'Interpreter','Latex','FontSize',9,'FontName','Times New Roman');        
% axis([0 4 0 60]);
% xlabel('Normalised Frequency $kL/\pi$','Interpreter','latex','FontName','Times New Roman');
% ylabel('Transmission Loss [dB]','Interpreter','latex','FontName','Times New Roman');
% set(hA,'FontUnits','points','FontWeight','normal', ...
%     'FontSize',9,'FontName','Times New Roman');
% axis square
% %title('TL of an Expansion Chamber Silencer');
% 
% 


%% Plot the results
figure
%plot(f,TL_CTR,ansys_CTR4(:,1),ansys_CTR4(:,2));
%p1h=plot(f,TL_CTR,'k-',ansys_CTR_Fig10(:,1),ansys_CTR_Fig10(:,2),'k.',f,TL_SEC,'k--');

xaxis= (2*pi*f/c_0)*L_expansion/pi;
p1h=plot(xaxis,TL_CTR,'k-',     ...
        (2*ansys_CTR_Fig10(:,1))/c_0*L_expansion,ansys_CTR_Fig10(:,2),'k.',     ...
        xaxis,TL_SEC,'k--');


grid

figno=gcf;
% Set the page size for the figure
set(figno, 'units', 'centimeters', 'pos', [0 0 10 8])

% Set some of the font properties
hA=gca;

fontname ='Times New Roman';
set(hA,'defaultaxesfontname',fontname);
set(hA,'defaulttextfontname',fontname);

% make gridlines more visible
set(hA,'GridAlpha',0.8);       

% Turn on the minor Tick marks
%set(hA,'XMinorTick','on','YMinorTick','on')
% hA.YAxis.MinorTickValues = [0:2:60];
% hA.XAxis.MinorTickValues = [0:0.1:4];
%set(hA,'ticklength',[0.0200    0.0250])

set(p1h,'linewidth',1);
% l1h=legend(  ['CTR: $N=',num2str(N_ratio_areas),'$'],   ...
%             ['SEC: $N=',num2str(N_ratio_areas),'$']);
l1h=legend('Theory CTR','ANSYS CTR','Theory SEC');
set(l1h,'Interpreter','Latex','FontSize',9,'FontName','Times New Roman');        
%axis([0 2000 0 45]);

axis([0 5 0 45]);

xlabel('Normalised frequency, $kL/\pi$','Interpreter','latex','FontName','Times New Roman');
%xlabel('Frequency (Hz)','Interpreter','latex','FontName','Times New Roman');
ylabel('Transmission Loss (dB)','Interpreter','latex','FontName','Times New Roman');
set(hA,'FontUnits','points','FontWeight','normal', ...
    'FontSize',9,'FontName','Times New Roman');
