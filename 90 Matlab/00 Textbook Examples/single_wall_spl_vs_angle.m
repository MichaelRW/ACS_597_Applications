%% Plot directivity of Single Wall Incoherent Plane Radiator
% See Bies and Hansen 5th Edition, Section 5.7
% 
% Carl Howard 14 October 2016
%

%% Define properties of air and the source
rho=1.21;
c=343;


%% Define coordinates of wall
H = 4;   % [m] width of wall
L = H;    % [m] length of wall

W=1;

%% plot the SPL in an arc around the wall
radius = 4;
theta = [-85:1:85];     %[degrees]

r=cos(theta*pi/180)*radius;
d=sin(theta*pi/180)*radius;
h=H/2;

alpha=H/L;
beta=h/L;

SPL=zeros(length(theta),1);

for aa=1:length(theta),

    gamma=r(aa)/L;
    delta=d(aa)/L;
        
    p_squared = rho * c * W / (2*pi*H*L) *  ...
        (       ...
        atan( (alpha-beta)*(delta+0.5) /  ...
        (gamma * sqrt( (alpha-beta)^2 + (delta+0.5)^2 + gamma^2 ) ) )...
        + atan( beta*(delta+0.5) / (gamma*sqrt( beta^2 + (delta+0.5)^2 + gamma^2 )  )  )    ...
        - atan( (alpha-beta)*(delta-0.5) / (gamma*sqrt( (alpha-beta)^2 + (delta-0.5)^2 +gamma^2)  )  ) ...
        - atan( beta*(delta-0.5) / (gamma*sqrt(beta^2 + (delta-0.5)^2 +gamma^2 )  ) )       ...
        );
    
    SPL(aa)=10*log10( abs(p_squared) / (20e-6)^2);
    
end;

%%%
% Polar plot of the sound radiation
% This is not very interesting as the directivity looks nearly like a
% monopole, and varies with distance. 
%
% This figure is not going to be included in the book.
figure(1);
% plot(theta,SPL)
% xlabel('Angle')
% ylabel('SPL')

p1h=polarplot(theta*pi/180,SPL,'k');

% Set the page size for the figure
set(1, 'units', 'centimeters', 'pos', [0 0 10 10])

set(p1h(1),'linewidth',2);

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

set(hA,'ThetaLim',[-90 90]);

% Set the limits for the radial axis
%set(hA,'RLim',[-60 0]);    

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


%% Plot a whole contour map
% Try to recreate Fig 8 in the paper by Hohenwarter

% Define the reference intensity level 1 pW/m^2
I_0 = 1e-12;    % [W/m^2]

% Define the power out of the rectangular opening.
% This comes from the caption on the graph in Fig 8 in the paper, which is
% written at $W/(4 \pi H L I_{0})=1$.
W=4*pi*H*L*I_0; % [W]

% Define the positions where the intensity is to be evaluated.
r=[0.1:0.1:30];     % [m] these are the coordinates in the horizontal axis.
d=0;                % [m] d=0 is along the centre of the 
h=[0:0.1:22];       % [m] these are the vertical axis

% Define some non-dimensional terms
alpha=H/L;          
delta=d/L;

% Create an empty matrix to speed up the for loops.
SPL=zeros(length(r),length(h));

for aa=1:length(r),
    for bb=1:length(h),
        
        gamma=r(aa)/L;
        beta=h(bb)/L;
        
        p_squared = rho * c * W / (2*pi*H*L) *  ...
            (       ...
            atan( (alpha-beta)*(delta+0.5) /  ...
            (gamma * sqrt( (alpha-beta)^2 + (delta+0.5)^2 + gamma^2 ) ) )...
            + atan( beta*(delta+0.5) / (gamma*sqrt( beta^2 + (delta+0.5)^2 + gamma^2 )  )  )    ...
            - atan( (alpha-beta)*(delta-0.5) / (gamma*sqrt( (alpha-beta)^2 + (delta-0.5)^2 +gamma^2)  )  ) ...
            - atan( beta*(delta-0.5) / (gamma*sqrt(beta^2 + (delta-0.5)^2 +gamma^2 )  ) )       ...
            );

        %SPL(aa,bb) = 10*log10( p_squared / (20e-6^2)  );
        % To plot the graph in the journal paper, we need to plot 
        % the Sound Intensity Level, so we will tweak the
        % p_squared matrix
        SPL(aa,bb) = 10*log10( abs(p_squared) / (rho*c*I_0)  );

    end;
end;


figure(2);
%[C,h]=contour(r,h,SPL.',[-18:2:4],'ShowText','on');
[C,hC]=contour(r,h,SPL.',[-18:2:4]);

% Set the page size for the figure
figno=gcf;
set(figno, 'units', 'centimeters', 'pos', [0 0 10 10])

% Change the colour of all the contour lines to black
set(hC,'LineColor','k')

hA=gca;
%set(hA,'XTick',[0:3:30]);
set(hA,'XTick',[0:2:30]);
%set(hA,'YTick',[0:2:22]);
set(hA,'YTick',[0:2:20]);

% ax2=gca;
% axislims=axis;
% ax2.YAxis.MinorTick = 'on';
% ax2.YAxis.MinorTickValues = [axislims(3):1:axislims(4)];
% ax2.YMinorGrid = 'on';
% grid
% set(ax2,'GridAlpha',1.0);
% set(ax2,'MinorGridAlpha',1.0);

%axis([0 30 0 22])
axis([0 30 0 20])


fontname ='Times New Roman';
set(hA,'FontName',fontname);
set(hA,'FontSize',9);
set(hA,'defaultaxesfontname',fontname);
set(hA,'defaulttextfontname',fontname);


axis equal


% To put manually labels on the contours, use the following command that
% will have to be tidied later.
%t = clabel(C,hC,'manual','LabelSpacing',200,'FontSize',9)
clabel(C,hC,'FontName',fontname)
clabel(C,hC,'FontSize',9)


xlabel('Distance $r$ [m]','Interpreter','Latex')
ylabel('Distance $h$ [m]','Interpreter','Latex')


%%%
% Plot a horizontal line at the centreline of the radiating area
hold on
h1a=plot([0 30],[H/2 H/2],'k-.');



%return;

%% Double check this result
% There is an Eq. (5) in the journal paper for the Intensity in the plane
% where d=0;

W=4*pi*L*H*I_0; % [W]

n=[0.1:0.1:30];
d=0;
h=[0:0.1:22];

intensity_v2=zeros(length(n),length(h));

for aa=1:length(n),
    for bb=1:length(h),
        
        intensity_v2(aa,bb) = W / (pi*H*L) *  ...
            (       ...
            atan( L / (2*n(aa))  ...
            * (H-h(bb)) / sqrt( (H-h(bb))^2 + (L^2)/4 + n(aa)^2 ) )...
            + atan( L / (2*n(aa))   ...
            * h(bb) / sqrt( h(bb)^2 + (L^2)/4 + n(aa)^2 ) )     ...
            );

    end;
end;

L_intensity_v2 = 10*log10( abs(intensity_v2) / (I_0)  );

% hold on
% [C2,h2]=contour(n,h,L_intensity_v2.',[-18:2:4],'LineStyle','--','color','r','ShowText','on');

% Put in the legend
% l2h=legend('B&H Eq (5.117)','Centreline','Hohenwarter Eq (5)');


%% Plot the Intensity Along the Centreline
% There is an Eq. (6) in the journal paper for the Intensity in the plane
% where d=0;

W=4*pi*L*H*I_0; % [W]


n=[0.1:0.1:30];
d=0;
h=[0];

intensity_v3=zeros(length(n),length(h));

for aa=1:length(n),
        intensity_v3(aa,1) = 2* W / (pi*H*L) *  ...
            (       ...
            atan( L * H/ (2*n(aa))  ...
            * 1 / sqrt( 4 *n(aa)^2 + (L^2) + (H)^2  ) ) ...
            );

end;

L_intensity_v3 = 10*log10( abs(intensity_v3) / (I_0)  );

%%%
% Calculate the result for a standard monopole
I_monopole = W*ones(size(n)) ./ (2*pi*n.^2);
L_monopole = 10*log10(I_monopole/I_0);



h=[0:0.1:22];
index=find(h==H/2);
if isempty(index),
    warning('There isn''t a column at the correct location in the centre.');
    disp('Instead will try to find the closest.');
    index=find(h<=H/2);
    if isempty(index),
        error('Sorry still could not find a location close to centre.');
    end;
    index=index(end);
    disp(['The location that was found was h=',num2str(h(index))]);
end;

    
figure
%h3=semilogx(n,L_monopole,r,SPL(:,index),'.',n,L_intensity_v3);
h3=plot(n,L_monopole,r,SPL(:,index),'.',n,L_intensity_v3);

ax3=gca;
axislims=axis;
ax3.YAxis.MinorTick = 'on';
ax3.YAxis.MinorTickValues = [axislims(3):2:axislims(4)];
ax3.YMinorGrid = 'on';

ax3.XAxis.MinorTick = 'on';
ax3.XAxis.MinorTickValues = [axislims(1):1:axislims(2)];
ax3.XMinorGrid = 'on';
grid
set(ax3,'GridAlpha',1.0);
set(ax3,'MinorGridAlpha',1.0);

xlabel('Distance $r$ [m]','Interpreter','Latex')
ylabel('Sound Intensity Level $L_{I}$ [dB re $10^{-12} \textrm{W} / \textrm{m}^{2}$]','Interpreter','Latex')


% Put in the legend
l3h=legend('Monopole','Hohenwarter Eq (5)','Hohenwarter Eq (6)');


%% Plot response vs Monopole
% Same as above, but only showing Hohenwarter Eq (6), and used to show how
% the reponse varies for change in ratio H / L .

% L=[1,2,3,4,5,10,15,20];
% H=L;

% L=1;
% H=[0.1 4 16 64];

% An alternative is to assume a constant area of H*L=1m^2
% and vary the value of alpha...
%alpha=[ 1 2 4 8 16 32];
alpha=[ 1, 2, 5, 10, 20 ];

H=sqrt(alpha);
L=ones(size(alpha))./H;


%W=4*pi*I_0*ones(1,length(L))./(L.*H); % [W]

% Assume a constant sound power from the radiating area and monopole.
W=I_0*ones(1,length(H));


%n=[0.1:0.5:30];
n=logspace(-2,2,1000);
d=0;
h=[0];

intensity_v4=zeros(length(n),1);

for aa=1:length(n),
    for bb=1:length(H),
        intensity_v4(aa,bb) = 2* W(bb)/ (pi*H(bb)*L(bb)) *  ...
            (       ...
            atan( L(bb) * H(bb) / (2*n(aa))  ...
            * 1 / sqrt( 4 *n(aa)^2 + (L(bb)^2) + (H(bb))^2  ) ) ...
            );
    end;
end;

L_intensity_v4 = 10*log10( abs(intensity_v4) / (I_0)  );

%%%
% Calculate the result for a standard monopole
I_monopole = W.'*ones(1,length(n)) ./ (2*pi*(ones(length(W),1)*n).^2);
L_monopole = 10*log10(I_monopole/I_0);
   
figure
%[C3,h3]=contour(n,h,L_intensity_v2.',[-18:2:4],'LineStyle','--','color','g','ShowText','on');
%ph4=semilogx(n,L_monopole,'k--',n,L_intensity_v4,'k-');
ph4=semilogx(n,L_monopole(1,:),'k--',n,L_intensity_v4,'-');

set(ph4(:),'LineWidth',1);

figno=gcf;
set(figno, 'units', 'centimeters', 'pos', [0 0 10 10])

hA=gca;

fontname ='Times New Roman';
set(hA,'FontName',fontname);
set(hA,'FontSize',9);

xlabel('Distance $r$ [m]','Interpreter','Latex')
ylabel('Sound Intensity Level $L_{I}$ [dB re $10^{-12} \textrm{W} / \textrm{m}^{2}$]','Interpreter','Latex')

axis ([1 100 -60 -5])

alpha=H./L;

%%% Put a legend on the graph
legend_text(1)={'monopole'};
for jj=1:length(alpha), 
    legend_text(jj+1)={['\alpha=',num2str(alpha(jj))]}; 
end;
lh4=legend(ph4,legend_text);


grid
% make gridlines more visible
set(hA,'GridAlpha',1.0);
set(hA,'MinorGridAlpha',1.0);
set(hA,'MinorGridLineStyle','-');

% Put some dot markers on the lines to enable drawing a legend
hold on
p4b=plot(n(50),L_monopole(1,50),'k.',n(50),L_intensity_v4(50,:),'k.');
set(p4b,'MarkerSize',16);

% Put in the legend
%l4h=legend('Monopole','Hohenwarter Eq (5)','Hohenwarter Eq (6)');


%% Plot Sound Intensity Vs Distance, normalised axes

% Same as above, but only showing Hohenwarter Eq (6), and used to show how
% the reponse varies for change in L=1,2,3,4,5.  H=constant

% L=2;
% H=[0.1,0.2,0.5,1.0]*L;

% L=1;
% H=[0.1 4 16 64];

% An alternative is to assume a constant area of H*L=1m^2
% and vary the value of alpha...
% alpha=[ 1 2 4 8 16 32];
alpha=[ 1, 2, 5, 10, 20 ];
H=sqrt(alpha);
L=ones(size(alpha))./H;

%%% 
% Re-calculate the value of alpha In case someone alters the code so 
% that H and L are explicitly defined, then it is necessary to 
% calculate alpha.
% ... but it should not be necessary as it was defined above.
alpha=H./L;


%W=4*pi*I_0*ones(1,length(L))./(L.*H); % [W]

% Assume a constant sound power from the radiating area and monopole.
W=I_0;

%%%
% define the distance from the wall the sound intensity will be calculated.
%n=[0.1:0.5:30];
n=logspace(-2,2,1000);
d=0;    % not used.
h=H/2;  % not used.

%%%
% define an empty matrix for the intensity to speed up the computations.
intensity_v5=zeros(length(n),1);

for aa=1:length(n),
    for bb=1:length(H),
        intensity_v5(aa,bb) = 2* W/ (pi*H(bb)*L(bb)) *  ...
            (       ...
            atan( L(bb) * H(bb) / (2*n(aa))  ...
            * 1 / sqrt( 4 *n(aa)^2 + (L(bb)^2) + (H(bb)^2)  ) ) ...
            )   ...
            ;
        normalised_xaxis(aa,bb)=n(aa)/sqrt(H(bb)*L(bb));
        %%% 
        % Calculate the diffence is the sound intensity (not LEVEL).
        % Not used.
        intensity_v6(aa,bb) = W / (2*pi*n(aa)^2) - intensity_v5(aa,bb) ;
    end;
end;

%%%
% The calculation of L_intensity_v5 needs to have the sound intensity of a
% monopole source subtracted from it when showing on a normalised scale
% relative to a monopole.
L_intensity_v5 = 10*log10( abs(intensity_v5) / (I_0)  );

% whereas L_intensity_v6 is the diffence between linear sound intensities.
% Not used.
L_intensity_v6 = 10*log10( abs(intensity_v6) / (I_0) );

%%%
% Calculate the result for a standard monopole
I_monopole = W*ones(length(n),1) ./ (2*pi*(n.').^2);
L_monopole = 10*log10(I_monopole/I_0);



figure
% ph5=semilogx(   n(1,:)/sqrt(1),L_monopole,'k--',    ...
%     repmat(n.',1,length(H))./repmat(sqrt(H*L),length(n),1),L_intensity_v5,'k-');

% This is the colour line version
% ph5=semilogx(   n(1,:),L_monopole-L_monopole,'k--',    ...
%       normalised_xaxis,L_intensity_v5-repmat(L_monopole,1,length(H)),'-');

% This is the black line version for the text book.
% Note that this plot does NOT have the dashed line for the monopole, as it
% it is fairly obvious that this is the 0 dB line.
ph5=semilogx(   ...
      normalised_xaxis,L_intensity_v5-repmat(L_monopole,1,length(H)),'k-');

% Not used
% ph5=semilogx(  normalised_xaxis,L_intensity_v6,'-');
  

% I think this is what is the original Fig 5.7 
% ph5=semilogx(   n(1,:),L_monopole-L_monopole(1),'k--',    ...
%      n(1,:),L_intensity_v5-L_monopole(1),'k-');

%%%
% define the range of the axes on the graph
axis ([0.1 100 -15 5])


%%% 
% Put a legend on the graph.
% Note that for the textbook version which is in black and white, the line
% for the monopole is not plotted, to reduce the number of lines shown on
% the graph.
% legend_text(1)={'monopole'};
% for jj=1:length(alpha), 
%     legend_text(jj+1)={['\alpha=',num2str(alpha(jj))]}; 
% end;
% lh5=legend(ph5,legend_text);


set(ph5(:),'LineWidth',1);

figno=gcf;
set(figno, 'units', 'centimeters', 'pos', [0 0 10 10])

hA=gca;

fontname ='Times New Roman';
set(hA,'FontName',fontname);
set(hA,'FontSize',9);

xlabel('Normalised Distance $r/\sqrt{HL}$','Interpreter','Latex')
%xlabel('Distance $r$','Interpreter','Latex')
%ylabel('Sound Intensity Level $L_{I}$ [dB re $10^{-12} \textrm{W} / \textrm{m}^{2}$]','Interpreter','Latex')
ylabel('Relative Level Compared to Monopole [dB]','Interpreter','Latex')

grid
% make gridlines more visible
set(hA,'GridAlpha',1.0);
set(hA,'MinorGridAlpha',1.0);
set(hA,'MinorGridLineStyle','-');


axislims=axis;
hA.YAxis.MinorTick = 'on';
hA.YAxis.MinorTickValues = [axislims(3):1:axislims(4)];
hA.YAxis.TickLength=[0.0300 0.0250];
hA.YMinorGrid = 'off';


% Put some dot markers on the lines to enable drawing a legend
%hold on
% p5a=plot(n(50),L_monopole(1,50),'k.',n(50),L_intensity_v5(50,:),'k.');
% set(p5a,'MarkerSize',10);

% Put in the legend
%l4h=legend('Monopole','Hohenwarter Eq (5)','Hohenwarter Eq (6)');
