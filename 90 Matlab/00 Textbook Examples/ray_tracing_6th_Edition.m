% This code includes two methods used for acoustic ray tracing in 2D. The code is designed for atmospheric acoustics such that the source and 
% receiver are always above ground level (z = 0) and the acoustic rays reflect at the ground surface (z = 0). However, it can be easily adapted 
% for underwater acoustics with some small tweaks.

% Method 2 (presented first) is used in the Bellhop code. Here, the Euler numerical discretisation method is used in place of the Runge-Kutta 
% one for simplicity. Method 1 is taught at TU Delft in the "Acoustic Remote Sensing and Sea Floor Mapping Course" in Lecture 3.

% More information about the procedure, including all relevant equations needed to implement Methods 1 and 2 is provided in the textbook, 
% “Engineering Noise Control: Theory and Practice Sixth Edition” Sections 5.3.4.4 and 5.3.4.5. 

clear all
close all

% wind speed profile - here we assume a linear profile to calculate the
% sound speed gradient
h_s = 100;
h_e = 5000;
c_s = 400;
c_e = 500;

psi0_deg = -5; % ray launch angle
psi0 = deg2rad(psi0_deg);
% deltar = 25; % step size along ray path

% initial conditions
x0 = 0;
h0 = h_s; % source height
c0 = ((c_s - c_e)/(h_s - h_e))*(h0 - h_e) + c_e; % sound speed at source height
grad_x0 = 0;
grad_z0 = (c_e - c_s)/(h_e - h_s); % sound speed gradient
xi0 = cos(psi0)/c0;
zeta0 = sin(psi0)/c0;

x(1) = x0;
z(1) = h0;
c(1) = c0;
s(1) = 0;
t(1) = 0;
grad_x(1) = grad_x0;
grad_z(1) = grad_z0;
xi(1) = xi0;
zeta(1) = zeta0;
dt = 0;
iterations = 300;

%% Method 2 using Euler

for ii = 1:iterations 
    deltar = 25; % step size along ray path
    xi(ii+1) = xi(ii) - (deltar)*((1/c(ii)^2)*grad_x(ii));  
    zeta(ii+1) = (zeta(ii) - (deltar)*((1/c(ii)^2)*grad_z(ii)));
    z(ii+1) = z(ii) + (deltar)*c(ii)*zeta(ii);
    x(ii+1) = x(ii) + (deltar)*c(ii)*xi(ii);
    c(ii+1) = ((c_s - c_e)/(h_s - h_e))*(z(ii+1) - h_e) + c_e;
    grad_z(ii+1) = grad_z(ii); 
    grad_x(ii+1) = grad_x(ii); 
    
    if z(ii+1)<=0 % point at which ray crosses horizontal ground plane
        deltar_temp = -z(ii)/(c(ii)*zeta(ii));
        xi(ii+1) = xi(ii) - (deltar_temp)*((1/c(ii)^2)*grad_x(ii));
        zeta(ii+1) = zeta(ii) - (deltar_temp)*((1/c(ii)^2)*grad_z(ii));
        z(ii+1) = z(ii) + (deltar_temp)*c(ii)*zeta(ii);
        x(ii+1) = x(ii) + (deltar_temp)*c(ii)*xi(ii);
        c(ii+1) = ((c_s - c_e)/(h_s - h_e))*(z(ii+1) - h_e) + c_e;
        grad_z(ii+1) = grad_z(ii);
        grad_x(ii+1) = grad_x(ii);
        zeta(ii+1) = -zeta(ii+1);
        deltar = deltar_temp;
    end
    s(ii+1) = s(ii) + deltar; % distance travelled by ray
    t(ii+1) = t(ii) + deltar*((1/c(ii) + 1/c(ii+1))/2); % travel time
end

%% Method 1
psi(1) = psi0; 
xx(1) = x0;
zz(1) = h0;
cc(1) = c0;
ss(1) = 0;
tt(1) = 0;
grad_z(1) = grad_z0;
incr = .05; % z-step size
if psi0 < 0
    deltaz = -incr;
else
    deltaz = incr;
end   

for ii = 1:11600
    Rc(ii) = (1/grad_z(ii)) * cc(ii)./cos(psi(ii)); % radius of curvature of segment
    psi(ii+1) = acos(cos(psi(ii)) + deltaz/Rc(ii));
    if deltaz < 0 
        psi(ii+1) = -psi(ii+1); % angle is negative when z-step is negative
    end
    dz(ii) = Rc(ii) * (cos(psi(ii+1)) - cos(psi(ii)));
    zz(ii+1) = zz(ii) + dz(ii);  
    if ~isreal(psi(ii+1))||zz(ii+1) < 0 % turning point of ray or point at which ray crosses horizontal ground plane 
        deltaz = -deltaz; % z-step changes direction
        if zz(ii+1) < 0 % case of ground crossing
            psi(ii+1) = acos(cos(psi(ii)) + deltaz/Rc(ii));% rays always travelling up after reflection
            psi(ii) = -psi(ii); % ray is now travelling in opposite direction
        else % case of turning point
            psi(ii+1) = -acos(cos(psi(ii)) + deltaz/Rc(ii));% rays always travelling down after turing point
        end
        dz(ii) = deltaz;
        zz(ii+1) = zz(ii) + dz(ii);
    end
    
    dx(ii) = -(Rc(ii) * (sin(psi(ii+1)) - sin(psi(ii))));
    xx(ii+1) = xx(ii) + dx(ii);
    cc(ii+1) = ((c_s - c_e)/(h_s - h_e))*(zz(ii+1) - h_e) + c_e;
    grad_z(ii+1) = (cc(ii+1) - cc(ii))/deltaz; 
    ss(ii+1) = ss(ii) + sqrt((xx(ii+1)-xx(ii))^2 + (zz(ii+1)-zz(ii))^2);
    arg1(ii+1) = atanh(tan(psi(ii+1)/2));
    arg2(ii+1) = atanh(tan(psi(ii)/2));
    tt(ii+1) = tt(ii) + abs((2/grad_z(ii+1))*(arg1(ii+1) - arg2(ii+1)));
end

% alternative method of calculating travel time (5.139)
tt2(1) = 0;
for ii = 1:length(xx)-1
    tt2(ii+1) = tt2(ii) + (ss(ii+1) - ss(ii))*((1/cc(ii+1) + 1/cc(ii))/2);
end

% alternative method of calculating travel time (5.139)
t2(1) = 0;
for ii = 1:length(x)-1
    t2(ii+1) = t2(ii) + (s(ii+1) - s(ii))*((1/c(ii+1) + 1/c(ii))/2);
end

figure
plot(x,z,'-b'); grid on

hold on
plot(xx,zz,'-r')

