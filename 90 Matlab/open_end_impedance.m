function Z = open_end_impedance(f, rho0, c, L, S, flanged)
% Only for circular cross sections
a = sqrt(S/pi);
k = 2*pi*f/c;
if flanged
L_0 = 8*a/(3*pi);
else
L_0 = 0.61*a;
end
L_e = L+L_0;
Z = 1j*rho0*c/S*tan(k*L_e);