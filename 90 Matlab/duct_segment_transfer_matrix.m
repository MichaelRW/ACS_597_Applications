function T = duct_segment_transfer_matrix(f, rho0, c, L, S)
k = 2*pi*f/c;
T = [cos(k*L), 1j*rho0*c/S*sin(k*L);...
1j*S/(rho0*c)*sin(k*L), cos(k*L)];