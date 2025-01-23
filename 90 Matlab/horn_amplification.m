

function A = horn_amplification(f, d_i, L_i, rho0, c)

% d_i is vector segment diameters (outlet to inlet)
% L_i is vector of segment lengths (outlet to inlet)
% A is amplification (inlet to outlet)


S_i = pi*d_i.^2/4;
flanged = false;
nFreq = length(f);
nSegments = length(L_i);
A = zeros(nFreq,1);
for iFreq = 1:nFreq
T_total = [1 0; 0 1]; % Start with identity matrix
for iSegment = 1:nSegments
L = L_i(iSegment);
S = S_i(iSegment);
T_segment = duct_segment_transfer_matrix(f(iFreq), rho0, c, L, S);
T_total = T_segment*T_total; % New segment goes on the left
end
T11 = T_total(1,1);
T12 = T_total(1,2);
Z = open_end_impedance(f(iFreq), rho0, c, 0, S(1), flanged);
A(iFreq) = -10*log10(abs(T11+T12/Z)^2);
end


