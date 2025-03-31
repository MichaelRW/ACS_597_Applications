function p = sum_of_monopoles(xyz_sources, Q_sources, xyz_receivers, f, rho0, c)
%sum_of_monopoles Compute sum of monopole sources at receiver locations.
%   p = sum_of_monopoles(xyz_sources, Q_sources, xyz_receivers, f, rho0, c)
%   computes the sum of monopole sources with strength Q_sources located at
%   xyz_sources, evaluated at receiver locations defined by xyz_receivers.
%   f is the frequency, rho0 is the fluid density, and c is the speed of
%   sound.
%
%   Inputs:
%       xyz_sources is a nSources-by-3 matrix of (x,y,z) source locations
%       Q_sources is a nSources-by-1 vector of complex source strengths
%       xyz_receivers is a nReceivers-by-3 matrix of (x,y,z) receiver
%       locations
%       f is a scaler frequency
%       rho0 is the fluid density
%       c is the speed of sound
%
%   Outputs:
%       p is a nReceivers-by-1 vector of complex pressure amplitudes at the
%       receiver locations

    arguments % This checks if the inputs are the same size as expected
        xyz_sources (:,3)
        Q_sources (:,1) {mustBeEqualFirstDim(xyz_sources, Q_sources)}
        xyz_receivers (:,3)
        f (1,1) % Only one frequency at a time
        rho0 (1,1)
        c (1,1)
    end

    % Define some constants
    w = 2*pi*f;
    k = w/c;
    nReceivers = size(xyz_receivers,1);
    nSources = size(xyz_sources,1);

    % Initialize pressure vector
    p = zeros(nReceivers,1,'like',1j);
    
    for iReceiver = 1:nReceivers
        xyz_R = xyz_receivers(iReceiver,:); % Receiver location
        for iSource = 1:nSources
            xyz_S = xyz_sources(iSource,:); % Source location
            r = sqrt(sum((xyz_R-xyz_S).^2)); % Distance from source to receiver
            Q = Q_sources(iSource);
            p_single_source = 1j*k*rho0*c/(4*pi*r)*Q*exp(-1i*k*r); % Pressure from single source
            p(iReceiver) = p(iReceiver) + p_single_source;
        end
    end
end

function mustBeEqualFirstDim(a,b)
        % Test for equal first dimension
    if ~isequal(size(a,1),size(b,1))
        eid = 'Size:notEqual';
        msg = 'First dimensions must match.';
        error(eid,msg)
    end
end