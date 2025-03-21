
function [ FRF ] = nDOF_direct_solution( masses, stiffnesses, dampings, freq, FRF_type )

% nDOF_direct_solution

% Calculate frequency response functions given system coefficients.


% FRF = nDOF_direct_solution( masses, stiffnesses, damping, freq, FRF_type )
 
% Inputs:
%
%       masses:  An nDOF vector of masses.
%       stiffnesses:  An nDOF+1 vector of stiffnesses, with the first and last entries being
%           connections to ground.
%       dampings:  An nDOF+1 vector of stiffnesses, with the first and last entries being
%       connections to ground.
%       freq:  An nFreq vector of frequencies.
%       FRF_type:  Is one of:  'admittance', 'mobility', 'accelerance', 'stiffness', 'impedance', or 'mass'.

% Output(s):
%
%       FRF:  An nFreq-by-nDOF-by-nDOF matrix containing the requested frequency response 
%           function between all pairs of degrees of freedom.


nDOF = length( masses );
    assert( length( stiffnesses ) == nDOF + 1, '*** Stiffnesses must be 1 longer than masses. ***' );
%
M = diag(masses);


Kdiag = diag( stiffnesses(1:end-1) ) + diag( stiffnesses(2:end) );
    Kupper = diag( stiffnesses(2:end-1), 1 );
    Klower = diag( stiffnesses(2:end-1), -1 );
        K = Kdiag + Kupper + Klower;


Cdiag = diag( dampings(1:end-1)) + diag( dampings(2:end) );
    Cupper = diag( dampings(2:end-1), 1 );
    Clower = diag( dampings(2:end-1), -1 );
        C = Cdiag + Cupper + Clower;


FRF_type = lower( FRF_type );


nFreq = length( freq );


FRF = zeros( nFreq, nDOF, nDOF, 'like', 1j );

for iFreq = 1:1:nFreq

    w = 2*pi*freq( iFreq );
        A = -w^2*M+1j*w*C+K;  % Dynamic stiffness; integrate and invert to get other forms.
    

    switch FRF_type

        case 'stiffness'
            FRF(iFreq,:,:) = A;

        case 'impedance'
            FRF(iFreq,:,:) = A/(1j*w);

        case 'mass'
            FRF(iFreq,:,:) = -A/w^2;

        case {'admittance','receptance','compliance'}
            FRF(iFreq,:,:) = A\eye(nDOF);

        case 'mobility'
            FRF(iFreq,:,:) = 1j*w*A\eye(nDOF);

        case {'accelerance','intertance'}
            FRF(iFreq,:,:) = -w^2*A\eye(nDOF);

    end

end


