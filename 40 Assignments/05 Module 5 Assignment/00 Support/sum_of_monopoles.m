
function [ p ] = sum_of_monopoles( xyz_sources, Q_sources,  xyz_receivers, f, rho0, c )

% Compute sum of monopole sources at receiver locations.


% Inputs:

%   xyz_sources:  Is an nSources-by-3 matrix of (x,y,z) source locations.
%   Q_sources:  Is an nSources-by-1 vector of complex-valued source strengths.
%   xyz_receivers:  Is am nReceivers-by-3 matrix of (x,y,z) receiver locations.
%   f:  Is a scaler frequency.
%   rho0:  Is the fluid density.
%   c:  Is the speed of sound.

%   Outputs:

%   p:  Is an nReceivers-by-1 vector of complex pressure amplitudes at the receiver locations.


arguments
    xyz_sources ( :, 3 )
    Q_sources ( :, 1 )  { mustBeEqualFirstDim( xyz_sources, Q_sources ) }
    xyz_receivers ( :, 3 )
    f ( 1, 1 )  % Only one frequency at a time.
    rho0 ( 1, 1 )
    c ( 1, 1 )
end


w = 2*pi*f;  % radians/s
    k = w / c;  % Wave number ( 1/m ) 


nSources = size( xyz_sources, 1 );
nReceivers = size( xyz_receivers, 1 );


p = zeros( nReceivers, 1, 'like', 1j );

    
for iReceiver = 1:1:nReceivers

        xyz_R = xyz_receivers( iReceiver, : );  % Receiver location.

        for iSource = 1:1:nSources

            xyz_S = xyz_sources( iSource, : );  % Source location.
                r = sqrt( sum( ( xyz_R - xyz_S ).^2 ) );  % Distance from source to receiver.

            Q = Q_sources( iSource );

            p_single_source = 1j*k*rho0*c/(4*pi*r)*Q*exp(-1i*k*r);  % Pressure from a single source.

            p( iReceiver ) = p( iReceiver ) + p_single_source;

        end

end


end  % End:  function [ p ] = sum_of_monopoles( xyz_sources, Q_sources,  xyz_receivers, f, rho0, c )





function mustBeEqualFirstDim( a, b )  % Test for equal first dimension.

    if ( ~isequal( size(a, 1), size(b, 1) ) )        
            error( 'Size:notEqual', '*** First dimensions must match. ***' );
    end

end


