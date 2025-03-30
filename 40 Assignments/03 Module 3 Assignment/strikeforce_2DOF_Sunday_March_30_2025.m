%% 2 DOF Example

% masses = [ 50  1 ];
% stiffnesses = [ 100e3  1800  0 ];
% 
% dampings = [ 0  0  0 ];
% % dampings = [ 0  (1 + 0.1*1j)  (1 + 0.1*1j) ];
% 
% frequencies = 0:0.01:40;
% 
% FRF = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );  % 4001-by-2-by-2
% 
% 
% admittance = zeros( numel( frequencies ), 1 );
% 
% for index = 1:1:numel( frequencies )
%     temp = diag( squeeze( FRF( index, :, : ) ) );
%         temp2 = temp(1) + temp(2)*1j;
%         admittance( index ) = abs( temp2 );
% end
% 
% clear temp1 temp2;
% 
% 
% % figure( ); ...
% %     semilogy( frequencies, admittance );  grid on;
% %     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
% 
% 
% figure( 'Name', 'Admittance' ); ...
%     semilogy( frequencies, abs( FRF( :, 1, 1 ) ) );  hold on;
%     semilogy( frequencies, abs( FRF( :, 1, 2 ) ) );
%     semilogy( frequencies, abs( FRF( :, 2, 1 ) ) );
%     semilogy( frequencies, abs( FRF( :, 2, 2 ) ) );  grid on;
%         legend( 'M1 Movement with F1', 'M1 Movement with F2', 'M2 Movement with F1', 'M2 Movement with F2', 'Location', 'North' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
%     axis( [ 5.5 10.5  1e-7 1e2 ] );


% figure( 'Name', 'Dynamic Stiffness (Reciprocal of Admittance)' ); ...
%     semilogy( frequencies, 1 ./ abs( FRF( :, 1, 1 ) ) );  hold on;
%     semilogy( frequencies, 1 ./ abs( FRF( :, 1, 2 ) ) );
%     semilogy( frequencies, 1 ./ abs( FRF( :, 2, 1 ) ) );
%     semilogy( frequencies, 1 ./ abs( FRF( :, 2, 2 ) ) );  grid on;
%         legend( 'M1 Movement with F1', 'M1 Movement with F2', 'M2 Movement with F1', 'M2 Movement with F2', 'Location', 'North' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );