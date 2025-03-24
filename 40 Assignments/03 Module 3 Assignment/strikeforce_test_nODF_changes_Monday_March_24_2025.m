%% 1 DOF Example

% masses = 50;
% stiffnesses = [ 100e3  0 ];
% dampings = [ 0  0 ];
%     % dampings = [ 0  (1 + 0.1*1j) ];
% 
% wo = sqrt( 100e3 / 50 );
%     fo = wo/(2*pi);  % 7.1 Hz
% 
% frequencies = 0:0.1:100;
% 
% 
% admittance = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'admittance' );  % 4001-by-2-by-2
% %
% figure( ); ...
%     loglog( frequencies, abs( admittance ) );  grid on;
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Admittance [$\frac{m}{N}$]' );
% 
% 
% impedance_original_function = nDOF_direct_solution_org( masses, stiffnesses, dampings, frequencies, 'impedance' );  % 4001-by-2-by-2
% impedance_updated_function = nDOF_direct_solution( masses, stiffnesses, dampings, frequencies, 'impedance' );  % 4001-by-2-by-2
% %
% figure( ); ...
%     loglog( frequencies, abs( impedance_original_function ) );  hold on;
%     loglog( frequencies, abs( impedance_updated_function ) );  grid on;
%         legend( 'Original Function', 'Updated Function', 'Location', 'North' );
%     xlabel( 'Frequency [Hz]' );  ylabel( 'Impedance [$\frac{N\cdots}{m}$]' );