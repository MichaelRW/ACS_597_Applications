% % Displacement Versus Frequency Ratio
% 
% Fo = 1;
% wf = 4;  % Force frequency, radians\s
% 
% lambda = 1;
%     k = (2*pi)/lambda;
% 
% phase_offset = 0;
% 
% w = 5;
%     r = wf / w;
% 
% time_indices = 0:1e-2:5;
% 
% h_x = @( Fo, k, wf, time_indices, phase_offset, r, epsilon ) ( Fo./k.*sin(4*time_indices - phase_offset) ) ./ ( sqrt( (1 - r^2)^2 + (2*r*epsilon)^2 ) );
% 
% figure( ); ...
%     epsilon = 0;
%         % plot( time_indices, h_x( Fo, k, wf, time_indices, phase_offset, r, epsilon ) );  hold on;
%     %
%     for wf = 10:-1:1
%         wf
%         epsilon = 0.25;
%         r = wf / w;
%             plot( h_x( Fo, k, wf, time_indices, phase_offset, r, epsilon ) );  hold on;
%         keyboard
%     end
%     %
%     legend( '', '' );
%     xlabel( 'Frequency Ratio [WU]' );  ylabel( 'Normalized Maximum Displacment [WU]' );