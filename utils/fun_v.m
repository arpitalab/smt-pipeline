function [output] = fun_v(x)
% output = fun_v(x)
%  Calculates the expected value of the MSD for generalized diffusion with
%  proper consideration of static and dynamic localization errors
% From Backlund et al. Phys. Rev E
% Equation 8 of the paper
% The expression in Backlund et al. assumes that the exposure time is the
% same as the sampling interval. 

delta_t = 0.2; % sampling time
exposure = 0.01; % exposure time
k = 1:25;
fk = (k+1).^(2+x(2))+(k-1).^(2+x(2))-2*k.^(2+x(2)); % lag time
output = (2*x(1)*delta_t.^x(2).*(fk)-2*x(1)*exposure.^(x(2)))./((x(2)+1).*(x(2)+2))+x(3);
output = (output');

end

