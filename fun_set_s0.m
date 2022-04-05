function [s0] = fun_set_s0(param)
% fun_set_s0
%
% Use a travelling wave method to estimate s0 for the requested mode
%
% INPUTS
% param         structure containing the parameters
% param.mode    The requested chamber mode 
%
% OUTPUTS
% s0            Estimate of the eigenvalue, s, for the requested mode

% Assume density is rhu upstream of x_f and rhd downstream of x_f and
% ignore heat release. Calculate time taken for a wave to travel from 
% x_f to x=0 and back, noting that the nondim speed of sound upstream 
% of x_f is 1
tau_u = 2*param.X_h/1;
% Calculate time taken for a wave to travel from x_f to x=1 and back,
% noting that the nondim speed of sound downstream of x_f is sqrt(rhu/rhd)
tau_d = 2*(1-param.X_h)/sqrt(param.rhu/param.rhd);
% Find the expected s of this mode
s0 = (2*param.mode*pi*1i + log(param.Ru * param.Rd))/(tau_u + tau_d);
% Display to screen
disp(['Initial guess for     s= ',num2str(s0,'%+.16f')])

end
