function [h,phpX_h,phpL_h] = fun_h(param,x)
% fun_h
%
% Set the nondimensional heat release envelope, h, and its partial 
% derivatives with respect to the parameters. 
%
% INPUTS: 
% param     structure containing the parameters
% x         column vector containing the gridpoints
%
% OUTPUTS:
% h(x)      heat release integral (integrates to 1 over the domain [0,1])
% phpX_h(x) partial derivative of h(x) w.r.t. X_h
% phpL_h(x) partial derivative of h(x) w.r.t. L_h

% Heat release envelope
h       = exp(-(x-param.X_h).^2/param.L_h^2)/sqrt(pi)/param.L_h;
% partial derivative w.r.t. X_h
phpX_h  = (2*(x-param.X_h)/param.L_h^2) .* h;
% partial derivative w.r.t. L_h
phpL_h  = (h .* (2*(x-param.X_h).^2) / param.L_h^3 - h/param.L_h);

end
