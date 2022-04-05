function [rh,prhpx,prhpX_h,prhpL_h] = fun_rh(param,x)
% fun_rh
%
% Set the nondimensional density profile and its partial derivatives with 
% respect to the co-ordinate and the parameters.
%
% INPUTS: 
% param      structure containing the parameters
% x          column vector containing the gridpoints
%
% OUTPUTS:
% rh(x)      density profile
% prhpx      partial derivative of rh(x) w.r.t. x
% prhpX_h(x) partial derivative of rh(x) w.r.t. X_h
% prhpL_h(x) partial derivative of rh(x) w.r.t. L_h

% tanh profile (at the flame) between upstream and downstream densities 
rh      =   param.rhu + (param.rhd-param.rhu)*0.5*(1 + tanh((x-param.X_h)/param.L_h));
prhpx   =               (param.rhd-param.rhu)*0.5*(sech((x-param.X_h)/param.L_h)).^2 / param.L_h;
prhpX_h =             - (param.rhd-param.rhu)*0.5*(sech((x-param.X_h)/param.L_h)).^2 / param.L_h;
prhpL_h =               (param.rhd-param.rhu)*0.5*(sech((x-param.X_h)/param.L_h)).^2 .* -(x-param.X_h)/param.L_h^2;

end
