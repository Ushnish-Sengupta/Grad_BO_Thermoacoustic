function [wr,pwrpx,pwrpX_h,pwrpL_h,pwrpX_w,pwrpL_w] = fun_wr(param,x)
% fun_wr
%
% Set the nondimensional measurement envelope, w, and its partial 
% derivatives with respect to the parameters and the coordinate
% 
% Then create (w/rho) and output this and its partial derivatives
%
% INPUTS: 
% param      structure containing the parameters
% x          column vector containing the gridpoints
%
% OUTPUTS:
% wr(x)      w/rho (x) 
% pwrpx(x)   partial derivative of w/rho(x) w.r.t. x
% pwrpX_h(x) partial derivative of w/rho(x) w.r.t. X_h
% pwrpL_h(x) partial derivative of w/rho(x) w.r.t. L_h
% pwrpX_w(x) partial derivative of w/rho(x) w.r.t. X_w
% pwrpL_w(x) partial derivative of w/rho(x) w.r.t. L_w


% Define the nondimensional measurement envelope, which integrates to 1 over the domain
% Then divide it by rho and output it and its partial derivatives wrt
% parameters. 

% Generate the measurement envelope, w(x), which integrates to 1, and its derivatives 
w = exp(-(x-param.X_w).^2/param.L_w^2)/sqrt(pi)/param.L_w;
% partial derivative wrt x
pwpx  = w .* (-2*(x-param.X_w)/param.L_w^2);
% partial derivative wrt x_m
pwpX_w = - pwpx;
% partial derivative wrt l_m
pwpL_w = w .* (2*(x-param.X_w).^2/param.L_w^3) - w / param.L_w;

% Read in the density, rh(x) and its derivatives w.r.t x, x_f, l_f
[rh,prhpx,prhpX_h,prhpL_h] = fun_rh(param,x);
% Generate the measurement envelope divided by density, ww/rh
wr = w./rh;
% partial derivative wrt x
pwrpx   = (pwpx   .* rh - prhpx   .* w)./(rh.^2);
% partial derivative wrt X_h
pwrpX_h = - prhpX_h .* w ./(rh.^2);
% partial derivative wrt L_h
pwrpL_h = - prhpL_h .* w ./(rh.^2);
% partial derivative wrt x_m
pwrpX_w = pwpX_w ./ rh;
% partial derivative wrt l_m
pwrpL_w = pwpL_w ./ rh;

end

