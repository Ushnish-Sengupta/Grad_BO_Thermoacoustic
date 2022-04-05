function [gam,zet,n,Ru,Rd,rhu,rhd,kus,kds,cu,cd] = unwrap_param(param)
% unwrap param
%
% INPUTS
% param     structure containing nondimensional parameters
%
% OUTPUTS
% gam
% zet
% n
% Ru
% Rd
% rhu
% rhd
% kus
% kds
% cu
% cd
%
% The other parameters (tau, X_w, L_w, X_h, L_h) are embedded within the
% matrices held within the mat structure.

% Gas parameters
gam = param.gam;
zet = param.zet;

% Heater parameters
n   = param.n;

% Boundary conditions
Ru  = param.Ru;
Rd  = param.Rd;
rhu = param.rhu;
rhd = param.rhd;

% Robin boundary condition coefficients, calculated from Ru and Rd
kus = param.kus;
kds = param.kds;

% Arbitrary boundary condition coefficients
cu  = param.cu;
cd  = param.cd;

end
