function [param,ref] = fun_nondim(param_dim)
% fun_nondim
%
% Convert the dimensional parameters to nondimensional parameters
%
% INPUTS: 
% param_dim structure containing the dimensional parameters
%
% OUTPUTS:
% param     structure containing the non-dimensional parameters
% ref       characteristic dimensions for length, pressure, and velocity

%% Calculate the reference scales
ref.l_dim = param_dim.X;
ref.p_dim = param_dim.pbar;
ref.u_dim = sqrt(param_dim.gam * param_dim.pbar / param_dim.rhobar);

%% Calculate the nondimensional parameters
% Gas parameters
param.gam = param_dim.gam;
param.zet = (param_dim.gam-1)/param_dim.gam;

% Heater parameters
param.n   = param_dim.n; % (This is dimensionless in 1D but not 2D or 3D)
param.tau = param_dim.tau / ref.l_dim * ref.u_dim;

% Heater position parameters
param.X_w = param_dim.X_w / ref.l_dim;
param.L_w = param_dim.L_w / ref.l_dim;
param.X_h = param_dim.X_h / ref.l_dim;
param.L_h = param_dim.L_h / ref.l_dim;

% Boundary conditions
param.Ru  = param_dim.Ru;
param.Rd  = param_dim.Rd;
param.rhu = param_dim.rhu * ref.u_dim^2 / ref.p_dim;
param.rhd = param_dim.rhd * ref.u_dim^2 / ref.p_dim;

% Robin boundary condition coefficients calculated from Ru and Rd
[kus,kds,~,~] = fun_kukd(param);
param.kus = kus;
param.kds = kds;

% Arbitrary boundary condition coefficients
param.cu  = param_dim.cu;
param.cd  = param_dim.cd;

% Starting point
param.mode = param_dim.mode;

end
