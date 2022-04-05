function [] = Sup_003
% Sup_003
%
% Check that the influence of the first hot mesh is the same as the second
% when they have the same power and are at the same location. 
%
% This outputs s (for the first hot wire) and ds (for the second). They
% should differ slightly because the heat release from the first hot wire
% is applied over the profile h(x), while the second is applied at a point.

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');
% Set the reflection coefficients to -1 in order to remove any damping
param_dim.Ru = -1;
param_dim.Rd = -1;
% Set the measurement position to be the same as the flame position
param_dim.X_w = param_dim.X_h;
% Set teh measurement length to be the same as the flame length
param_dim.l_m = param_dim.L_h;

%% Calculate the reference scales and the nondimensional parameters
[param,ref] = fun_nondim(param_dim);

%% Set the numerical scheme, the max number of iterations, and the starting value of s
scheme.N     = 101;
scheme.s0    = fun_set_s0(param);

% Calculate the eigenmode and the internal sensitivities
[emode,ds_int,~] = fun_Helm('FEW','DA','nonlin',param,scheme);
% Calculate the dimensional s (in radians/sec)
s_dim = emode.s * ref.u_dim / ref.l_dim;

% Set the power of the secondary heater to be the same as the first
Q = param_dim.Q;
% Calculate nw by using the same FTF and U as the main heater
nw_dim = Q * param_dim.FTF/param_dim.U;
% Set the tube area 
A_dim = param_dim.S_c; % tube x-sect area
tau_dim = param_dim.tau;

% Extract the non-dimensional feedback sensitivity from velocity into heat
% input (in continuous nondim form)
dsdqpu_nondim = ds_int.qpu/emode.M00;

% Calculate ds
ds = (ref.u_dim/ref.l_dim) * nw_dim / (ref.p_dim * A_dim) * exp(-s_dim*tau_dim) * dsdqpu_nondim;

% Find the value of ds at the primary heater location
ind = dsearchn(emode.x0,param.X_h);

% extract the value of ds there
ds_secondary = ds(ind);

% Compare the two values
disp([' s = ',num2str(real(s_dim),4),' rad/sec'])
disp(['ds = ',num2str(real(ds_secondary),4),' rad/sec'])
disp('If the feedback sensitivitiy is correct, s should be close to ds.')

end


