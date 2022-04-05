function [ds_ext] = fun_ext_int_CA(ds_int,param,x0,x1,N,M00,M11)
% fun_ext_int_CA
%
% Evaluate the external sensitivities in the continuous framework
%
% INPUTS
% ds_int    structure containing the internal sensitivities (see fun_ds_CA.m)
% param     structure containing the parameters
% x0        positions of the gridpoints x0_i (x0 = x for FD and WD)
% x1        positions of the gridpoints x1_i (x1 = x for FD and WD)
% N         N+1 is the number of gridpoints in FD. N is the number of elements in FE
% M00       mass matrix: <f0,g0> = f0^H * M00 * g0 (M00 = M for FD and WD)
% M11       mass matrix: <f1,g1> = f1^H * M11 * g1 (M00 = M for FD and WD)
%
% OUTPUTS
% ds_ext    structure containing the external sensitivities

% NOTE
% This is equivalent to fun_ext_int_DA.m with M00 = M11 = I 

%% Read in deriviatives of internal parameters w.r.t external parameters
[~  ,     dhdX_h,  dhdL_h                      ] = fun_h(param,x1);
if length(ds_int.wr) == N
    [~  , ~, dwrdX_h, dwrdL_h, dwrdX_w, dwrdL_w] = fun_wr(param,x0);
    Mwr = M00;
else
    [~  , ~, dwrdX_h, dwrdL_h, dwrdX_w, dwrdL_w] = fun_wr(param,x1);
    Mwr = M11;
end
if length(ds_int.v) == N
    [rh, ~, drhdX_h, drhdL_h                   ] = fun_rh(param,x0);
    Mv = M00;
else
    [rh, ~, drhdX_h, drhdL_h                   ] = fun_rh(param,x1);
    Mv = M11;
end
[~  , ~                    , dkudRu  , dkddRd  ] = fun_kukd(param);
% Create dv/d* from d(rho)/d*
dvdX_h = (-1./rh.^2 .* drhdX_h); 
dvdL_h = (-1./rh.^2 .* drhdL_h); 

%% Initialize ds_ext to contain the same fields as param
% (we need to include all fields, even if zero, for the Taylor Test)
ds_ext = param;
names = fieldnames(ds_ext);
for nn = 1:length(names)
    f = getfield(ds_ext,names{nn});
    ds_ext = setfield(ds_ext,names{nn},0*f');
end

%% Calculate ds/d* for the external parameters
ds_ext.n   = ds_int.n;
ds_ext.X_w = ds_int.wr * Mwr * dwrdX_w;
ds_ext.L_w = ds_int.wr * Mwr * dwrdL_w;
ds_ext.X_h = ds_int.h  * M11 * dhdX_h ...
           + ds_int.v  * Mv  * dvdX_h ...
           + ds_int.wr * Mwr * dwrdX_h;
ds_ext.L_h = ds_int.h  * M11 * dhdL_h ...
           + ds_int.v  * Mv  * dvdL_h ...
           + ds_int.wr * Mwr * dwrdL_h;
ds_ext.tau = ds_int.t  * M11 * ones(N+1,1);

% If Ru is -1 or +1 then Dirichlet or Neumann conditions will be applied and ds/d(Ru) is not relevant
if param.Ru == +1 || param.Ru == -1
    ds_ext.Ru = 0;
else
    ds_ext.Ru  = ds_int.ku * dkudRu;
end

% If Rd is -1 or +1 then Dirichlet or Neumann conditions will be applied and ds/d(Rd) is not relevant
if param.Rd == +1 || param.Rd == -1
    ds_ext.Rd = 0;
else
    ds_ext.Rd  = ds_int.kd * dkddRd;
end

end
