function [ds_int] = fun_ds_CA(mat,param,T,N,pl,s,pr,discretization)
% fun_ds_CA
%
% Evaluate the internal sensitivities in the continuous adjoint framework
%
% INPUTS
% mat            structure containing the FD or FE matrices (see mat_**.m)
% param          structure containing the parameters (see fun_nondim.m) 
% T              structure containing the transformation matrices (see mat_AC_**_**.m) 
% N              N+1 is the number of gridpoints in FD. N is the number of elements in FE
% pl             the adjoint pressure eigenvector
% s              the eigenvalue
% pr             the direct pressure eigenvector
% discretization the discretization (FD, FE, or WD)
%
% OUTPUTS
% ds_int         structure containing the internal sensitivities

%% Unwrap param
[gam,zet,n,~,~,rhu,rhd,kus,kds,~,~] = unwrap_param(param);

switch discretization
    case 'FEW'
        % Unwrap mat_FE
        [D01,M11,M00,~,~,t,h,wr,~] = unwrap_mat_FE(mat);
    case {'FDS','FDW'}
        % Unwrap mat_FD and, in order to use the code below, rename 
        % D as D01, M as M11, M as M00.
        [D01,M11,~,t,h,wr,~] = unwrap_mat_FD(mat); M00 = M11;
        % Relabel T in order to be able to use code below
        T = fun_relabel_T(T);
    case 'SBP'
        % Unwrap mat_SBP and, in order to use the code below, rename
        % D1 as D01, M  as M11, M  as M00.
        [D01,~  ,M11,~,~,~,~,~,~,~,~,t,h,wr,~] = unwrap_mat_SBP(mat); M00 = M11;
        % Relabel T in order to be able to use code below
        T = fun_relabel_T(T);
end

%% Calculate the inner product, which is the same for all sensitivities
% Calculate dAds in continuous form without boundary conditions
dAds = - zet * n * diag(-t.*exp(-s*t)) * h * wr.' * M00 * D01;
% Calculate dGds in continuous form without boundary conditions
dGds = dAds - 2*s*eye(N+1)/gam; % n.b. eye(N+1), not M11, required for CA. 
% Calculate dGds, including the Robin boundary conditions, which have k(s)
% switch discretization
ip = -(pl'      * M11 * dGds * pr) ...      % Bulk term
     -(pl(N+1)' * kus/rhu    * pr(N+1)) ... % upstream boundary term
     -(pl(  1)' * kds/rhd    * pr(1  ));    % downstream boundary term
 
% Base state sensitivities in continuous form
ds_int.n   = - pl' * M11 * zet         * diag(exp(-s*t)) * h    * (wr.' * M00 *  D01 * pr)   / ip;
ds_int.t   =   pl' .*     (zet * n * s * diag(exp(-s*t)) * h).' * (wr.' * M00 *  D01 * pr)   / ip;
ds_int.h   = - pl'       * zet * n     * diag(exp(-s*t))        * (wr.' * M00 *  D01 * pr)   / ip;
ds_int.wr  = - pl' * M11 * zet * n     * diag(exp(-s*t)) * h                  * (D01 * pr).' / ip;
ds_int.v   = - pl' * D01'                                                    .* (D01 * pr).' / ip;
ds_int.ku  =   pl(N+1)'  *  s/rhu  *  pr(N+1) / ip; 
ds_int.kd  =   pl(1)'    *  s/rhd  *  pr(1)   / ip; 
% Feedback sensitivities in continuous form
ds_int.mru =  (pl' * T.MR10) .* (T.U01 * pr).' / ip; % row .* row sensitivity to feedback from u into mass eq.
ds_int.mrp =  (pl' * T.MR11) .* (T.P11 * pr).' / ip; % row .* row sensitivity to feedback from p into mass eq.
ds_int.fru =  (pl' * T.FR10) .* (T.U01 * pr).' / ip; % row .* row sensitivity to feedback from u into momentum eq.
ds_int.frp =  (pl' * T.FR10) .* (T.P01 * pr).' / ip; % row .* row sensitivity to feedback from p into momentum eq.
ds_int.qpu =  (pl' * T.QP10) .* (T.U01 * pr).' / ip; % row .* row sensitivity to feedback from u into energy eq.
ds_int.qpp =  (pl' * T.QP11) .* (T.P11 * pr).' / ip; % row .* row sensitivity to feedback from p into energy eq.

end

function T = fun_relabel_T(T)
% fun_relabel_T
%
% relabel structure variables so that routines can be shared between FD, FE, SBP, and WD
%
% INPUTS
% T                     structure containing variables to be relabeled
% target_discretization target discretization
%
% OUTPUTS
% T                     structure containing relabeled variables

% Relabel FD, SBP & WD structure as FE structure
T.P11 = T.P;
T.P01 = T.P;
T.U01 = T.U;
T.MR10 = T.MR;
T.MR11 = T.MR;
T.FR10 = T.FR;
T.QP10 = T.QP;
T.QP11 = T.QP;
T = rmfield(T,{'U','P','MR','FR','QP'});

end

