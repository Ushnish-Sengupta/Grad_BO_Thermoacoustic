function [A,C,dAds,T,dA,ds_int] = mat_AC_FDW_DA(mat,param,N,s)
% mat_AC_FDW_DA
%
% Generate A and C matrices for Discrete Adjoint Finite Difference in the weak form 
%
% INPUTS: 
% mat.D     first derivative matrix: dp/dx = D*p
% mat.M     mass matrix: <f,g> = f^H * M * g
% mat.V     matrix containing the specific volume v(x_i) along the diagonal
% mat.t     column vector containing time delay t(x_i)
% mat.h     column vector containing heat release envelope h(x_i) 
% mat.wr    column vector containing measurement envelope / density wr(x_i)
% mat.dwrdx column vector containing d(wr)/d(x) at x_i
% param     structure containing the parameters
% N         N+1 is the number of gridpoints
% s         value of s at which to calculate the matrices
%
% OUTPUTS:
% A         matrix representing acoustics and heat release, with boundary conditions applied 
% C         matrix representing 1/gamma, with boundary conditions applied
% dA/ds     dA/ds
% T.P       matrix that maps p to p
% T.U       matrix that maps p to u
% T.MR      matrix that maps adjp to the Receptivity of the mass equation to injection of mdot/rho
% T.FR      matrix that maps adjp to the Receptivity of the momentum equation to injectrion of f/rho
% T.QP      matrix that maps adjp to the Receptivity of the energy equation to injection of q/p
% dA.*      structure containing all the sensitivities: dA/d(param)
% ds_int.*  initialized structure for ds of all the internal parameters

%% Unwrap mat 
[D,M,V,t,h,wr,~] = unwrap_mat_FD(mat);

%% Unwrap param
[gam,zet,n,Ru,Rd,rhu,rhd,kus,kds,cu,cd] = unwrap_param(param);

%% Generate A, C, and dAds matrices
A    = - D'*M*V*D - M * zet * n * diag(exp(-s*t))     * h * wr.' * M * D;
C    =   M/gam;
dAds =            - M * zet * n * diag(-t.*exp(-s*t)) * h * wr.' * M * D;

%% Apply the boundary conditions to A, C, and dAds
[A,C,dAds] = fun_bcs_weak(A,C,dAds,N,Ru,cu,kus,rhu,Rd,cd,kds,rhd,s);

%% Define the conversion / transformation matrices
if nargout >= 4
    % Conversion from pr to pr
    T.P = eye(N+1);
    % Conversion from pr to ur
    T.U = - V * D / s;
    % Conversion from pl to Receptivity of mass equation to injection of mdot/rho
    T.MR = s * M;
    % Conversion from pl to Receptivity of momentum equation to injectrion of f/rho
    T.FR = D' * M;
    % Conversion from pl to Receptivity of energy equation to injection of q/p
    T.QP = s * zet * M;
end

%% Generate all the sensitivities of matrix A: dA.* = dA.*l  x  d.*  x  dA.*r
if nargout >= 5
    % initialize the ds structure
    ds_int.n    = 0;
    ds_int.t    = zeros(1,N+1);
    ds_int.h    = zeros(1,N+1);
    ds_int.wr   = zeros(1,N+1);
    ds_int.v    = zeros(1,N+1);
    ds_int.ku   = 0;
    ds_int.kd   = 0;
    ds_int.mru = zeros(1,N+1);
    ds_int.mrp = zeros(1,N+1);
    ds_int.fru = zeros(1,N+1);
    ds_int.frp = zeros(1,N+1);
    ds_int.qpu = zeros(1,N+1);
    ds_int.qpp = zeros(1,N+1);

    % Generate sensitivities in the bulk (before boundary conditions)
    % Sensitivity to n
    dA.n   = - M * zet         * diag(exp(-s*t)) * h * wr.' * M * D;
    % Sensitivity to t
    dA.tl  = - M;
    dA.tr  =     - zet * n * s * diag(exp(-s*t)) * h * wr.' * M * D;
    % Sensitivity to h
    dA.hl  = - M * zet * n     * diag(exp(-s*t));
    dA.hr  =                                           wr.' * M * D;
    % Sensitivity to wr
    dA.wrl = - M * zet * n     * diag(exp(-s*t)) * h;
    dA.wrr =                                                  M * D;
    % Sensitivity to v (1/rho)
    dA.vl  = - D'*M;
    dA.vr  =         D;

    % Sensitivity to feedback from u into mass/rho equation
    dA.mrul = T.MR; 
    dA.mrur = T.U; 
    % Sensitivity to feedback from p into mass/rho equation
    dA.mrpl = T.MR; 
    dA.mrpr = T.P;
    % Sensitivity to feedback from u into momentum/rho equation
    dA.frul = T.FR; 
    dA.frur = T.U; 
    % Sensitivity to feedback from p into momentum/rho equation
    dA.frpl = T.FR;
    dA.frpr = T.P;
    % Sensitivity to feedback from u into energy/p equation 
    dA.qpul = T.QP;
    dA.qpur = T.U; 
    % Sensitivity to feedback from p into energy/p equation 
    dA.qppl = T.QP;
    dA.qppr = T.P;

    % Apply downstream boundary conditions
    if Rd == -1
        % homogenous Dirichlet on full matrices
        dA.n(1,:)      = zeros(1,N+1); % top row of dA.n    -> 0
        dA.n(:,1)      = zeros(N+1,1); % lef col of dA.n    -> 0
        % homogenous dirichlet on split matrices
        dA.tl(1,:)     = zeros(1,N+1); % top row of dA.tl   -> 0
        dA.tr(:,1)     = zeros(N+1,1); % lef col of dA.tr   -> 0
        dA.hl(1,:)     = zeros(1,N+1); % top row of dA.hl   -> 0
        dA.hr(1)       = 0;            % lef col of dA.hr   -> 0
        dA.wrl(1)      = 0;            % top row of dA.wrl  -> 0
        dA.wrr(:,1)    = zeros(N+1,1); % lef col of dA.wrr  -> 0
        dA.vl(1,:)     = zeros(1,N+1); % top row of dA.vl   -> 0
        dA.vr(:,1)     = zeros(N+1,1); % lef col of dA.vr   -> 0
        dA.kd = 0;
        dA.mrul(1,:)   = zeros(1,N+1); % top row of dA.mrul -> 0
        dA.mrur(:,1)   = zeros(N+1,1); % lef col of dA.mrur -> 0
        dA.mrpl(1,:)   = zeros(1,N+1); % top row of dA.mrpl -> 0
        dA.mrpr(:,1)   = zeros(N+1,1); % lef col of dA.mrpr -> 0
        dA.frul(1,:)   = zeros(1,N+1); % top row of dA.frul -> 0
        dA.frur(:,1)   = zeros(N+1,1); % lef col of dA.frur -> 0
        dA.frpl(1,:)   = zeros(1,N+1); % top row of dA.frpl -> 0
        dA.frpr(:,1)   = zeros(N+1,1); % lef col of dA.frpr -> 0
        dA.qpul(1,:)   = zeros(1,N+1); % top row of dA.qpul -> 0
        dA.qpur(:,1)   = zeros(N+1,1); % lef col of dA.qpur -> 0
        dA.qppl(1,:)   = zeros(1,N+1); % top row of dA.qppl -> 0
        dA.qppr(:,1)   = zeros(N+1,1); % lef col of dA.qppr -> 0
    elseif Rd == +1
        % homogenous neumann (no action required)
        dA.kd = 0;
    else
        % Robin (no action required for t, v, w, r)
        dA.kd = +s/rhd;
    end
    %
    % Apply upstream boundary conditions
    if Ru == -1
        % homogenous dirichlet on full matrices
        dA.n(N+1,:)    = zeros(1,N+1); % bot row of dA.n    -> 0
        dA.n(:,N+1)    = zeros(N+1,1); % rit col of dA.n    -> 0
        % homogenous Dirichlet on split matrices
        dA.tl(N+1,:)   = zeros(1,N+1); % bot row of dA.tl   -> 0
        dA.tr(:,N+1)   = zeros(N+1,1); % rit col of dA.tr   -> 0
        dA.hl(N+1,:)   = zeros(1,N+1); % bot row of dA.hl   -> 0
        dA.hr(N+1)     = 0;            % rit col of dA.hr   -> 0
        dA.wrl(N+1)    = 0;            % bot row of dA.wrl  -> 0
        dA.wrr(:,N+1)  = zeros(N+1,1); % rit col of dA.wrr  -> 0
        dA.vl(N+1,:)   = zeros(1,N+1); % bot row of dA.vl   -> 0
        dA.vr(:,N+1)   = zeros(N+1,1); % rit col of dA.vr   -> 0
        dA.ku = 0;
        dA.mrul(N+1,:) = zeros(1,N+1); % bot row of dA.mrul -> 0
        dA.mrur(:,N+1) = zeros(N+1,1); % rit col of dA.mrur -> 0
        dA.mrpl(N+1,:) = zeros(1,N+1); % bot row of dA.mrpl -> 0
        dA.mrpr(:,N+1) = zeros(N+1,1); % rit col of dA.mrpr -> 0
        dA.frul(N+1,:) = zeros(1,N+1); % bot row of dA.frul -> 0
        dA.frur(:,N+1) = zeros(N+1,1); % rit col of dA.frur -> 0
        dA.frpl(N+1,:) = zeros(1,N+1); % bot row of dA.frpl -> 0
        dA.frpr(:,N+1) = zeros(N+1,1); % rit col of dA.frpr -> 0
        dA.qpul(N+1,:) = zeros(1,N+1); % bot row of dA.qpul -> 0
        dA.qpur(:,N+1) = zeros(N+1,1); % rit col of dA.qpur -> 0
        dA.qppl(N+1,:) = zeros(1,N+1); % bot row of dA.qppl -> 0
        dA.qppr(:,N+1) = zeros(N+1,1); % rit col of dA.qppr -> 0
    elseif Ru == +1
        % homogenous Neumann (no action required)
        dA.ku = 0;
    else
        % Robin (no action required for t, v, w, r)
        dA.ku = s/rhu;
    end
end

end
