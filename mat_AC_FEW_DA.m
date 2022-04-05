function [A11,C11,dA11ds,T,dA11,ds_int] = mat_AC_FEW_DA(mat,param,N,s)
% mat_AC_FEW_DA
%
% Generate A and C matrices for Discrete Adjoint Finite Element Weak form
%
% INPUTS: 
% mat.D01    first derivative matrix: (dp/dx)0 = D01*p1
% mat.M00    mass matrix: <f0,g0> = f0^H * M00 * g0
% mat.M01    mass matrix: <f0,g1> = f0^H * M01 * g1
% mat.M11    mass matrix: <f1,g1> = f1^H * M11 * g1
% mat.v0    matrix containing the specific volume v(x0_i) along the diagonal
% mat.t     column vector containing time delay t(x1_i)
% mat.h     column vector containing heat release envelope h(x1_i) 
% mat.wr    column vector containing measurement envelope / density wr(x0_i)
% mat.dwrdx1 column vector containing d(wr)/d(x) at x1_i
% param      structure containing the parameters
% N          N+1 is the number of gridpoints
% conjs      value of conj(s) at which to calculate the matrices
%
% OUTPUTS:
% A11        matrix representing acoustics and heat release, with boundary conditions applied 
% C11        matrix representing 1/gamma, with boundary conditions applied
% dA11/ds    dA/ds
% T.P11      matrix that maps p1 to p1
% T.U01      matrix that maps p1 to u0
% T.MR10     matrix that maps adjp1 to the Receptivity of the mass equation to injection of mdot/rho0
% T.MR11     matrix that maps adjp1 to the Receptivity of the mass equation to injection of mdot/rho1
% T.FR10     matrix that maps adjp1 to the Receptivity of the momentum equation to injectrion of f/rho0 
% T.QP10     matrix that maps adjp1 to the Receptivity of the energy equation to injection of q/p0
% T.QP11     matrix that maps adjp1 to the Receptivity of the energy equation to injection of q/p0
% dA11.*     structure containing all the sensitivities: dA11/d(param)
% ds_int.*   initialized structure for ds of all the internal parameters

%% Unwrap mat 
[D01,M11,M00,M01,V00,t,h,wr,~] = unwrap_mat_FE(mat);

%% Unwrap param
[gam,zet,n,Ru,Rd,rhu,rhd,kus,kds,cu,cd] = unwrap_param(param);

%% Generate A, C, and dAds matrices
A11    = - D01'*M00*V00*D01 - M11 * zet * n * diag(exp(-s*t))     * h * wr.' * M00 * D01;
C11    =   M11/gam;
dA11ds =                    - M11 * zet * n * diag(-t.*exp(-s*t)) * h * wr.' * M00 * D01;

%% Apply the boundary conditions to A, C, and dAds
[A11,C11,dA11ds] = fun_bcs_weak(A11,C11,dA11ds,N,Ru,cu,kus,rhu,Rd,cd,kds,rhd,s);

%% Define the conversion / transformation matrices
if nargout >= 4
    % Conversion from pr1 to pr1
    T.P11 = eye(N+1);
    % Conversion from pr1 to ur0
    T.U01 = - V00 * D01 / s;
    % Conversion from pl to Receptivity of mass equation to injection of mdot/rho
    T.MR10 = s * M01';
    T.MR11 = s * M11;
    % Conversion from pl to Receptivity of momentum equation to injectrion of f/rho
    T.FR10 = D01' * M00;
    % Conversion from pl to Receptivity of energy equation to injection of q/p
    T.QP10 = zet * s * M01';
    T.QP11 = zet * s * M11;
end

%% Generate all the sensitivities of matrix A: dA.* = dA.*l  x  d.*  x  dA.*r
if nargout >= 5
    % initialize the ds structure
    ds_int.n   = 0;
    ds_int.t   = zeros(1,N+1);
    ds_int.h   = zeros(1,N+1);
    ds_int.wr  = zeros(1,N);
    ds_int.v   = zeros(1,N);
    ds_int.ku  = 0;
    ds_int.kd  = 0;
    ds_int.mru = zeros(1,N);
    ds_int.mrp = zeros(1,N+1);
    ds_int.fru = zeros(1,N);
    ds_int.frp = zeros(1,N);
    ds_int.qpu = zeros(1,N);
    ds_int.qpp = zeros(1,N+1);
    
    % Generate sensitivities in the bulk (before boundary conditions)
    % Sensitivity to n
    dA11.n   = - M11 * zet         * diag(exp(-s*t)) * h * wr.' * M00 * D01;
    % Sensitivity to t1
    dA11.tl  = - M11;
    dA11.tr  =       - zet * n * s * diag(exp(-s*t)) * h * wr.' * M00 * D01;
    % Sensitivity to vp1
    dA11.hl  = - M11 * zet * n     * diag(exp(-s*t));
    dA11.hr  =                                             wr.' * M00 * D01;
    % Sensitivity to wr0
    dA11.wrl = - M11 * zet * n     * diag(exp(-s*t)) * h;
    dA11.wrr =                                                    M00 * D01;
    % Sensitivity to v0 (1/rho)
    dA11.vl  = - D01'*M00;
    dA11.vr  =             D01;
    
    % Sensitivity to feedback from u into mass/rho equation
    dA11.mrul = T.MR10;
    dA11.mrur = T.U01;
    % Sensitivity to feedback from p into mass/rho equation
    dA11.mrpl = T.MR11;
    dA11.mrpr = T.P11;
    % Sensitivity to feedback from u into momentum/rho equation
    dA11.frul = T.FR10;
    dA11.frur = T.U01;
    % Sensitivity to feedback from p into momentum/rho equation
    dA11.frpl = D01'; 
    dA11.frpr = M01;  
    % Sensitivity to feedback from u into energy/p equation
    dA11.qpul = T.QP10;
    dA11.qpur = T.U01;
    % Sensitivity to feedback from p into energy/p equation
    dA11.qppl = T.QP11;
    dA11.qppr = T.P11;
    
    % Apply downstream boundary conditions
    if Rd == -1
        % homogenous Dirichlet on full matrices
        dA11.n(1,:)      = zeros(1,N+1); % top row of dA11.n      -> 0
        dA11.n(:,1)      = zeros(N+1,1); % lef col of dA11.n      -> 0
        % homogenous dirichlet on split matrices
        dA11.tl(1,:)     = zeros(1,N+1); % top row of dA11.tl     -> 0
        dA11.tr(:,1)     = zeros(N+1,1); % lef col of dA11.tr     -> 0
        dA11.hl(1,:)     = zeros(1,N+1); % top row of dA11.hl     -> 0
        dA11.hr(1)       = 0;            % lef col of dA11.hr     -> 0
        dA11.wrl(1)      = 0;            % top row of dA11.wrl    -> 0
        dA11.wrr(:,1)    = zeros(N,1);   % lef col of dA11.wrr    -> 0
        dA11.vl(1,:)     = zeros(1,N);   % top row of dA11.vl     -> 0
        dA11.vr(:,1)     = zeros(N,1);   % lef col of dA11.vr     -> 0
        dA11.kd          = 0;
        dA11.mrul(1,:)   = zeros(1,N);   % top row of dA11.mrul10 -> 0
        dA11.mrur(:,1)   = zeros(N,1);   % lef col of dA11.mrur01 -> 0
        dA11.mrpl(1,:)   = zeros(1,N+1); % top row of dA11.mrpl11 -> 0
        dA11.mrpr(:,1)   = zeros(N+1,1); % lef col of dA11.mrpr11 -> 0
        dA11.frul(1,:)   = zeros(1,N);   % top row of dA11.frul10 -> 0
        dA11.frur(:,1)   = zeros(N,1);   % lef col of dA11.frur01 -> 0
        dA11.frpl(1,:)   = zeros(1,N);   % top row of dA11.frpl10 -> 0
        dA11.frpr(:,1)   = zeros(N,1);   % lef col of dA11.frpr01 -> 0
        dA11.qpul(1,:)   = zeros(1,N);   % top row of dA11.qpul10 -> 0
        dA11.qpur(:,1)   = zeros(N,1);   % lef col of dA11.qpur01 -> 0
        dA11.qppl(1,:)   = zeros(1,N+1); % top row of dA11.qppl11 -> 0
        dA11.qppr(:,1)   = zeros(N+1,1); % lef col of dA11.qppr11 -> 0
    elseif Rd == +1
        % homogenous neumann (no action required)
        dA11.kd = 0;
        % Note: the actual sensitivity to kd is +s/rhu (as for Robin b.c.'s)
        % but kd is irrelevant when Rd = 1 so I hard-wire this term to zero.
    else
        % Robin (no action required for t, v, w, r)
        dA11.kd = +s/rhd;
    end
    %
    % Apply upstream boundary conditions
    if Ru == -1
        % homogenous dirichlet on full matrices
        dA11.n(N+1,:)    = zeros(1,N+1); % bot row of dA11.n      -> 0
        dA11.n(:,N+1)    = zeros(N+1,1); % rit col of dA11.n      -> 0
        % homogenous Dirichlet on split matrices
        dA11.tl(N+1,:)   = zeros(1,N+1); % bot row of dA11.tl     -> 0
        dA11.tr(:,N+1)   = zeros(N+1,1); % rit col of dA11.tr     -> 0
        dA11.hl(N+1,:)   = zeros(1,N+1); % bot row of dA11.hl     -> 0
        dA11.hr(N+1)     = 0;            % rit col of dA11.hr     -> 0
        dA11.wrl(N+1)    = 0;            % bot row of dA11.wrl    -> 0
        dA11.wrr(:,N+1)  = zeros(N,1);   % rit col of dA11.wrr    -> 0
        dA11.vl(N+1,:)   = zeros(1,N);   % bot row of dA11.vl     -> 0
        dA11.vr(:,N+1)   = zeros(N,1);   % rit col of dA11.vr     -> 0
        dA11.ku          = 0;
        dA11.mrul(N+1,:) = zeros(1,N);   % bot row of dA11.mrul10 -> 0
        dA11.mrur(:,N+1) = zeros(N,1);   % rit col of dA11.mrur01 -> 0
        dA11.mrpl(N+1,:) = zeros(1,N+1); % bot row of dA11.mrpl11 -> 0
        dA11.mrpr(:,N+1) = zeros(N+1,1); % rit col of dA11.mrpr11 -> 0
        dA11.frul(N+1,:) = zeros(1,N);   % bot row of dA11.frul10 -> 0
        dA11.frur(:,N+1) = zeros(N,1);   % rit col of dA11.frur01 -> 0
        dA11.frpl(N+1,:) = zeros(1,N);   % bot row of dA11.frpl10 -> 0
        dA11.frpr(:,N+1) = zeros(N,1);   % rit col of dA11.frpr01 -> 0
        dA11.qpul(N+1,:) = zeros(1,N);   % bot row of dA11.qpul10 -> 0
        dA11.qpur(:,N+1) = zeros(N,1);   % rit col of dA11.qpur01 -> 0
        dA11.qppl(N+1,:) = zeros(1,N+1); % bot row of dA11.qppl11 -> 0
        dA11.qppr(:,N+1) = zeros(N+1,1); % rit col of dA11.qppr11 -> 0
    elseif Ru == +1
        % homogenous Neumann (no action required)
        dA11.ku = 0;
        % Note: the actual sensitivity to ku is s/rhu (as for Robin b.c.'s)
        % but ku is irrelevant when Ru = 1 so I hard-wire this term to zero.
    else
        % Robin (no action required for t, v, w, r)
        dA11.ku = s/rhu;
    end
end


end
