function [A,C,dAds,T,dA,ds_int] = mat_AC_SBP_DA(mat,param,N,s)
% mat_AC_SBP_DA
%
% Generate A and C matrices for Discrete Adjoint Finite Difference in the
% strong form in the Summation by Parts formulation with the Simultaneous
% Approximation Term. 
%
% INPUTS: 
% mat.D1     first derivative matrix: dp/dx = D1*p
% mat.D2     second derivative matrix: d2p/dx2 = D2*p
% mat.D2v    second derivative matrix: d/dx(v * d/dx)p = D2v*p
% mat.M      mass matrix: <f,g> = f^H * M * g
% mat.BS     boundary derivative operator
% mat.E0     matrix labelling the boundary at x(1)
% mat.EN     matrix labelling the boundary at x(N+1)
% mat.V      matrix containing the specific volume v(x_i) along the diagonal
% mat.t      column vector containing time delay t(x_i)
% mat.h      column vector containing heat release envelope h(x_i) 
% mat.wr     column vector containing measurement envelope / density wr(x_i)
% mat.dwrdx  column vector containing d(wr)/d(x) at x_i
% param      structure containing the parameters
% N          N+1 is the number of gridpoints
% s          value of s at which to calculate the matrices
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
[D1,D01,M,M00,M01,V00,V11,S,Ed,Eu,I,t,h,wr,~] = unwrap_mat_SBP(mat);

%% Unwrap param
[gam,zet,n,~,~,rhu,rhd,kus,kds,~,~] = unwrap_param(param);

%% Create a matrix with -v(1) at (1,1) and v(N+1) at (N+1,N+1)
Vtilde = zeros(N+1); Vtilde(1,1) = -1/rhd; Vtilde(N+1,N+1) = 1/rhu;

%% Generate A, C, and dAds matrices
A    = - D01' * M00 * V00 * D01 ...   % d/dx (v(x) * d/dx) in the bulk, in SBP format
       + Vtilde*S ...                 % d/dx (v(x) * d/dx) at the boundaries, in SBP format (these are akin to boundary conditions)
       - M * zet * n * diag(exp(-s*t))     * h * wr.' * M * D1; % heat release term
C    =   M/gam;
dAds = - M * zet * n * diag(-t.*exp(-s*t)) * h * wr.' * M * D1;

%% Apply boundary conditions to A, C, and dAds
% (Neumann is obtained by setting kus or kds to zero)
% (Dirichlet is obtained by setting kus or kds to a large number*)
% This is a computationally inefficient way to implement the boundary
% conditions (it alters only 6 elements of A and 2 elements of dAds) but it
% is easy to follow and allows extension to higher order difference
% operators, as in Appendix B of Nystrand's thesis.
% 
% * see Eq. (3.8) of Nystrand's thesis for a weak implementation of the
% Dirichlet boundary condition. When |sig_0| >> 1 and |sig_1| >> 1, 
% Eq. (3.8) seems to be equivalent to a Robin boundary condition with 
% |s * kus| >> 1 and |s * kds| >> 1. Therefore it is sensible to
% implement Dirichlet simply by setting kus and kds to large numbers. 

% Robin/Neumann condition on the upstream boundary
A    = A    - Eu*(-(s * kus)/rhu*I+Vtilde*S); % Application of SAT b'conds
dAds = dAds - Eu*(-(    kus)/rhu*I   );

% Robin/Neumann condition on the downstream boundary
A    = A    - Ed*(-(s * kds)/rhd*I+Vtilde*S); % Application of SAT b'conds
dAds = dAds - Ed*(-(    kds)/rhd*I   );

% % Neumann condition on the upstream boundary
% A    = A    - Eu*(Vbar*S); % Application of SAT b'conds
% % Neumann condition on the downstream boundary
% A    = A    - Ed*(Vbar*S); % Application of SAT b'conds
% % Dirichlet condition on the upstream boundary
% kus = 1000; A = A - (Eu*Vbar*S) - kus * Eu/rhu;
% % Dirichlet condition on the downstream boundary
% kds = -1000; A = A - (Ed*Vbar*S) + kds * Ed/rhu;

%% Define the conversion / transformation matrices
if nargout >= 4
    % Conversion from pr to pr
    T.P  = eye(N+1);
    % Conversion from pr to ur
    T.U  = - V11 * D1 / s;
    % Conversion from pl to Receptivity of mass equation to injection of mdot/rho
    T.MR = s * M;
    % Conversion from pl to Receptivity of momentum equation to injection of f/rho
    T.FR = D01' * M00;
    % Conversion from pl to Receptivity of energy equation to injection of q/p
    T.QP = s * zet * M;
end

%% Generate all the sensitivities of matrix A: dA.* = dA.*l  x  d.*  x  dA.*r
if nargout >= 5
    % initialize the ds structure
    ds_int.n   = 0;
    ds_int.t   = zeros(1,N+1);
    ds_int.h   = zeros(1,N+1);
    ds_int.wr  = zeros(1,N+1);
    ds_int.v   = zeros(1,N);
    ds_int.ku  = 0;
    ds_int.kd  = 0;
    ds_int.mru = zeros(1,N+1);
    ds_int.mrp = zeros(1,N+1);
    ds_int.fru = zeros(1,N);
    ds_int.frp = zeros(1,N);
    ds_int.qpu = zeros(1,N+1);
    ds_int.qpp = zeros(1,N+1);
    
    % Generate sensitivities in the bulk (before boundary conditions)
    % Sensitivity to n
    dA.n   = - M * zet         * diag(exp(-s*t)) * h * wr.' * M * D1;
    % Sensitivity to t
    dA.tl  = - M;
    dA.tr  =     - zet * n * s * diag(exp(-s*t)) * h * wr.' * M * D1;
    % Sensitivity to h
    dA.hl  = - M * zet * n     * diag(exp(-s*t));
    dA.hr  =                                           wr.' * M * D1;
    % Sensitivity to wr
    dA.wrl = - M * zet * n     * diag(exp(-s*t)) * h;
    dA.wrr =                                                  M * D1;
    % Sensitivity to v (1/rho)
    dA.vl  = -D01' * M00;
    dA.vr  =  D01;

    % Sensitivity to feedback from u into mass/rho equation
    dA.mrul = T.MR;
    dA.mrur = T.U;
    % Sensitivity to feedback from p into mass/rho equation
    dA.mrpl = T.MR;
    dA.mrpr = T.P;
    % Sensitivity to feedback from u into momentum/rho equation
    dA.frul = T.FR;
    dA.frur = - V00 * D01 / s;
    % Sensitivity to feedback from p into momentum/rho equation
    dA.frpl = D01';
    dA.frpr = M01; 
    % Sensitivity to feedback from u into energy/p equation
    dA.qpul = T.QP;
    dA.qpur = T.U;
    % Sensitivity to feedback from p into energy/p equation
    dA.qppl = T.QP;
    dA.qppr = T.P;

    % Apply downstream Robin boundary conditions
    dA.kd = s/rhd;
    
    % Apply upstream Robin boundary conditions
    dA.ku = s/rhu;
end
