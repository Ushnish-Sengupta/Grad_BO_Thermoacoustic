function [A,C,dAds,T] = mat_AC_SBP_CA(mat,param,N,conjs)
% mat_AC_SBP_CA
%
% Generate A and C matrices for Continuous Adjoint Finite Difference in the
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

%% Unwrap mat 
[D1,D01,M,M00,~,V00,V11,S,Ed,Eu,I,t,h,~,dwrdx] = unwrap_mat_SBP(mat);

%% Unwrap param
[gam,zet,n,~,~,rhu,rhd,kus,kds,~,~] = unwrap_param(param);

%% Create a matrix with -v(1) at (1,1) and v(N+1) at (N+1,N+1)
Vbar = zeros(N+1); Vbar(1,1) = -1/rhd; Vbar(N+1,N+1) = 1/rhu;

%% Generate A, C, and dAds matrices
A    = - D01' * M00 * V00 * D01 ... % d/dx (v(x) * d/dx) in the bulk, in SBP format
       + Vbar*S...                  % d/dx (v(x) * d/dx) at the boundaries, in SBP format (these are akin to boundary conditions)
       + M * zet * n * dwrdx * h' * diag(          exp(-conjs*conj(t))) * M; % heat release term
C    =   M/gam;
dAds = + M * zet * n * dwrdx * h' * diag(-conj(t).*exp(-conjs*conj(t))) * M;

%% Apply the boundary conditions to A, C, and dAds
% (Neumann is obtained by setting kus or kds to zero)
% (Dirichlet is obtained by setting kus or kds to a large number)
% This is a computationally inefficient way to implement the boundary
% conditions (it alters only 6 elements of A and 2 elements of dAds) but it
% is easy to follow and allows extension to higher order difference
% operators, as in Appendix B of Nystrand's thesis.

% Robin condition on the upstream boundary
A    = A    - Eu*(-(conjs * conj(kus))/rhu*I+Vbar*S); % Application of SAT b'conds
dAds = dAds - Eu*(-(        conj(kus))/rhu*I       );

% Robin condition on the downstream boundary
A    = A    - Ed*(-(conjs * conj(kds))/rhd*I+Vbar*S); % Application of SAT b'conds
dAds = dAds - Ed*(-(        conj(kds))/rhd*I       );

%% Define the conversion / transformation matrices
if nargout >= 4
    % Calculate the direct eigenvalue
    s = conj(conjs);
    % Conversion from direct p1 mode to direct p1 mode (used for continuous sensitivities) 
    T.P  = eye(N+1);
    % Conversion from direct p1 mode to direct u0 mode (used for continuous sensitivities) 
    T.U  = - V11 * D1 / s;
    % Conversion from pr to Receptivity of mass equation to injection of mdot/rho
    T.MR = s * eye(N+1);
    % Conversion from pr to Receptivity of momentum equation to injectrion of f/rho
    T.FR = D1';
    % Conversion from pr to Receptivity of energy equation to injection of q/p
    T.QP = zet * s * eye(N+1);
end

end
