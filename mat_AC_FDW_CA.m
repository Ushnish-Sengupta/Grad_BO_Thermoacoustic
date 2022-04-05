function [A,C,dAds,T] = mat_AC_FDW_CA(mat,param,N,conjs)
% mat_AC_FDW_CA
%
% Generate A and C matrices for Continuous Adjoint Finite Difference in the weak form 
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
% conjs     value of conj(s) at which to calculate the matrices
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

%% Unwrap mat 
[D,M,V,t,h,~,dwrdx] = unwrap_mat_FD(mat);

%% Unwrap param
[gam,zet,n,Ru,Rd,rhu,rhd,kus,kds,cu,cd] = unwrap_param(param);

%% Generate A, C, and dAds matrices
A    = - D'*M*V*D + M * zet * n * dwrdx * h' * diag(          exp(-conjs*conj(t))) * M;
C    = M/gam;
dAds =              M * zet * n * dwrdx * h' * diag(-conj(t).*exp(-conjs*conj(t))) * M;

%% Apply the boundary conditions to A, C, and dAds
[A,C,dAds] = fun_bcs_weak(A,C,dAds,N,Ru,cu,conj(kus),rhu,Rd,cd,conj(kds),rhd,conjs);

%% Define the conversion / transformation matrices
if nargout >= 4
    % Calculate the direct eigenvalue
    s = conj(conjs);
    % Conversion from direct p mode to direct p mode (used for continuous sensitivities) 
    T.P = eye(N+1);
    % Conversion from direct p mode to direct u mode (used for continuous sensitivities) 
    T.U = - V * D / s;
    % Conversion from pr to Receptivity of mass equation to injection of mdot/rho
    T.MR = s * eye(N+1);
    % Conversion from pr to Receptivity of momentum equation to injectrion of f/rho
    T.FR = D';
    % Conversion from pr to Receptivity of energy equation to injection of q/p
    T.QP = zet * s * eye(N+1);
end

end
