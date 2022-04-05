function [A11,C11,dA11ds,T] = mat_AC_FEW_CA(mat,param,N,conjs)
% mat_AC_FEW_CA
%
% Generate A and C matrices for Continuous Adjoint Finite Element
%
% INPUTS: 
% mat.D01    first derivative matrix: (dp/dx)0 = D01*p1
% mat.M00    mass matrix: <f0,g0> = f0^H * M00 * g0
% mat.M01    mass matrix: <f0,g1> = f0^H * M01 * g1
% mat.M11    mass matrix: <f1,g1> = f1^H * M11 * g1
% mat.V00    matrix containing the specific volume v(x0_i) along the diagonal
% mat.t1     column vector containing time delay t(x1_i)
% mat.h1     column vector containing heat release envelope h(x1_i) 
% mat.wr0    column vector containing measurement envelope / density wr(x0_i)
% mat.dwrdx1 column vector containing d(wr)/d(x) at x1_i
% param      structure containing the parameters
% N          N+1 is the number of gridpoints
% conjs      value of conj(s) at which to calculate the matrices
%
% OUTPUTS:
% A11        matrix representing acoustics and heat release, with boundary conditions applied 
% C11        matrix representing 1/gamma, with boundary conditions applied
% dA11/ds    dA/ds
% T.P01      matrix that maps p1 to p0
% T.P11      matrix that maps p1 to p1
% T.U01      matrix that maps p1 to u0
% T.MR10     matrix that maps adjp1 to the Receptivity of the mass equation to injection of mdot/rho0
% T.MR11     matrix that maps adjp1 to the Receptivity of the mass equation to injection of mdot/rho1
% T.FR10     matrix that maps adjp1 to the Receptivity of the momentum equation to injectrion of f/rho0 
% T.QP10     matrix that maps adjp1 to the Receptivity of the energy equation to injection of q/p0
% T.QP11     matrix that maps adjp1 to the Receptivity of the energy equation to injection of q/p0

%% Unwrap mat 
[D01,M11,M00,~,V00,t1,h1,~,dwrdx1] = unwrap_mat_FE(mat);

%% Unwrap param
[gam,zet,n,Ru,Rd,rhu,rhd,kus,kds,cu,cd] = unwrap_param(param);

%% Generate A11, C11, and dA11ds matrices
A11    = - D01'*M00*V00*D01 + M11 * zet * n * dwrdx1 * h1' * diag(exp(-conjs*conj(t1))) * M11;
C11    = M11/gam;
dA11ds = M11 * zet * n * dwrdx1 * h1' * diag(-conj(t1).*exp(-conjs*conj(t1))) * M11;

%% Apply the boundary conditions to A, C, and dAds
[A11,C11,dA11ds] = fun_bcs_weak(A11,C11,dA11ds,N,Ru,cu,conj(kus),rhu,Rd,cd,conj(kds),rhd,conjs); 

%% Define the conversion / transformation matrices
if nargout >= 4
    % Calculate the direct eigenvalue
    s = conj(conjs);
    % Define averaging matrix
    M01bar = eye(N+1) + diag(ones(1,N),+1); M01bar = M01bar(1:N,:)/2;
    % Conversion from pr to Receptivity of mass equation to injection of mdot/rho
    T.MR10 = s * M01bar';
    T.MR11 = s * eye(N+1);
    % Conversion from pr to Receptivity of momentum equation to injectrion of f/rho
    T.FR10 = D01';
    % Conversion from pr to Receptivity of energy equation to injection of q/p
    T.QP10 = zet * s * M01bar';
    T.QP11 = zet * s * eye(N+1);
    % Conversion from direct p1 mode to direct p1 mode (used for continuous sensitivities) 
    T.P11  = eye(N+1);
    % Conversion from direct p1 mode to direct p0 mode (used for continuous sensitivities) 
    T.P01  = M01bar;
    % Conversion from direct p1 mode to direct u0 mode (used for continuous sensitivities) 
    T.U01  = - V00 * D01 / s;
end

end
