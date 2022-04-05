function [pos,mat] = mat_FE(param,N)
% mat_FE
%
% Generate matrices for the finite element scheme
%
% INPUTS: 
% param      structure containing the parameters
% N          N+1 is the number of gridpoints
%
% OUTPUTS:
% pos.x0     positions of the gridpoints x0_i
% pos.x1     positions of the gridpoints x1_i
% mat.D01    first derivative matrix: (dp/dx)0 = D01*p1
% mat.M00    mass matrix: <f0,g0> = f0^H * M00 * g0
% mat.M01    mass matrix: <f0,g1> = f0^H * M01 * g1
% mat.M11    mass matrix: <f1,g1> = f1^H * M11 * g1
% mat.V00    matrix containing the specific volume v(x0_i) along the diagonal
% mat.t      column vector containing time delay t(x1_i)
% mat.h      column vector containing heat release envelope h(x1_i) 
% mat.wr     column vector containing measurement envelope / density wr(x0_i)
% mat.dwrdx  column vector containing d(wr)/d(x) at x1_i
%
% NOTES
% *0 are col vectors representing P0 functions. They contain N elements
% *1 are col vectors representing P1 functions. They contain N+1 elements
% *00 are matrices that map P0 to P0. They contain (N,N) elements
% *01 are matrices that map P1 to P0. They contain (N+1,N) elements
% *11 are matrices that map P1 to P1. They contain (N+1,N+1) elements

%% Generate the discretization points and the building block matrices
% Generate the positions of the gridpoints for P1 functions, x1
x1 = linspace(+1,0,N+1)'; dx = 1/N;
% Generate the positions of the gridpoints for P0 functions, x0
x0 = conv2(x1,[0.5;0.5],'valid');
% Generate 1st order difference matrix for a P1 function (makes a P0 function)
D01 = eye(N+1) - diag(ones(1,N),+1); D01 = D01(1:N,:); D01 = D01/dx;
% Generate the mass matrix for two P1 functions, M11
M11 = 4*eye(N+1) + diag(ones(1,N),+1) + diag(ones(1,N),-1); M11(1,1) = 2; M11(N+1,N+1) = 2; M11 = M11*(dx/6);
% Generate the mass matrix for two P0 functions, M00
M00 = eye(N)*dx;
% Generate the mass matrix for a P0 function * a P1 function, M01
M01 = eye(N+1) + diag(ones(1,N),+1); M01 = M01(1:N,:)/2*dx; 

%% Load in the density profile and create the density matrix, R00
[rh0,~,~,~] = fun_rh(param,x0); V00 = diag(1./rh0);

%% Load in the flame profile and generate the envelope functions
% Generate the heat release envelope divided by pressure, vp(x), which integrates to 1
h = fun_h(param,x1);
% Generate the measurement envelope divided by density, wr(x), and its derivative w.r.t. x  
[wr,~] = fun_wr(param,x0); [~,dwrdx] = fun_wr(param,x1);
% Generate the time delay vector (uniform in space for this problem)
t = param.tau*ones(N+1,1);

%% Wrap into mat structure
pos.x0     = x0;
pos.x1     = x1;
mat.D01    = D01;
mat.M11    = M11;
mat.M00    = M00;
mat.M01    = M01;
mat.V00    = V00;
mat.t      = t;
mat.h      = h;
mat.wr     = wr;
mat.dwrdx  = dwrdx;

end