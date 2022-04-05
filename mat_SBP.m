function [pos,mat] = mat_SBP(param,N)
% mat_SBP
%
% Generate matrices for the finite difference summation by parts scheme
% with a simultaneous approximation term for the boundary condition penalty
%
% This version can handle a non-uniform density. See Appendix B of
% Nystrand's 2014 thesis. In Nystrand's thesis, the spatially-varying
% diffusion coefficient is called b. Its equivalent here is the specific
% volume, v.
%
% In this version, the (N+1,N+1) matrix representing d/dx ( v(x) * d/dx) is
% given by D01' * V0 * D01, where D01 is a (N,N+1) matrix that is identical
% to that used in the FE method. V0 is a (N,N) matrix containing the
% specific volume along the diagonal. V0 can be calculated in two ways,
% which are nearly identical:
% (i) v_i = 0.5*v(x_i) + 0.5*v(x_(i+1)); i.e. the mean of the density
% evaluated at the gridpoints. This gives a matrix representing 
% d/dx (v(x) * d/dx) that is identical to that in Nystrand's thesis.
% (ii) v_i = v(0.5*x_i + 0.5*x_(i+1)); i.e. the density evaluated at the
% mean of the gridpoints. This gives a matrix representing 
% d/dx (v(x) * d/dx) that is identical to that in the FE method, apart from
% the boundary conditions. 
%
% The M00 matrix contains dx along the diagonal. It is included so that D01
% can be divided by dx (as it is in the FE method). This highlights the
% similarity between the SBP-SAT method and the FE method for second order
% central differencing. 
%
% INPUTS: 
% param      structure containing the parameters
% N          N+1 is the number of gridpoints
%
% OUTPUTS:
% pos.x1     positions of the gridpoints x1_i
% pos.x0     positions of the gridpoints x0_i (mid-way between x1_i)
% mat.D1     first derivative matrix: dp/dx = D1*p
% mat.D01    first order difference matrix, used to create d/dx(v*d/dx)
% mat.M      mass matrix for x1 points: <f1,g1> = f1^H * M * g1 (not the same as M11 in FE) 
% mat.M00    mass matrix for x0 points: <f0,g0> = f0^H * M00 * g0 (the same as M00 in FE) 
% mat.M01    mass matrix for x0*x1 pts: <f0,g1> = f0^H * M01 * g1 (the same as M01 in FE) 
% mat.V00    diagonal matrix containing specific volume at x0 points: v(x0)
% mat.V11    diagonal matrix containing specific volume at x1 points: v(x1)
% mat.S      second order accurate first derivative operator at the boundaries
% mat.Ed     matrix that labels the downstream boundary
% mat.Eu     matrix that labels the upstream boundary 
% mat.I      (N+1,N+1) identity matrix
% mat.t      column vector containing time delay t(x_i)
% mat.h      column vector containing heat release envelope h(x_i) 
% mat.wr     column vector containing measurement envelope / density wr(x_i)
% mat.dwrdx  column vector containing d(wr)/d(x) at x_i

%% Generate the positions of the gridpoints, x
x1 = linspace(+1,0,N+1)'; dx = 1/N;

%% Generate the mass matrix, M (known as H in Nystrand's thesis)
m = ones(N+1,1); m(1) = 0.5; m(N+1) = 0.5; m = m*dx; M = diag(m);

%% Generate the homogenous D1 and S matrices (D2 matrix commented out)
% Generate 2nd order accurate central difference 1st derivative matrix, D1
D1 = diag(0.5*ones(1,N),-1) + diag(-0.5*ones(1,N),+1); D1(1,1:2) = [1,-1]; D1(N+1,N:N+1) = [1 -1]; D1 = D1/dx;
% One can check that D1 is a SBP 1st derivative operator. 1st derivative 
% SBP operators are given by D1 = inv(M)*Q, where Q + Q' = B and 
% B = diag(+1,0,...,0,-1) when x is arranged from +1 to 0.
% Q = M*D1; B = Q + Q'; % Check that B(1,1) = 1, B(N+1,N+1) = -1; B(else,else) = 0
% Generate 2nd order accurate central difference 2nd derivative matrix, D2
% (This is never used, but it is instructive to leave it in the comments.)
% D2 = diag(ones(1,N),-1) + diag(ones(1,N),+1) + diag(-2*ones(1,N+1)); D2(1,:) = D2(2,:); D2(N+1,:) = D2(N,:); D2 = D2/dx^2;
% Generate the S matrix, which approximates the 1st derivative operator at the boundaries
% This differs from Nystrand's S matrix, in that it contains 0 rather than
% 1 along the leading diagonal, except at the top-left and bottom-right
% elements. This difference does not matter because S is always multiplied
% by a matrix containing 0 at all points except the top-left and
% bottom-right elements.  
S = zeros(N+1); S(1,1:3) = [-1.5 2 -0.5]; S(N+1,N-1:N+1) = [0.5 -2 1.5]; S = S/dx;

%% Generate the building block matrices for the inhomogenous D2 matrix
% The inhomogenous D2 matrix is labelled D2v and represents 
% d/dx ( v(x) * d/dx ). It is given by D2v = inv(M) * (-Mv + Vbar*S),
% where Vbar = diag(-v(1),0,...,0,v(N+1)), and Mv is Nystrand's M^(b)
% matrix. Mv can be coded directly (see comments at the end of this file)
% or it can be created from D01 and V0 via Mv = D01' * M00 * V0 * D01,
% where M00 is dx*eye(N) and V0 contains the mean density between adjacent
% pairs of gridpoints. V0 can be created in two different ways, which both
% give nearly identical results (see below). 

% Generate a (N,N) matrix with dx along the diagonal
M00 = eye(N)*dx;

% Generate the 1st order difference matrix that generates the gradient (at
% the x0 points) of functions evaluated at the x1 points (same as FE)
D01 = eye(N+1) - diag(ones(1,N),+1); D01 = D01(1:N,:); D01 = D01/dx;

% Choose how to evaluate the density at the centrepoints
density_method = 'centrepoint'; % 'mean' || 'centrepoint'
switch density_method
    case 'mean'
        % This is the version for Mv in Nystrand's thesis. The values of
        % the density are evaluated at the x1 gridpoints. Then the means of
        % adjacent values are used in the calculation.
        % Load in the density profile and create the specific volume
        % vector, v1, at the gridpoints x1.
        [rh1,~,~,~] = fun_rh(param,x1); v1 = (1./rh1); V11 = diag(v1);
        % Find the mean of adjacent values of v1
        v0 = conv2(v1,[0.5;0.5],'valid'); V00 = diag(v0);
    case 'centrepoint'
        % This is the version most like FE, and is almost identical to Mv
        % in Nystrand's thesis. The values are the density are evaluated at
        % the x0 gridpoints, which are the mid-points between the x1
        % gridpoints.
        % Generate the positions of the gridpoints for P0 functions, x0
        x0 = conv2(x1,[0.5;0.5],'valid');
        % Load in the density profile at the x0 points
        [rh0,~,~,~] = fun_rh(param,x0); V00 = diag(1./rh0);
        % Load in the denstiy profile at the x1 points; required for U(P)
        [rh1,~,~,~] = fun_rh(param,x1); V11 = diag(1./rh1);
end

% Generate the mass matrix for a P0 function * a P1 function, M01
M01 = eye(N+1) + diag(ones(1,N),+1); M01 = M01(1:N,:)/2*dx; 

%% Create the vectors and matrices that are zero except at the boundaries
% Create the downstream boundary label vector
ed = zeros(1,N+1); ed(1)   = 1; Ed = diag(ed);
% Create the upstream boundary label vector 
eu = zeros(1,N+1); eu(N+1) = 1; Eu = diag(eu);

%% Load in the flame profile and generate the envelope functions
% Generate the heat release envelope divided by pressure, vp(x), which integrates to 1
h = fun_h(param,x1);
% Generate the measurement envelope divided by density, wr(x), and its derivative w.r.t. x  
[wr,dwrdx] = fun_wr(param,x1);
% Generate the time delay vector (uniform in space for this problem)
t = param.tau*ones(N+1,1);

%% Wrap into mat structure
pos.x1     = x1;  % x1 gridpoints (same as FE)
pos.x0     = x0;  % x0 gridpoints (same as FE); midpoints between x1 gridpoints
mat.D1     = D1;  % d/dx from x1 points to x1 points
mat.D01    = D01; % d/dx from x1 points to x0 points
mat.M      = M;   % mass matrix for x1 points (n.b. not the same as M11 in FE)
mat.M00    = M00; % mass matrix for x0 points (same as M00 in FE)
mat.M01    = M01;
mat.V00    = V00; % specific volume at x0 points: v(x0)
mat.V11    = V11; % specific volume at x1 points: v(x1)
mat.S      = S;   % second order accurate first derivative operator at the boundaries
mat.Ed     = Ed;  % matrix that labels the downstream boundary
mat.Eu     = Eu;  % matrix that labels the upstream boundary 
mat.I      = eye(N+1);
mat.t      = t;
mat.h      = h;
mat.wr     = wr;
mat.dwrdx  = dwrdx;

end

% Direct implementation of Nystrand's thesis (Appendix B)
%
% S = eye(N+1); S(1,1:3) = [-1.5 2 -0.5]; S(N+1,N-1:N+1) = [0.5 -2 1.5]; S = S/dx;
%
%         % Create the left boundary closure of M^(b), given by a 3x3 matrix
%         Ml = [ 0.5*v(1) + 0.5*v(2) , -0.5*v(1) - 0.5*v(2)        ,  0                          ; ...
%             -0.5*v(1) - 0.5*v(2) ,  0.5*v(1) + v(2) + 0.5*v(3) , -0.5*v(2) - 0.5*v(3)        ; ...
%             0                   , -0.5*v(2) - 0.5*v(3)        ,  0.5*v(2) + v(3) + 0.5*v(4) ];
%         
%         % Create the right boundary closure of M^(b), given by a 3x3 matrix
%         Mr = [ 0.5*v(N) + v(N-1) + 0.5*v(N-2), -0.5*v(N) - 0.5*v(N-1)          ,  0                         ; ...
%             -0.5*v(N) - 0.5*v(N-1)         ,  0.5*v(N+1) + v(N) + 0.5*v(N-1) , -0.5*v(N+1) - 0.5*v(N) ; ...
%             0                             , -0.5*v(N+1) - 0.5*v(N)          ,  0.5*v(N+1) + 0.5*v(N)  ];
%         
%         % Create the interior stencil of M^(b), given by a row vector
%         Mi = zeros(N+1);
%         for nn = 2:N % n.b. Nystrand makes a mistake in writing 4:N-3 here
%             Mi(nn,nn-1) = -0.5*v(nn-1) - 0.5*v(nn)              ;
%             Mi(nn,nn)   =  0.5*v(nn-1) +     v(nn) + 0.5*v(nn+1);
%             Mi(nn,nn+1) =              - 0.5*v(nn) - 0.5*v(nn+1);
%         end
%         
%         % Construct the M^(b) matrix
%         Mb = Mi; Mb(1:3,1:3) = Ml ; Mb(N-1:N+1,N-1:N+1) = Mr;
%         % Nystrand forgot to divide by dx, but this is required. 
%         Mb = Mb/dx;
%         
%         % Construct the D2(b) operator; D2b = inv(H) * (-Mb + Vbar*S);
%         D2v = M\(-Mb + Vbar*S);
%
% Create the bar(B) matrix (known as bar(B) in Nystrand's thesis
% Vbar = zeros(N+1); Vbar(1,1) = -v1(1); Vbar(N+1,N+1) = v1(N+1);
