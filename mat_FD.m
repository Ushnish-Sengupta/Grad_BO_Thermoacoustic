function [pos,mat] = mat_FD(param,N)
% mat_FD
%
% Generate matrices for the finite difference scheme
%
% INPUTS: 
% param     structure containing the parameters
% N         N+1 is the number of gridpoints
%
% OUTPUTS:
% pos.x     positions of the gridpoints x_i
% mat.D     first derivative matrix: dp/dx = D*p
% mat.M     mass matrix: <f,g> = f^H * M * g
% mat.V     matrix containing the specific volume v(x_i) along the diagonal
% mat.t     column vector containing time delay t(x_i)
% mat.h     column vector containing heat release envelope h(x_i) 
% mat.wr    column vector containing measurement envelope / density wr(x_i)
% mat.dwrdx column vector containing d(wr)/d(x) at x_i

%% Generate the points and the building block matrices
% Create collocation points from x = 0 to +1 and mass col vector, m
[x,m] = fun_clencurt(N);
% Create 1st order Chebyshev differentiation matrix from 0 to 1, D
D = fun_cheb(N);
% Create the mass matrix, M
M = diag(m);

%% Load in the density profile and create the density matrix, R
[rh,~,~,~] = fun_rh(param,x); V = diag(1./rh);

%% Load in the flame profile
% Generate the heat release envelope divided by pressure, vp(x), which integrates to 1
h = fun_h(param,x);
% Generate the measurement envelope divided by density, wr(x), and its derivative w.r.t. x 
[wr,dwrdx] = fun_wr(param,x);
% Generate the time delay vector (uniform in space for this problem)
t = param.tau*ones(N+1,1);

%% Wrap into mat structure
pos.x     = x;
mat.D     = D;
mat.M     = M;
mat.V     = V;
mat.t     = t;
mat.h     = h;
mat.wr    = wr;
mat.dwrdx = dwrdx;

end

function [D] = fun_cheb(N)
% CHEB  compute D = differentiation matrix, x = Chebyshev grid
% (adapted from Spectral Methods in Matlab by Nick Trefethen)

if N==0, D=0; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D,2));                 % diagonal entries
D  = D*2; % convert to points from x = 0 to 1

end

function [x,m] = fun_clencurt(N)
% CLENCURT  nodes x (Chebyshev points) and weights m for Clenshaw-Curtis quadrature
% (adapted from Spectral Methods in Matlab by Nick Trefethen)

theta = pi*(0:N)'/N; x = cos(theta);
m = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
if mod(N,2)==0
    m(1) = 1/(N^2-1); m(N+1) = m(1);
    for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    v = v - cos(N*theta(ii))/(N^2-1);
else
    m(1) = 1/N^2; m(N+1) = m(1);
    for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
end
m(ii) = 2*v/N;
x = (x+1.0)/2; % Convert to points from x = 0 to +1
m = m/2;       % Convert weights to points from x = 0 to +1

end

  