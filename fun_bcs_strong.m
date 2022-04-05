function [A,C,dAds] = fun_bcs_strong(D,A,C,dAds,N,Ru,cu,kus,Rd,cd,kds,s)
% fun_bcs_strong
%
% Apply boundary conditions to A and C matrices for the strong form equations 
%
% INPUTS: 
% D     first derivative matrix: dp/dx = D*p
% A     matrix representing acoustics and heat release without b.c.'s
% C     matrix representing 1/gamma without b.c.'s
% dAds  dA/ds
% N     N+1 is the number of gridpoints
% Ru    upstream reflection coefficient (used here as a switch)
% cu    arbitrary weighting for upstream boundary condition
% kus   Robin upstream boundary coefficient is ku; kus = ku/s
% Rd    downstream reflection coefficient (used here as a switch)
% cd    arbitrary weighting for downstream boundary condition
% kds   Robin downstream boundary coefficient is kd; kds = kd/s
% s     The value of s at which to evaluate the boundary condition
%
% OUTPUTS:
% A     matrix A with boundary conditions applied
% C     matrix C with boundary conditions applied
% dAds  matrix dAds with boundary conditions applied

if Ru == -1
    % Homogenous Dirichlet condition on the upstream boundary (x = 0)
    A(N+1,:)      = zeros(1,N+1);
    A(:,N+1)      = zeros(N+1,1);
    A(N+1,N+1)    = cu;
    C(N+1,N+1)    = 0;
    dAds(N+1,:)   = zeros(1,N+1);
    dAds(:,N+1)   = zeros(N+1,1);
    dAds(N+1,N+1) = 0;
elseif Ru == +1
    % Homogenous Neumann condition on the upstream boundary
    A(N+1,:)      = D(N+1,:) * cu;
    C(N+1,N+1)    = 0;
    dAds(N+1,:)   = zeros(1,N+1);
else
    % Robin condition on the upstream boundary
    A(N+1,:)      = D(N+1,:) * cu;
    C(N+1,N+1)    = 0;
    dAds(N+1,:)   = zeros(1,N+1);
    % Robin coefficient cu*dp/dx = - cu*kus*s*p where kus = (Ru-1)/(Ru+1) * sqrt(rhu/gam);
    % (note that ku = s * kus and kus is not a function of s)
    A(N+1,N+1)    = A(N+1,N+1)    + s * kus * cu;
    dAds(N+1,N+1) = dAds(N+1,N+1) +     kus * cu;
end

if Rd == -1
    % Homogenous Dirichlet condition on the downstream boundary (x = +1)
    A(1,:)        = zeros(1,N+1);
    A(:,1)        = zeros(N+1,1);
    A(1,1)        = cd;
    dAds(1,:)     = zeros(1,N+1);
    dAds(:,1)     = zeros(N+1,1);
    dAds(1,1)     = 0;
    C(1,1)        = 0;
elseif Rd == +1
    % Homogenous Neumann condition on the downstream boundary
    A(1,:)        = D(1,:) * cd;
    C(1,1)        = 0;
    dAds(1,:)     = zeros(1,N+1);
else
    % Robin condition on the downstream boundary
    A(1,:)        = D(1,:) * cd;
    C(1,1)        = 0;
    % Robin coefficient (cd * dp/dx = + cd*kds*s*p) where kds = (Rd-1)/(Rd+1) * sqrt(rhd/gam); 
    % (note that kd = s * kds and kds is not a function of s)
    A(1,1)        = A(1,1)        - s * kds * cd;
    dAds(1,1)     = dAds(1,1)     -     kds * cd;
end

end
