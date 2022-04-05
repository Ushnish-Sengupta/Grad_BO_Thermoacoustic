function [A,C,dAds] = fun_bcs_weak(A,C,dAds,N,Ru,cu,kus,rhu,Rd,cd,kds,rhd,s)
% fun_bcs_weak
%
% Apply boundary conditions to A and C matrices for the weak form equations 
%
% INPUTS: 
% A     matrix representing acoustics and heat release without b.c.'s
% C     matrix representing 1/gamma without b.c.'s
% dAds  dA/ds
% N     N+1 is the number of gridpoints or the number of finite elements
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
    % Homogenous Dirichlet condition on the upstream boundary
    A(N+1,:)   = zeros(1,N+1); % set bott  row of A0  to zero
    A(:,N+1)   = zeros(N+1,1); % set right col of A0  to zero
    A(N+1,N+1) = cu;            % set bott-right element of A0 to cu
    dAds(N+1,:) = zeros(1,N+1);
    dAds(:,N+1) = zeros(N+1,1);
    dAds(N+1,N+1) = 0;
    C(N+1,:) = zeros(1,N+1);    % set bott  row of C  to zero
    C(:,N+1) = zeros(N+1,1);    % set right col of C  to zero
elseif Ru == +1
    % Homogenous Neumann condition on the upstream boundary
    % (no action required)
else
    % Robin condition on the upstream boundary dp/dx = - kus/rhu*s*p
    % (note that ku = s * kus and kus is not a function of s)
    A(N+1,N+1)    = A(N+1,N+1)    + s * kus / rhu;
    dAds(N+1,N+1) = dAds(N+1,N+1) +     kus / rhu;
end

if Rd == -1
    % Homogenous Dirichlet condition on the downstream boundary
    A(1,:)  = zeros(1,N+1); % set top  row of A0  to zero
    A(:,1)  = zeros(N+1,1); % set left col of A0  to zero
    A(1,1)  = cd;            % set top-left element of A0 to cd
    dAds(1,:) = zeros(1,N+1);
    dAds(:,1) = zeros(N+1,1);
    dAds(1,1) = 0;
    C(1,:)   = zeros(1,N+1); % set top  row of C  to zero
    C(:,1)   = zeros(N+1,1); % set left col of C  to zero
elseif Rd == +1
    % Homogenous Neumann condition on the downstream boundary
    % (no action required)
else
    % Robin condition on the downstream boundary dp/dx = kds/rhd*s*p
    % (note that kd = s * kds and kds is not a function of s)
    A(1,1)     = A(1,1)     + s * kds / rhd;
    dAds(1,1)  = dAds(1,1)  +     kds / rhd;
end

end



