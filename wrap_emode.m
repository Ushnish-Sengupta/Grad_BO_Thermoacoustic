function [emode] = wrap_emode(pos,mat,T,pl,s,pr,discretization,type_of_adjoint)
% wrap_emode
%
% wrap eigenmode properties into a structure
%
% INPUTS
% pos               structure containing x or x1 and x0
% mat               structure containing the FD, FE, or SBP matrices
% T                 structure containing the transformation matrices
% pl                left eigenvector
% s                 eigenvalue
% pr                right eigenvector
% discretization    discretization (FDS, FEW, FDW, or SBP)
% type_of_adjoint   type of adjoint (DA or CA)
%
% OUTPUTS
% emode             structure containing the above quantities

% Wrap the direct or adjoint eigenvalue
emode.s   = s;

% Normalize pr and pl and rename T.* for FEW
switch discretization
    case 'FEW'
        pr = fun_normalize(pr,mat.M11); 
        pl = fun_normalize(pl,mat.M11);
        T.MR = T.MR11;
        T.FR = T.FR10;
        T.QP = T.QP11;
        T.U  = T.U01;
    case {'FDS','FDW','SBP'}
        pr = fun_normalize(pr,mat.M); 
        pl = fun_normalize(pl,mat.M);
end

% Wrap the discretization points
switch discretization
    case {'FEW','SBP'}
        emode.x1 = pos.x1;   
        emode.x0 = pos.x0;
    case {'FDS','FDW'}
        emode.x = pos.x;
end

% Wrap the eigenvectors and receptitivities
switch type_of_adjoint
    case 'DA'
        emode.pr = pr;               % right p eigenvector & direct  p eigenfunction
        emode.pl = pl;               % left  p eigenvector & adjoint p eigenfunction
        emode.ur = T.U * pr;         % right u eigenvector & direct  u eigenfunction
        emode.mr = pl' * T.MR;       % Discrete Receptivity of mass equation to injection of mdot/rho
        emode.fr = pl' * T.FR;       % Discrete Receptivity of momentum equation to injection of f/rho
        emode.qp = pl' * T.QP;       % Discrete Receptivity of energy equation to injection of q/p
    case 'CA'
        emode.pr = pr;               % right p^dag eigenvector & adjoint pressure eigenfunction
        emode.mr = pr' * T.MR;       % Continuous Receptivity of mass equation to injection of mdot/rho
        emode.fr = pr' * T.FR;       % Continuous Receptivity of momentum equation to injectrion of f/rho
        emode.qp = pr' * T.QP;       % Continuous Receptivity of energy equation to injection of q/p
end

% Wrap the mass matrices for DA 
switch type_of_adjoint
    case 'DA'
        switch discretization
            case 'FEW'
                emode.M11 = mat.M11; % Mass matrices, required to convert from Discrete to Continuous Receptivity
                emode.M00 = mat.M00; % Mass matrices, required to convert from Discrete to Continuous Receptivity
            case {'FDS','FDW'}
                emode.M   = mat.M;   % Mass matrix, required to convert from Discrete to Continuous Receptivity
            case 'SBP'
                emode.M   = mat.M;   % Mass matrix on x1 points, required to convert from Discrete to Continuous Receptivity
                emode.M00 = mat.M00; % Mass matrix on x0 points, required to convert from Discrete to Continuous Receptivity
        end
end

end

