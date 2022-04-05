function [kus,kds,dkusdRu,dkdsdRd] = fun_kukd(param)
% fun_kukd
%
% Calculate the Robin boundary coefficients, ku and kd
%
% The Robin boundary coefficients, ku and kd, are related to the reflection
% coefficients, Ru and Rd, by:
%
% k = s*(rho/gamma)^0.5 * (R-1)/(R+1)
%
% This function returns ku/s and kd/s given Ru and Rd. These are constants
% (determined only by Ru and Rd) and are multiplied by s when the boundary
% conditions are applied in fun_bcs_* 
%
% INPUTS
% param.Ru  upstream reflection coefficient
% param.Rd  downstream reflection coefficient
% param.rhu upstream density
% param.rhd downstream density
% param.gam ratio of the specific heat capacities
%
% OUTPUTS
% kus       the upstream Robin boundary coefficient is ku, where kus = ku/s
% kds       the downstream Robin boundary coefficient is kd, where kds = kd/s
% dkusdRu   d(kus)/d(Ru)
% dkdsdRd   d(kds)/d(Rd)

if param.Ru == -1 
    % Dirichlet boundary condition. kus is irrelevant except for SBP
    kus = 1e12; % set to a large number for SBP-SAT b'conds
    dkusdRu = 0;
elseif param.Ru == +1
    % Neumann boundary condition. kus is irrelevant except for SBP
    kus = 0;
    dkusdRu = 0;
else
    kus = (param.Ru-1)/(param.Ru+1) * sqrt(param.rhu/param.gam);
    dkusdRu    = 2 / (param.Ru+1)^2 * sqrt(param.rhu/param.gam);
end

if param.Rd == -1
    % Dirichlet boundary condition. kds is irrelevant except for SBP
    kds = 1e12; % set to a large number for SBP-SAT b'conds
    dkdsdRd = 0;
elseif param.Rd == +1
    % Neumann boundary condition. kds is irrelevant except for SBP
    kds = 0;
    dkdsdRd = 0;
else
    kds = (param.Rd-1)/(param.Rd+1) * sqrt(param.rhd/param.gam);
    dkdsdRd    = 2 / (param.Rd+1)^2 * sqrt(param.rhd/param.gam);
end

end
