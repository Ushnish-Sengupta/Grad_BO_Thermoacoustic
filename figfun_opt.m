function [pv,sv] = figfun_opt(discretization,type_of_adjoint,lin_or_nonlin,vary,goal,maxstep,NN)

% A crude optimization algorithm using base state sensitivites
%
% discretization:
% FD        finite difference of the strong form
% FE        finite element of the weak form
% WD        finite difference of the weak form
% 
% type_of_adjoint:
% DA        discrete adjoint
% CA        continuous adjoint
%
% lin_or_nonlin:
% lin       solve as a series of linear eigenvalue problems
% nonlin    solve as a nonlinear eigenvalue problem using a Newton method
%
% vary:
% X_w       measurement point
% X_h       heat release point
% heater    measurement and heat release points simultaneously
%
% goal:
% maX_hreq
% min_freq
% max_growth
% min_growth
%
% maxstep   maximum step size
%
% NN        Number of steps in optimization

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');
% Move the heater for the purpose of this demonstration
param_dim.X_w    = 0.25;    % position of measurement zone (metres)
param_dim.X_h    = 0.35;    % position of heat release region (metres)

%% Calculate the reference scales and the nondimensional parameters
param = fun_nondim(param_dim);

%% Set the numerical scheme, the max number of iterations, and the starting value of s
scheme.N     = 101;
scheme.itmax = 10;
scheme.s0    = fun_set_s0(param);
switch type_of_adjoint
    case 'DA'
        % Discrete adjoint; no action required
    case 'CA'
        % Continuous adjoint; start from the c.c. of s0
        scheme.s0 = conj(scheme.s0);
end

% Set the direction to step in
switch goal
    case 'max_growth'
        direction = +1;
    case 'min_growth'
        direction = -1;
    case 'maX_hreq'
        direction = -1i;
    case 'min_freq'
        direction = +1i;
    otherwise
        disp('goal not recognized')
end

% figure(1); clf; hold on

%% Pre-allocate
sv = zeros(NN,1);
pv = zeros(NN,1);

%% Calculate the gradients and move in the desired direction NN times
for nn = 1:NN
    % Calculate s and d(s)/d(all parameters)
    [emode,~,ds_ext] = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme);
    switch type_of_adjoint
        case 'DA'
            s = emode.s;
        case 'CA'
            s = conj(emode.s);
    end
    sv(nn) = s;
    % Update scheme.s0
    scheme.s0 = s;
    switch vary
        case 'X_w'
            % Vary the position of the measurement envelope
            param.X_w = param.X_w + maxstep*tanh(real(direction*(ds_ext.X_w)));
            disp(['X_w = ',num2str(param.X_w),', s = ',num2str(s)])
            pv(nn) = param.X_w;
        case 'X_h'
            % Vary the position of the heat release envelope
            param.X_h = param.X_h + maxstep*tanh(real(direction*(ds_ext.X_h)));
            disp(['X_h = ',num2str(param.X_h),', s = ',num2str(s)])
            pv(nn) = param.X_h;
        case 'heater'
            % Vary the positions of the measurement envelope and heat release envelope
            % (keeping their relative positions the same)
            param.X_w = param.X_w + maxstep*tanh(real(direction*((ds_ext.X_h + ds_ext.X_w)/2)));
            param.X_h = param.X_h + maxstep*tanh(real(direction*((ds_ext.X_h + ds_ext.X_w)/2)));
            disp(['heater posn = ',num2str((param.X_h + param.X_w)/2),', s = ',num2str(s)])
            pv(nn) = (param.X_h + param.X_w)/2;
        otherwise
            disp(['The script cannot optimize for vary = ',vary])
            return
    end
end

end
