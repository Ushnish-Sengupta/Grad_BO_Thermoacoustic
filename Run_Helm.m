function [] = Run_Helm(discretization,type_of_adjoint,lin_or_nonlin,plot_type,TT_int,TT_ext)
% Run_Helm
% 
% call fun_Helm and plot requested results
%
% discretization:
% FDS       finite difference of the strong form
% FEW       finite element of the weak form
% FDW       finite difference of the weak form
% SBP       finite difference of the strong form in Summation by Parts form with the Simultaneous Approximation Term
% 
% type_of_adjoint:
% DA        discrete adjoint
% CA        continuous adjoint
%
% lin_or_nonlin
% lin       solve as a series of linear eigenvalue problems
% nonlin    solve as a nonlinear eigenvalue problem using a Newton method
%
% plot_type (optional)
% 
% TT_int (optional) Perform a Taylor Test on internal parameters
% {}        none
% {'all'}   all internal parameters {'t','h','wr','v','ku','kd'}
% {'*'}     user-defined internal parameters, e.g. {'t','h'}
%
% TT_ext (optional) Perform a Taylor Test on external parameters
% {}        none
% {'all'}   all external parameters {'n','X_w','L_w','X_h','X_h','tau','Ru','Rd'};
% {'*'}     user-defined external parameters, e.g. {'n','X_w'}
%
% INTERNAL PARAMETERS
% t         time delay field
% h         flame envelope field
% wr        measurement field divided by rho
% v         specific volume field
% ku        upstream Robin boundary coefficient
% kd        downstream Robin boundary coefficient
% all_int   all of the above simultaneously
%
% EXTERNAL PARAMETERS
% X_h etc.

if nargin == 3; plot_type = ''; TT_int = {}; TT_ext = {}; end
if nargin == 4;                 TT_int = {}; TT_ext = {}; end

%% Set the dimensional parameters
param_dim = fun_param_dim;

%% Calculate the reference scales and the nondimensional parameters
[param,ref] = fun_nondim(param_dim);

%% Set the numerical scheme, the max number of iterations, and the starting value of s
scheme.N     = 100;
scheme.itmax = 10;
scheme.s0    = fun_set_s0(param);
switch type_of_adjoint
    case 'CA'
        % Continuous adjoint; start from the c.c. of s0
        scheme.s0 = conj(scheme.s0);
end

%% Set figure number
switch discretization
    case 'FDS'
        fignum = 100;
    case 'FEW'
        fignum = 200;
    case 'FDW'
        fignum = 300;
    case 'SBP'
        fignum = 400;
    otherwise
        disp('discretization should be FDS, FEW, FDW, or SBP'); beep
        return
end

switch type_of_adjoint
    case 'DA'
        fignum = fignum + 10;
    case 'CA'
        fignum = fignum + 20;
end

switch lin_or_nonlin
    case {'nonlin','nonlinear'}
        fignum = fignum + 1;
    case 'linear'
        fignum = fignum + 2;
end

%% Calculate the eigenvalue with the linear or nonlinear method
switch plot_type
    case {'','emode','rec'}
        if nargin == 6 % TT requested
            [emode,~,~] = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme,TT_int,TT_ext);
        else
             emode = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme);
        end
    case {'bs_sens','fb_sens'}
        if nargin == 6 % TT requested
            [emode,ds_int,ds_ext] = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme,TT_int,TT_ext);
        else
            [emode,ds_int,ds_ext] = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme);
        end
end

%% Plots
switch plot_type
    case ''
        % Exit without plotting
    case 'emode'
        % Plot the direct and adjoint eigenvectors
        figure(fignum); clf
        switch type_of_adjoint
            case 'DA'
                subplot(3,1,1); hold on; ylabel('p_r')
                title([discretization,'\_',type_of_adjoint,'\_',lin_or_nonlin,'; eigenvalue s = ',num2str(emode.s)])
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x ,emode.pr)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x1,emode.pr)
                end
                subplot(3,1,2); hold on; ylabel('u_r')
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x ,emode.ur)
                    case 'FEW'
                        fun_plot_fun(emode.x0,emode.ur)
                    case 'SBP'
                        fun_plot_fun(emode.x1,emode.ur)
                end
                subplot(3,1,3); hold on; ylabel('p_l'); xlabel('x')
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x ,emode.pl)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x1,emode.pl)
                end
            case 'CA'
                subplot(3,1,3); hold on; ylabel('$p^\dag_r$','Interpreter','Latex'); xlabel('x')
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,emode.pr)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x1,emode.pr)
                end
        end
    case 'bs_sens'
        % Write some of the base state sensitivities to screen
        disp(['Base state sensitivity to n  = ',num2str(ds_int.n,8)])
        disp(['Base state sensitivity to ku = ',num2str(ds_int.ku,8)])
        disp(['Base state sensitivity to kd = ',num2str(ds_int.kd,8)])
        % Plot the base state sensitivities
        figure(fignum)
        clf
        subplot(4,1,1); hold on
        title([discretization,'\_',type_of_adjoint,'\_',lin_or_nonlin,'; eigenvalue s = ',num2str(emode.s)])
        ylabel('$\partial s/\partial t$','Interpreter','Latex')
        subplot(4,1,2); hold on
        ylabel('$\partial s/\partial h$','Interpreter','Latex')
        subplot(4,1,3); hold on
        ylabel('$\partial s/\partial w_\rho$','Interpreter','Latex')
        subplot(4,1,4); hold on
        ylabel('$\partial s/\partial v$','Interpreter','Latex')
        xlabel('$x$','Interpreter','Latex')
        switch type_of_adjoint
            case 'DA'
                subplot(4,1,1)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x ,ds_int.t /emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x1,ds_int.t/emode.M11)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.t/emode.M)
                end
                subplot(4,1,2)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x, ds_int.h /emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x1,ds_int.h/emode.M11)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.h/emode.M)
                end
                subplot(4,1,3)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.wr/emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.wr/emode.M00)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.wr/emode.M)
                end
                subplot(4,1,4)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.v/emode.M)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x0,ds_int.v/emode.M00)
                end
            case 'CA'
                subplot(4,1,1)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.t)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x1,ds_int.t)
                end
                subplot(4,1,2)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.h)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x1,ds_int.h)
                end
                subplot(4,1,3)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.wr)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.wr)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.wr)
                end
                subplot(4,1,4)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.v)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.v)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.v)
                end
        end
        % Display external sensitivities to screen
        disp('External sensitivities are:')
        disp(ds_ext)
    case 'fb_sens'
        % Plot the feedback sensitivities
        figure(fignum)
        clf
        subplot(3,2,1); hold on
        title('from $u$','Interpreter','Latex')
        ylabel('into mass eq.','Interpreter','Latex')
        subplot(3,2,2); hold on
        title('from $p$','Interpreter','Latex')
        subplot(3,2,3); hold on
        ylabel('into momentum eq.','Interpreter','Latex')
        subplot(3,2,4); hold on
        subplot(3,2,5); hold on
        ylabel('into energy eq.','Interpreter','Latex')
        xlabel('$x$','Interpreter','Latex')
        subplot(3,2,6); hold on
        xlabel('$x$','Interpreter','Latex')
        switch type_of_adjoint
            case 'DA'
                subplot(3,2,1)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.mru/emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.mru/emode.M00)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.mru/emode.M)
                end
                subplot(3,2,2)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.mrp/emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x1,ds_int.mrp/emode.M11)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.mrp/emode.M)
                end
                subplot(3,2,3)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.fru/emode.M)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x0,ds_int.fru/emode.M00)
                end
                subplot(3,2,4)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.frp/emode.M)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x0,ds_int.frp/emode.M00)
                end
                subplot(3,2,5)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.qpu/emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.qpu/emode.M00)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.qpu/emode.M)
                end
                subplot(3,2,6)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.qpp/emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x1,ds_int.qpp/emode.M11)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.qpp/emode.M)
                end
            case 'CA'
                subplot(3,2,1)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x ,ds_int.mru)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.mru)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.mru)
                end
                subplot(3,2,2)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x ,ds_int.mrp)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x1,ds_int.mrp)
                end
                subplot(3,2,3)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x ,ds_int.fru)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.fru)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.fru)
                end
                subplot(3,2,4)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x ,ds_int.frp)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.frp)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.frp)
                end
                subplot(3,2,5)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.qpu)
                    case 'FEW'
                        fun_plot_fun(emode.x0,ds_int.qpu)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.qpu)
                end
                subplot(3,2,6)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,ds_int.qpp)
                    case 'FEW'
                        fun_plot_fun(emode.x1,ds_int.qpp)
                    case 'SBP'
                        fun_plot_fun(emode.x1,ds_int.qpp)
                end
        end
    case 'rec'
        % Receptivities of the mass, momentum, and energy equations
        figure(fignum)
        clf
        subplot(3,1,1); hold on
        ylabel('$\dot{m}_\rho^\dag$','Interpreter','Latex')
        title([discretization,'\_',type_of_adjoint,'\_',lin_or_nonlin,'; eigenvalue s = ',num2str(emode.s)])
        subplot(3,1,2); hold on
        ylabel('$\dot{f}_\rho^\dag$','Interpreter','Latex')
        subplot(3,1,3); hold on
        ylabel('$\dot{q}_p^\dag$','Interpreter','Latex')
        xlabel('$x$','Interpreter','Latex')
        switch type_of_adjoint
            case 'DA'
                subplot(3,1,1)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,emode.mr/emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x1,emode.mr/emode.M11)
                    case 'SBP'
                        fun_plot_fun(emode.x1,emode.mr/emode.M)
                end
                subplot(3,1,2)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,emode.fr/emode.M)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x0,emode.fr/emode.M00)
                end
                subplot(3,1,3)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,emode.qp/emode.M)
                    case 'FEW'
                        fun_plot_fun(emode.x1,emode.qp/emode.M11)
                    case 'SBP'
                        fun_plot_fun(emode.x1,emode.qp/emode.M)
                end
            case 'CA'
                subplot(3,1,1)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,emode.mr)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x1,emode.mr)
                end
                subplot(3,1,2)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,emode.fr)
                    case 'FEW'
                        fun_plot_fun(emode.x0,emode.fr)
                    case 'SBP'
                        fun_plot_fun(emode.x1,emode.fr)
                end
                subplot(3,1,3)
                switch discretization
                    case {'FDS','FDW'}
                        fun_plot_fun(emode.x,emode.qp)
                    case {'FEW','SBP'}
                        fun_plot_fun(emode.x1,emode.qp)
                end
        end
    otherwise
        disp(['Plot type "',plot_type,'" not recognized'])
end

%% Convert to dimensional
s_dim = emode.s * ref.u_dim / ref.l_dim;
disp(['g.rate(Hz), freq(Hz) = ',num2str(s_dim/2/pi,16)])
disp(['g.rate(rad/sec), freq(rad/sec) = ',num2str(s_dim,16)])

end

function [] = fun_plot_fun(x,z)

plot(x,abs(z),'k',x,-abs(z),'k')
plot(x,real(z),'r')
plot(x,imag(z),'b')

end

