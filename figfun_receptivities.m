function [] = figfun_receptivities(model,N,lin_or_nonlin,figname)
% figfun_receptivities
%
% INPUTS
% model         model ('Rijke' or 'GT')
% N             N+1 = number of gridpoints (FD or WD) or number of elements (FE)
% lin_or_nonlin iteration method
% figname       [] to plot to screen or the filename to print to file

%% Set the plot parameters
colmap = figfun_colmap;

%% Set the dimensional parameters
param_dim = fun_param_dim(model);

%% Calculate the reference scales and the nondimensional parameters
param = fun_nondim(param_dim);

%% Set the starting value of s, the numerical scheme and max number of iterations
scheme.s0    = fun_set_s0(param);
scheme.N     = N;
scheme.itmax = 10;

%% Calculate the eigenvalues with different methods
% Finite Difference strong form (nonlinear method) DA
[emode_FDS_DA] = fun_Helm('FDS','DA',lin_or_nonlin,param,scheme);
% Finite Element (nonlinear method) DA
[emode_FEW_DA] = fun_Helm('FEW','DA',lin_or_nonlin,param,scheme);
% Finite Difference weak form (nonlinear method) DA
[emode_FDW_DA] = fun_Helm('FDW','DA',lin_or_nonlin,param,scheme);
% Finite Difference strong form with SBP-SAT (nonlinear method) DA
[emode_SBP_DA] = fun_Helm('SBP','DA',lin_or_nonlin,param,scheme);
% Take the c.c. of scheme.s0 as starting point for the CA algorithms
scheme.s0 = conj(scheme.s0);
% Calculate the eigenvalues and CA sensitivities with different methods
% Finite Element (nonlinear method) DA
[emode_FEW_CA] = fun_Helm('FEW','CA',lin_or_nonlin,param,scheme);
% Finite Difference strong form (nonlinear method) DA
[emode_FDS_CA] = fun_Helm('FDS','CA',lin_or_nonlin,param,scheme);
% Finite Difference weak form (nonlinear method) DA
[emode_FDW_CA] = fun_Helm('FDW','CA',lin_or_nonlin,param,scheme);
% Finite Difference strong form with SBP-SAT (nonlinear method) CA
[emode_SBP_CA] = fun_Helm('SBP','CA',lin_or_nonlin,param,scheme);

%% Plot the results on top of each other
figure(1); clf
subplot(3,1,1); hold on
plot(emode_FDS_DA.x , abs(emode_FDS_DA.mr/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(emode_FDS_CA.mr                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, abs(emode_FEW_DA.mr/emode_FEW_DA.M11),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, abs(emode_FEW_CA.mr                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(emode_FDW_DA.mr/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(emode_FDW_CA.mr                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(emode_SBP_DA.mr/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(emode_SBP_CA.mr                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$|\dot{m}_{\overline{\rho}}^\dag|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +5.00])
    case 'GT'
        set(gca,'YLim',[ 0.00 +3.00])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +15.0])
end
subplot(3,1,2); hold on
plot(emode_FDS_DA.x , abs(emode_FDS_DA.fr/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(emode_FDS_CA.fr                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x0, abs(emode_FEW_DA.fr/emode_FEW_DA.M00),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x0, abs(emode_FEW_CA.fr                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(emode_FDW_DA.fr/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(emode_FDW_CA.fr                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x0, abs(emode_SBP_DA.fr/emode_SBP_DA.M00),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(emode_SBP_CA.fr                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$|\dot{f}_{\overline{\rho}}^\dag|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +7.00])
    case 'GT'
        set(gca,'YLim',[ 0.00 +10.0])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +150.])
end
subplot(3,1,3); hold on
plot(emode_FDS_DA.x , abs(emode_FDS_DA.qp/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(emode_FDS_CA.qp                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, abs(emode_FEW_DA.qp/emode_FEW_DA.M11),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, abs(emode_FEW_CA.qp                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(emode_FDW_DA.qp/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(emode_FDW_CA.qp                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(emode_SBP_DA.qp/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(emode_SBP_CA.qp                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format('$x$','$|\dot{q}_{\overline{p}}^\dag|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +1.50])
    case 'GT'
        set(gca,'YLim',[ 0.00 +1.00])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +4.00])
end

%% Print the figure to file
if figname
    % Change the paper position
    set(gcf,'PaperPosition',[0.6350    6.3500   20.3200   20.3200])
    eval(['print(''-depsc2'',''figures/',figname,'.eps'')'])
    % Create the legend
    labels = { ...
        'FDS\_DA'
        'FDS\_CA'
        'FEW\_DA'
        'FEW\_CA'
        'FDW\_DA'
        'FDW\_CA'
        'SBP\_DA'
        'SBP\_CA'
        };
    figfun_legend(colmap,labels,figname);
end

end
