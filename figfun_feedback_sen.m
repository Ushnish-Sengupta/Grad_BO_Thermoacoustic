function [] = figfun_feedback_sen(model,N,lin_or_nonlin,figname)
% figfun_feedback_sen
%
% INPUTS
% model     model ('Rijke' or 'GT')
% N         N+1 = number of gridpoints (FD or WD) or number of elements (FE)
% figname   '' to plot to screen or the figure name to print to file

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
[emode_FDS_DA,ds_FDS_DA] = fun_Helm('FDS','DA',lin_or_nonlin,param,scheme);
% Finite Element (nonlinear method) DA
[emode_FEW_DA,ds_FEW_DA] = fun_Helm('FEW','DA',lin_or_nonlin,param,scheme);
% Finite Difference weak form (nonlinear method) DA
[emode_FDW_DA,ds_FDW_DA] = fun_Helm('FDW','DA',lin_or_nonlin,param,scheme);
% Finite Difference strong form with SBP-SAT (nonlinear method) DA
[emode_SBP_DA,ds_SBP_DA] = fun_Helm('SBP','DA',lin_or_nonlin,param,scheme);
% Take the c.c. of scheme.s0 as starting point for the CA algorithms
scheme.s0 = conj(scheme.s0);
% Finite Difference strong form (nonlinear method) CA
[emode_FDS_CA,ds_FDS_CA] = fun_Helm('FDS','CA',lin_or_nonlin,param,scheme);
% Finite Element (nonlinear method) CA
[emode_FEW_CA,ds_FEW_CA] = fun_Helm('FEW','CA',lin_or_nonlin,param,scheme);
% Finite Difference weak form (nonlinear method) CA
[emode_FDW_CA,ds_FDW_CA] = fun_Helm('FDW','CA',lin_or_nonlin,param,scheme);
% Finite Difference strong form with SBP-SAT (nonlinear method) CA
[emode_SBP_CA,ds_SBP_CA] = fun_Helm('SBP','CA',lin_or_nonlin,param,scheme);

%% Plot the results on top of each other
figure(1); clf
FontSize = 18;
subplot(3,2,1); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.mru/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.mru                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x0, abs(ds_FEW_DA.mru/emode_FEW_DA.M00),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x0, abs(ds_FEW_CA.mru                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.mru/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.mru                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(ds_SBP_DA.mru/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.mru                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'to mass eq.','from $u$','type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +0.80])
    case 'GT'
        set(gca,'YLim',[ 0.00 +0.60])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +1.40])
end
yc = ylim; text(0.05,0.9*yc(2),'$|\partial s/ \partial \dot{m}_{\overline{\rho},u}|$','Interpreter','Latex','FontSize',FontSize)
subplot(3,2,2); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.mrp/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.mrp                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, abs(ds_FEW_DA.mrp/emode_FEW_DA.M11),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, abs(ds_FEW_CA.mrp                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.mrp/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.mrp                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(ds_SBP_DA.mrp/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.mrp                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],[],'from $p$','type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +2.00])
    case 'GT'
        set(gca,'YLim',[ 0.00 +1.40])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +6.00])
end
yc = ylim; text(0.05,0.9*yc(2),'$|\partial s/ \partial \dot{m}_{\overline{\rho},p}|$','Interpreter','Latex','FontSize',FontSize)
subplot(3,2,3); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.fru/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.fru                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x0, abs(ds_FEW_DA.fru/emode_FEW_DA.M00),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x0, abs(ds_FEW_CA.fru                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.fru/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.fru                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x0, abs(ds_SBP_DA.fru/emode_SBP_DA.M00),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.fru                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'to mv eq.',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +1.50])
    case 'GT'
        set(gca,'YLim',[ 0.00 +2.00])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +10.0])
end
yc = ylim; text(0.05,0.9*yc(2),'$|\partial s/ \partial {f}_{\overline{\rho},u}|$','Interpreter','Latex','FontSize',FontSize)
subplot(3,2,4); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.frp/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.frp                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x0, abs(ds_FEW_DA.frp/emode_FEW_DA.M00),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x0, abs(ds_FEW_CA.frp                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.frp/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.frp                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x0, abs(ds_SBP_DA.frp/emode_SBP_DA.M00),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.frp                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],[],[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +1.00])
    case 'GT'
        set(gca,'YLim',[ 0.00 +1.00])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +50.0])
end
yc = ylim; text(0.05,0.9*yc(2),'$|\partial s/ \partial {f}_{\overline{\rho},p}|$','Interpreter','Latex','FontSize',FontSize)
subplot(3,2,5); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.qpu/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.qpu                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x0, abs(ds_FEW_DA.qpu/emode_FEW_DA.M00),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x0, abs(ds_FEW_CA.qpu                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.qpu/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.qpu                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(ds_SBP_DA.qpu/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.qpu                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format('$x$','to energy eq.',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +0.20])
    case 'GT'
        set(gca,'YLim',[ 0.00 +0.20])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +0.40])
end
yc = ylim; text(0.05,0.9*yc(2),'$|\partial s/ \partial \dot{q}_{\overline{p},u}|$','Interpreter','Latex','FontSize',FontSize)
subplot(3,2,6); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.qpp/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.qpp                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, abs(ds_FEW_DA.qpp/emode_FEW_DA.M11),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, abs(ds_FEW_CA.qpp                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.qpp/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.qpp                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(ds_SBP_DA.qpp/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.qpp                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format('$x$',[],[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +0.50])
    case 'GT'
        set(gca,'YLim',[ 0.00 +0.40])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +2.00])
end
yc = ylim; text(0.05,0.9*yc(2),'$|\partial s/ \partial \dot{q}_{\overline{p},p}|$','Interpreter','Latex','FontSize',FontSize)

%% Print the figure to file
if figname
    % Change the paper position
    set(gcf,'PaperPosition',[0.6350    6.3500   20.3200   20.3200])
    eval(['print(''-depsc2'',''figures/',figname,'.eps'')'])
    % Create the legend
    labels = { ...
        'FDS\_DA\_nonlin'
        'FDS\_CA\_nonlin'
        'FEW\_DA\_nonlin'
        'FEW\_CA\_nonlin'
        'FDW\_DA\_nonlin'
        'FDW\_CA\_nonlin'
        'SBP\_DA\_nonlin'
        'SBP\_CA\_nonlin'
        };
    figfun_legend(colmap,labels,figname);
end

end
