function [] = figfun_base_sen_z(model,N,lin_or_nonlin,figname)
% figfun_base_sen_z
%
% INPUTS
% model     model ('Rijke' or 'GT')
% N         N+1 = number of gridpoints (FD or WD) or number of elements (FE)
% itmax     maximum number of iterations for active iteration method
% figname   '' to plot to screen or the figure name to print to file

%% Set the plot parameters
colmap = figfun_colmap;
col_imag       = colmap(1,:);
col_real       = colmap(8,:);
col_imag_SA    = colmap(3,:);
col_real_SA    = colmap(7,:);

%% Set the dimensional parameters
param_dim = fun_param_dim(model);

%% Calculate the reference scales and the nondimensional parameters
param = fun_nondim(param_dim);

%% Set the starting value of s, the numerical scheme and max number of iterations
scheme.s0    = fun_set_s0(param);
scheme.N     = N;
scheme.itmax = 10;

%% Calculate the eigenvalues and DA sensitivities with FE DA 
% Finite Element (nonlinear method) DA
[emode_FEW_DA,ds_FEW_DA] = fun_Helm('FEW','DA',lin_or_nonlin,param,scheme);

%% Calculate the sensitivities using the direct pressure instead of the adj
% Extract the eigenvalue and right eigenvector
s  = emode_FEW_DA.s;
pr = emode_FEW_DA.pr;
% Load in the FE matrices
[~,mat] = mat_FE(param,N);
% Calculate dAds and dA at the final value of s and initialize ds_int
[~,C,dAds,~,dA,ds_int] = mat_AC_FEW_DA(mat,param,N,s);
% Calculate the denominator that is common to all sensitivities
dGds = dAds - 2*s*C; ip = -(pr' * dGds * pr);
% Calculate internal base state sensitivities and feedback sensitivities
ds_FEW_DA_SA = fun_ds_DA(ds_int,1,pr,dA,pr,ip,N);

%% Plot the results on top of each other
figure(1); clf
subplot(4,1,1); hold on
plot(emode_FEW_DA.x1, imag(ds_FEW_DA_SA.t/emode_FEW_DA.M11),'Color',col_imag_SA,'LineWidth',4)
plot(emode_FEW_DA.x1, real(ds_FEW_DA_SA.t/emode_FEW_DA.M11),'Color',col_real_SA,'LineWidth',4)
plot(emode_FEW_DA.x1, imag(ds_FEW_DA.t   /emode_FEW_DA.M11),'Color',col_imag   ,'LineWidth',2)
plot(emode_FEW_DA.x1, real(ds_FEW_DA.t   /emode_FEW_DA.M11),'Color',col_real   ,'LineWidth',2)
figfun_format([],'$|\partial s / \partial \tau|$',[],'type_1')
switch model
    case 'Rijke'
        set(gca,'YLim',[-3.00 +3.00])
    case 'GT'
        set(gca,'YLim',[-0.20 +0.20])
    case 'rocket'
        set(gca,'YLim',[-30.0 +30.0])
end
subplot(4,1,2); hold on
plot(emode_FEW_DA.x1, imag(ds_FEW_DA_SA.h/emode_FEW_DA.M11),'Color',col_imag_SA,'LineWidth',4)
plot(emode_FEW_DA.x1, real(ds_FEW_DA_SA.h/emode_FEW_DA.M11),'Color',col_real_SA,'LineWidth',4)
plot(emode_FEW_DA.x1, imag(ds_FEW_DA.h   /emode_FEW_DA.M11),'Color',col_imag   ,'LineWidth',2)
plot(emode_FEW_DA.x1, real(ds_FEW_DA.h   /emode_FEW_DA.M11),'Color',col_real   ,'LineWidth',2)
figfun_format([],'$|\partial s / \partial h|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[-0.04 +0.04])
    case 'GT'
        set(gca,'YLim',[-0.05 +0.05])
    case 'rocket'
        set(gca,'YLim',[-10.0 +10.0])
end
subplot(4,1,3); hold on
plot(emode_FEW_DA.x0, imag(ds_FEW_DA_SA.wr/emode_FEW_DA.M00),'Color',col_imag_SA,'LineWidth',4)
plot(emode_FEW_DA.x0, real(ds_FEW_DA_SA.wr/emode_FEW_DA.M00),'Color',col_real_SA,'LineWidth',4)
plot(emode_FEW_DA.x0, imag(ds_FEW_DA.wr   /emode_FEW_DA.M00),'Color',col_imag   ,'LineWidth',2)
plot(emode_FEW_DA.x0, real(ds_FEW_DA.wr   /emode_FEW_DA.M00),'Color',col_real   ,'LineWidth',2)
figfun_format([],'$|\partial s / \partial w_{\overline{\rho}}|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[-0.06 +0.06])
    case 'GT'
        set(gca,'YLim',[-0.05 +0.05])
    case 'rocket'
        set(gca,'YLim',[-10.0 +10.0])
end
subplot(4,1,4); hold on
plot(emode_FEW_DA.x0, imag(ds_FEW_DA_SA.v/emode_FEW_DA.M00),'Color',col_imag_SA,'LineWidth',4)
plot(emode_FEW_DA.x0, real(ds_FEW_DA_SA.v/emode_FEW_DA.M00),'Color',col_real_SA,'LineWidth',4)
plot(emode_FEW_DA.x0, imag(ds_FEW_DA.v   /emode_FEW_DA.M00),'Color',col_imag   ,'LineWidth',2)
plot(emode_FEW_DA.x0, real(ds_FEW_DA.v   /emode_FEW_DA.M00),'Color',col_real   ,'LineWidth',2)
figfun_format('$x$','$|\partial s / \partial \overline{v}|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[-6.00 +06.0])
    case 'GT'
        set(gca,'YLim',[-6.00 +06.0])
    case 'rocket'
        set(gca,'YLim',[-10.0 +10.0])
end

%% Print the figure to file
if figname
    % Change the paper position
    set(gcf,'PaperPosition',[0.6350    6.3500   20.3200   20.3200])
    eval(['print(''-depsc2'',''figures/',figname,'.eps'')'])
    % Create the legend
    figure(2); clf
    subplot(7,1,1); hold on
    axis off
    plot([0 0.8],[3 3],'-','Color',col_real,'LineWidth',2)
    text(0.9,3,'real component (influence on growth rate)','Interpreter','Latex')
    plot([0 0.8],[2 2],'-','Color',col_imag,'LineWidth',2)
    text(0.9,2,'imaginary componenet (influence on frequency)','Interpreter','Latex')
    plot([0 0.8],[1 1],'-','Color',col_real_SA,'LineWidth',4)
    text(0.9,1,'real component formed with $p$ instead of $p^\dag$','Interpreter','Latex')
    plot([0 0.8],[0 0],'-','Color',col_imag_SA,'LineWidth',4)
    text(0.9,0,'imaginary componenet formed with $p$ instead of $p^\dag$','Interpreter','Latex')
    axis([0 4 0 3])
    eval(['print(''-depsc2'',''figures/',figname,'_legend.eps'')'])
end

end
