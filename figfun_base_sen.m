function [] = figfun_base_sen(model,N,lin_or_nonlin,figname)
% figfun_base_sen
%
% INPUTS
% model     model ('Rijke' or 'GT')
% N         N+1 = number of gridpoints (FD or WD) or number of elements (FE)
% itmax     maximum number of iterations for active iteration method
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

%% Calculate the eigenvalues and DA sensitivities with different methods
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
% Calculate the eigenvalues and CA sensitivities with different methods
% Finite Element (nonlinear method) DA
[emode_FEW_CA,ds_FEW_CA] = fun_Helm('FEW','CA',lin_or_nonlin,param,scheme);
% Finite Difference strong form (nonlinear method) DA
[emode_FDS_CA,ds_FDS_CA] = fun_Helm('FDS','CA',lin_or_nonlin,param,scheme);
% Finite Difference weak form (nonlinear method) DA
[emode_FDW_CA,ds_FDW_CA] = fun_Helm('FDW','CA',lin_or_nonlin,param,scheme);
% Finite Difference strong form with SBP-SAT (nonlinear method) DA
[emode_SBP_CA,ds_SBP_CA] = fun_Helm('SBP','CA',lin_or_nonlin,param,scheme);

%% Plot the results on top of each other
figure(1); clf
subplot(4,1,1); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.t/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.t                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, abs(ds_FEW_DA.t/emode_FEW_DA.M11),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, abs(ds_FEW_CA.t                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.t/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.t                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(ds_SBP_DA.t/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.t                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$|\partial s / \partial \tau|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +3.00])
    case 'GT'
        set(gca,'YLim',[ 0.00 +0.20])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +40.0])
end
subplot(4,1,2); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.h/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.h                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, abs(ds_FEW_DA.h/emode_FEW_DA.M11),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, abs(ds_FEW_CA.h                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.h/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.h                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(ds_SBP_DA.h/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.h                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$|\partial s / \partial h|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[-0.00 +0.04])
    case 'GT'
        set(gca,'YLim',[-0.00 +0.05])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +10.0])
end
subplot(4,1,3); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.wr/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.wr                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x0, abs(ds_FEW_DA.wr/emode_FEW_DA.M00),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x0, abs(ds_FEW_CA.wr                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.wr/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.wr                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, abs(ds_SBP_DA.wr/emode_SBP_DA.M  ),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.wr                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$|\partial s / \partial w_{\overline{\rho}}|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +0.06])
    case 'GT'
        set(gca,'YLim',[+0.00 +0.05])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +10.0])
end
subplot(4,1,4); hold on
plot(emode_FDS_DA.x , abs(ds_FDS_DA.v/emode_FDS_DA.M  ),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , abs(ds_FDS_CA.v                 ),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x0, abs(ds_FEW_DA.v/emode_FEW_DA.M00),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x0, abs(ds_FEW_CA.v                 ),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , abs(ds_FDW_DA.v/emode_FDW_DA.M  ),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , abs(ds_FDW_CA.v                 ),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x0, abs(ds_SBP_DA.v/emode_SBP_DA.M00),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, abs(ds_SBP_CA.v                 ),'Color',colmap(8,:),'LineWidth',1)
figfun_format('$x$','$|\partial s / \partial \overline{v}|$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +06.0])
    case 'GT'
        set(gca,'YLim',[ 0.00 +06.0])
    case 'rocket'
        set(gca,'YLim',[ 0.00 +50.0])
end

%% Print the figure to file
if figname
    % Change the paper position
    set(gcf,'PaperPosition',[0.6350    6.3500   20.3200   20.3200])
    eval(['print(''-depsc2'',''figures/',figname,'.eps'')'])
    % Create the text to go under the figure
    eval(['fid = fopen(''figures/',figname,'.tex'',''w'');'])
    fprintf(fid,' & $N = %u$ \\\\ \n',scheme.N);
    fprintf(fid,'\\texttt{FDS\\_DA} & $\\partial s / \\partial n = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDS_DA.n), imag(ds_FDS_DA.n)]);
    fprintf(fid,'\\texttt{FDS\\_CA} & $\\partial s / \\partial n = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDS_CA.n), imag(ds_FDS_CA.n)]);
    fprintf(fid,'\\texttt{FEW\\_DA} & $\\partial s / \\partial n = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FEW_DA.n), imag(ds_FEW_DA.n)]);
    fprintf(fid,'\\texttt{FEW\\_CA} & $\\partial s / \\partial n = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FEW_CA.n), imag(ds_FEW_CA.n)]);
    fprintf(fid,'\\texttt{FDW\\_DA} & $\\partial s / \\partial n = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDW_DA.n), imag(ds_FDW_DA.n)]);
    fprintf(fid,'\\texttt{FDW\\_CA} & $\\partial s / \\partial n = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDW_CA.n), imag(ds_FDW_CA.n)]);
    fprintf(fid,'\\texttt{SBP\\_DA} & $\\partial s / \\partial n = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_SBP_DA.n), imag(ds_SBP_DA.n)]);
    fprintf(fid,'\\texttt{SBP\\_CA} & $\\partial s / \\partial n = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_SBP_CA.n), imag(ds_SBP_CA.n)]);
    fprintf(fid,'\\hline');
    fprintf(fid,'\\texttt{FDS\\_DA} & $\\partial s / \\partial k_u = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDS_DA.ku), imag(ds_FDS_DA.ku)]);
    fprintf(fid,'\\texttt{FDS\\_CA} & $\\partial s / \\partial k_u = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDS_CA.ku), imag(ds_FDS_CA.ku)]);
    fprintf(fid,'\\texttt{FEW\\_DA} & $\\partial s / \\partial k_u = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FEW_DA.ku), imag(ds_FEW_DA.ku)]);
    fprintf(fid,'\\texttt{FEW\\_CA} & $\\partial s / \\partial k_u = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FEW_CA.ku), imag(ds_FEW_CA.ku)]);
    fprintf(fid,'\\texttt{FDW\\_DA} & $\\partial s / \\partial k_u = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDW_DA.ku), imag(ds_FDW_DA.ku)]);
    fprintf(fid,'\\texttt{FDW\\_CA} & $\\partial s / \\partial k_u = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDW_CA.ku), imag(ds_FDW_CA.ku)]);
    fprintf(fid,'\\texttt{SBP\\_DA} & $\\partial s / \\partial k_u = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_SBP_DA.ku), imag(ds_SBP_DA.ku)]);
    fprintf(fid,'\\texttt{SBP\\_CA} & $\\partial s / \\partial k_u = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_SBP_CA.ku), imag(ds_SBP_CA.ku)]);
    fprintf(fid,'\\hline');
    fprintf(fid,'\\texttt{FDS\\_DA} & $\\partial s / \\partial k_d = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDS_DA.kd), imag(ds_FDS_DA.kd)]);
    fprintf(fid,'\\texttt{FDS\\_CA} & $\\partial s / \\partial k_d = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDS_CA.kd), imag(ds_FDS_CA.kd)]);
    fprintf(fid,'\\texttt{FEW\\_DA} & $\\partial s / \\partial k_d = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FEW_DA.kd), imag(ds_FEW_DA.kd)]);
    fprintf(fid,'\\texttt{FEW\\_CA} & $\\partial s / \\partial k_d = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FEW_CA.kd), imag(ds_FEW_CA.kd)]);
    fprintf(fid,'\\texttt{FDW\\_DA} & $\\partial s / \\partial k_d = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDW_DA.kd), imag(ds_FDW_DA.kd)]);
    fprintf(fid,'\\texttt{FDW\\_CA} & $\\partial s / \\partial k_d = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_FDW_CA.kd), imag(ds_FDW_CA.kd)]);
    fprintf(fid,'\\texttt{SBP\\_DA} & $\\partial s / \\partial k_d = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_SBP_DA.kd), imag(ds_SBP_DA.kd)]);
    fprintf(fid,'\\texttt{SBP\\_CA} & $\\partial s / \\partial k_d = %+10.8f %+10.8f \\ori$ \\\\ \n',[real(ds_SBP_CA.kd), imag(ds_SBP_CA.kd)]);
    fclose(fid);
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
