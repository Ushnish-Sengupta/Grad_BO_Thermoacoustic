function [] = figfun_adjoint_mode(model,N,lin_or_nonlin,figname)
% figfun_adjoint_mode
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
% Finite Difference strong form (nonlinear method) CA
[emode_FDS_CA] = fun_Helm('FDS','CA',lin_or_nonlin,param,scheme);
% Finite Element (nonlinear method) CA
[emode_FEW_CA] = fun_Helm('FEW','CA',lin_or_nonlin,param,scheme);
% Finite Difference weak form (nonlinear method) CA
[emode_FDW_CA] = fun_Helm('FDW','CA',lin_or_nonlin,param,scheme);
% Finite Difference strong form with SBP-SAT (nonlinear method) CA
[emode_SBP_CA] = fun_Helm('SBP','CA',lin_or_nonlin,param,scheme);

%% Plot the results on top of each other
figure(1); clf
subplot(4,1,1); hold on
plot(emode_FDS_DA.x , real(emode_FDS_DA.pl),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , real(emode_FDS_CA.pr),'Color',colmap(2,:),'LineWidth',7) %
plot(emode_FEW_DA.x1, real(emode_FEW_DA.pl),'Color',colmap(3,:),'LineWidth',6) %
plot(emode_FEW_CA.x1, real(emode_FEW_CA.pr),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , real(emode_FDW_DA.pl),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , real(emode_FDW_CA.pr),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, real(emode_SBP_DA.pl),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, real(emode_SBP_CA.pr),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$p^\dag_r$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[ 0.00 +2.00])
    case 'GT'
        set(gca,'YLim',[ 0.00 +2.00])
    case 'rocket'
        set(gca,'YLim',[-3.00 +3.00])
end
subplot(4,1,2); hold on
plot(emode_FDS_DA.x , imag(emode_FDS_DA.pl),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , imag(emode_FDS_CA.pr),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, imag(emode_FEW_DA.pl),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, imag(emode_FEW_CA.pr),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , imag(emode_FDW_DA.pl),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , imag(emode_FDW_CA.pr),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, imag(emode_SBP_DA.pl),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, imag(emode_SBP_CA.pr),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$p^\dag_i$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'YLim',[-0.05 +0.05])
    case 'GT'
        set(gca,'YLim',[-0.10 +0.10])
    case 'rocket'
        set(gca,'YLim',[-3.00 +3.00])
end
subplot(4,1,3); hold on
plot(emode_FDS_DA.x , real(emode_FDS_DA.pl),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , real(emode_FDS_CA.pr),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, real(emode_FEW_DA.pl),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, real(emode_FEW_CA.pr),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , real(emode_FDW_DA.pl),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , real(emode_FDW_CA.pr),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, real(emode_SBP_DA.pl),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, real(emode_SBP_CA.pr),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$p^\dag_r$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'XLim',[ 0.00 +0.01])
        set(gca,'YLim',[ 0.00 +0.10])
    case 'GT'
        set(gca,'XLim',[+0.99 +1.00])
        set(gca,'YLim',[+1.00 +1.50])
    case 'rocket'
        set(gca,'XLim',[ 0.00 +0.01])
        set(gca,'YLim',[-3.50 -2.50])
end
subplot(4,1,4); hold on
plot(emode_FDS_DA.x , imag(emode_FDS_DA.pl),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_CA.x , imag(emode_FDS_CA.pr),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA.x1, imag(emode_FEW_DA.pl),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_CA.x1, imag(emode_FEW_CA.pr),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA.x , imag(emode_FDW_DA.pl),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_CA.x , imag(emode_FDW_CA.pr),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA.x1, imag(emode_SBP_DA.pl),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_CA.x1, imag(emode_SBP_CA.pr),'Color',colmap(8,:),'LineWidth',1)
figfun_format('$x$','$p^\dag_i$',[],'type_1');
switch model
    case 'Rijke'
        set(gca,'XLim',[ 0.00 +0.01])
        set(gca,'YLim',[ 0.00 +0.04])
    case 'GT'
        set(gca,'XLim',[+0.99 +1.00])
        set(gca,'YLim',[-0.00 +0.05])
    case 'rocket'
        set(gca,'XLim',[ 0.00 +0.01])
        set(gca,'YLim',[-2.50 -1.50])
end

%% Print the figure to file
if figname
    % Change the paper position
    set(gcf,'PaperPosition',[0.6350    6.3500   20.3200   20.3200])
    eval(['print(''-depsc2'',''figures/',figname,'.eps'')'])
    % Create the text to go under the figure
    eval(['fid = fopen(''figures/',figname,'.tex'',''w'');'])
    fprintf(fid,' & $N = %u$ \\\\ \n',scheme.N);
    fprintf(fid,'\\texttt{FDS\\_DA} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FDS_DA.s), imag(emode_FDS_DA.s)]);
    fprintf(fid,'\\texttt{FDS\\_CA} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FDS_CA.s), imag(emode_FDS_CA.s)]);
    fprintf(fid,'\\texttt{FEW\\_DA} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FEW_DA.s), imag(emode_FEW_DA.s)]);
    fprintf(fid,'\\texttt{FEW\\_CA} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FEW_CA.s), imag(emode_FEW_CA.s)]);
    fprintf(fid,'\\texttt{FDW\\_DA} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FDW_DA.s), imag(emode_FDW_DA.s)]);
    fprintf(fid,'\\texttt{FDW\\_CA} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FDW_CA.s), imag(emode_FDW_CA.s)]);
    fprintf(fid,'\\texttt{SBP\\_DA} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_SBP_DA.s), imag(emode_SBP_DA.s)]);
    fprintf(fid,'\\texttt{SBP\\_CA} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_SBP_CA.s), imag(emode_SBP_CA.s)]);
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
