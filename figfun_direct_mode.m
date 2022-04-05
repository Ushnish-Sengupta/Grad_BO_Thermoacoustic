function [] = figfun_direct_mode(model,N,itmax,figname)
% figfun_direct_mode
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
scheme.itmax = itmax;

%% Calculate the eigenvalues with different methods
% Finite Difference strong form (linear method)
[emode_FDS_DA_linear] = fun_Helm('FDS','DA','linear',param,scheme);
% Finite Difference strong form (nonlinear method)
[emode_FDS_DA_nonlin] = fun_Helm('FDS','DA','nonlin',param,scheme);
% Finite Element (linear method)
[emode_FEW_DA_linear] = fun_Helm('FEW','DA','linear',param,scheme);
% Finite Element (nonlinear method)
[emode_FEW_DA_nonlin] = fun_Helm('FEW','DA','nonlin',param,scheme);
% Finite Difference weak form (linear method)
[emode_FDW_DA_linear] = fun_Helm('FDW','DA','linear',param,scheme);
% Finite Difference weak form (nonlinear method)
[emode_FDW_DA_nonlin] = fun_Helm('FDW','DA','nonlin',param,scheme);
% Finite Difference SBP-SAT (linear method)
[emode_SBP_DA_linear] = fun_Helm('SBP','DA','linear',param,scheme);
% Finite Difference SBP-SAT (nonlinear method)
[emode_SBP_DA_nonlin] = fun_Helm('SBP','DA','nonlin',param,scheme);

%% Plot the results on top of each other
figure(1); clf
subplot(4,1,1); hold on
plot(emode_FDS_DA_nonlin.x , real(emode_FDS_DA_nonlin.pr),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_DA_linear.x , real(emode_FDS_DA_linear.pr),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA_nonlin.x1, real(emode_FEW_DA_nonlin.pr),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_DA_linear.x1, real(emode_FEW_DA_linear.pr),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA_nonlin.x , real(emode_FDW_DA_nonlin.pr),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_DA_linear.x , real(emode_FDW_DA_linear.pr),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA_nonlin.x1, real(emode_SBP_DA_nonlin.pr),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_DA_linear.x1, real(emode_SBP_DA_linear.pr),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$p_r$',[],'type_1');
subplot(4,1,2); hold on
plot(emode_FDS_DA_nonlin.x , imag(emode_FDS_DA_nonlin.pr),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_DA_linear.x , imag(emode_FDS_DA_linear.pr),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA_nonlin.x1, imag(emode_FEW_DA_nonlin.pr),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_DA_linear.x1, imag(emode_FEW_DA_linear.pr),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA_nonlin.x , imag(emode_FDW_DA_nonlin.pr),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_DA_linear.x , imag(emode_FDW_DA_linear.pr),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA_nonlin.x1, imag(emode_SBP_DA_nonlin.pr),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_DA_linear.x1, imag(emode_SBP_DA_linear.pr),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$p_i$',[],'type_1');
subplot(4,1,3); hold on
plot(emode_FDS_DA_nonlin.x , real(emode_FDS_DA_nonlin.ur),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_DA_linear.x , real(emode_FDS_DA_linear.ur),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA_nonlin.x0, real(emode_FEW_DA_nonlin.ur),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_DA_linear.x0, real(emode_FEW_DA_linear.ur),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA_nonlin.x , real(emode_FDW_DA_nonlin.ur),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_DA_linear.x , real(emode_FDW_DA_linear.ur),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA_nonlin.x1, real(emode_SBP_DA_nonlin.ur),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_DA_linear.x1, real(emode_SBP_DA_linear.ur),'Color',colmap(8,:),'LineWidth',1)
figfun_format([],'$u_r$',[],'type_1');
subplot(4,1,4); hold on
plot(emode_FDS_DA_nonlin.x , imag(emode_FDS_DA_nonlin.ur),'Color',colmap(1,:),'LineWidth',8)
plot(emode_FDS_DA_linear.x , imag(emode_FDS_DA_linear.ur),'Color',colmap(2,:),'LineWidth',7)
plot(emode_FEW_DA_nonlin.x0, imag(emode_FEW_DA_nonlin.ur),'Color',colmap(3,:),'LineWidth',6)
plot(emode_FEW_DA_linear.x0, imag(emode_FEW_DA_linear.ur),'Color',colmap(4,:),'LineWidth',5)
plot(emode_FDW_DA_nonlin.x , imag(emode_FDW_DA_nonlin.ur),'Color',colmap(5,:),'LineWidth',4)
plot(emode_FDW_DA_linear.x , imag(emode_FDW_DA_linear.ur),'Color',colmap(6,:),'LineWidth',3)
plot(emode_SBP_DA_nonlin.x1, imag(emode_SBP_DA_nonlin.ur),'Color',colmap(7,:),'LineWidth',2)
plot(emode_SBP_DA_linear.x1, imag(emode_SBP_DA_linear.ur),'Color',colmap(8,:),'LineWidth',1)
figfun_format('$x$','$u_i$',[],'type_1');

%% Print the figure to file
if figname
    % Change the paper position
    set(gcf,'PaperPosition',[0.6350    6.3500   20.3200   20.3200])
    eval(['print(''-depsc2'',''figures/',figname,'.eps'')'])
    % Create the text to go under the figure
    eval(['fid = fopen(''figures/',figname,'.tex'',''w'');'])
    fprintf(fid,' & $N = %u$ \\\\ \n',scheme.N);
    fprintf(fid,'\\texttt{FDS\\_DA\\_nonlin} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FDS_DA_nonlin.s), imag(emode_FDS_DA_nonlin.s)]);
    fprintf(fid,'\\texttt{FDS\\_DA\\_linear} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FDS_DA_linear.s), imag(emode_FDS_DA_linear.s)]);
    fprintf(fid,'\\texttt{FEW\\_DA\\_nonlin} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FEW_DA_nonlin.s), imag(emode_FEW_DA_nonlin.s)]);
    fprintf(fid,'\\texttt{FEW\\_DA\\_linear} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FEW_DA_linear.s), imag(emode_FEW_DA_linear.s)]);
    fprintf(fid,'\\texttt{FDW\\_DA\\_nonlin} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FDW_DA_nonlin.s), imag(emode_FDW_DA_nonlin.s)]);
    fprintf(fid,'\\texttt{FDW\\_DA\\_linear} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_FDW_DA_linear.s), imag(emode_FDW_DA_linear.s)]);
    fprintf(fid,'\\texttt{SBP\\_DA\\_nonlin} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_SBP_DA_nonlin.s), imag(emode_SBP_DA_nonlin.s)]);
    fprintf(fid,'\\texttt{SBP\\_DA\\_linear} & $s = %10.8f %+10.8f \\ori$ \\\\ \n',[real(emode_SBP_DA_linear.s), imag(emode_SBP_DA_linear.s)]);
    fclose(fid);
    % Create the legend
    labels = { ...
        'FDS\_DA\_nonlin'
        'FDS\_DA\_linear'
        'FEW\_DA\_nonlin'
        'FEW\_DA\_linear'
        'FDW\_DA\_nonlin'
        'FDW\_DA\_linear'
        'SBP\_DA\_nonlin'
        'SBP\_DA\_linear'
        };
    figfun_legend(colmap,labels,figname);
end


end
