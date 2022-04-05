function [] = tabfun_baseflow(model,filename)
% tabfun_baseflow
%
% Create tables of the dimensional and nondimensional parameters for model

% Load in the dimensional parameters
[param_dim] = fun_param_dim(model);

% Convert to nondim parameters
[param,ref] = fun_nondim(param_dim);

%% Write the table
% Open the table file
eval(['fid = fopen(''figures/',filename,'_dim.tex'',''w'');']);
fprintf(fid,'$\\mean{\\rho}$      & (kgm$^{-3})$ & $ %5.2f $ \\\\ \n',param_dim.rhobar);
fprintf(fid,'$\\mean{p}$          & (Pa)         & $ %3.2e $ \\\\ \n',param_dim.pbar);
fprintf(fid,'$\\length_\\tube$    & (m)          & $ %5.2f $ \\\\ \n',param_dim.X);
fprintf(fid,'$\\area_\\tube$      & (m$^2$)      & $ %3.2e $ \\\\ \n',param_dim.S_c);
fprintf(fid,'$\\nold$             & (Jm$^{-1}$)  & $ %3.2e $ \\\\ \n',param_dim.eta);
fprintf(fid,'$\\tau$              & (s)          & $ %3.2e $ \\\\ \n',param_dim.tau);
fprintf(fid,'$\\loc_\\meas$       & (m)          & $ %6.3f $ \\\\ \n',param_dim.X_w);
fprintf(fid,'$\\length_\\meas$    & (m)          & $ %6.3f $ \\\\ \n',param_dim.L_w);
fprintf(fid,'$\\loc_\\hrel$       & (m)          & $ %6.3f $ \\\\ \n',param_dim.X_h);
fprintf(fid,'$\\length_\\hrel$    & (m)          & $ %6.3f $ \\\\ \n',param_dim.L_h);
fprintf(fid,'$\\overline{\\rho_u}$ & (kgm$^{-3})$ & $ %6.3f $ \\\\ \n ',param_dim.rhu);
fprintf(fid,'$\\overline{\\rho_d}$ & (kgm$^{-3})$ & $ %6.3f $ \\\\ \n ',param_dim.rhd);
fclose(fid);
% Reference values
eval(['fid = fopen(''figures/',filename,'_ref.tex'',''w'');']);
fprintf(fid,'$L_{ref}$       & (m)          & $ %5.2f $ \\\\ \n',ref.l_dim);
fprintf(fid,'$p_{ref}$       & (Pa)         & $ %3.2e $ \\\\ \n',ref.p_dim);
fprintf(fid,'$u_{ref}$       & (ms$^{-1}$)  & $ %5.0f $ \\\\ \n',ref.u_dim);
fclose(fid);
% nondim values
eval(['fid = fopen(''figures/',filename,'_nondim.tex'',''w'');'])
fprintf(fid,'$\\gamma$         &  --          & $ %5.2f $ \\\\ \n',param.gam);
fprintf(fid,'$R_u$             &  --          & $ %+6.3f %+6.3f \\ori $ \\\\ \n ',[real(param.Ru),imag(param.Ru)]);
fprintf(fid,'$R_d$             &  --          & $ %+6.3f %+6.3f \\ori $ \\\\ \n ',[real(param.Rd),imag(param.Rd)]);
fprintf(fid,'$n$               &  --          & $ %6.3f $ \\\\ \n',param.n);
fprintf(fid,'$\\tau$           &  --          & $ %6.3f $ \\\\ \n',param.tau);
fprintf(fid,'$\\loc_\\meas$    &  --          & $ %6.3f $ \\\\ \n',param.X_w);
fprintf(fid,'$\\length_\\meas$ &  --          & $ %6.3f $ \\\\ \n',param.L_w);
fprintf(fid,'$\\loc_\\hrel$    &  --          & $ %6.3f $ \\\\ \n',param.X_h);
fprintf(fid,'$\\length_\\hrel$ &  --          & $ %6.3f $ \\\\ \n',param.L_h);
fprintf(fid,'$\\rho_u$         &  --          & $ %6.3f $ \\\\ \n',param.rhu);
fprintf(fid,'$\\rho_d$         &  --          & $ %6.3f $ \\\\ \n',param.rhd);
fclose(fid);

% Return if you do not want to plot the distributions
% return

%% Read in the acoustic and heat release rate profiles
xb = 0:0.001:1;
% Density 
rh = fun_rh(param,xb);
% Heat release envelope divided by pressure, vp(x), which integrates to 1
h  = fun_h(param,xb);
% Measurement envelope divided by density, wr(x)
wr = fun_wr(param,xb);

%% Load in the colormap
[colmap] = figfun_colmap();

%% Plot the results on top of each other
figure(1); clf
FontSize = 18;
LineWidth = 3;
subplot(3,1,1); hold on
plot(xb,rh,'LineWidth',LineWidth,'Color',colmap(1,:))
ylabel('$\overline{\rho}$','Interpreter','Latex','FontSize',FontSize)
set(gca,'FontSize',FontSize,'FontName','Times')
box on; grid on
subplot(3,1,2); hold on
plot(xb,h,'LineWidth',LineWidth,'Color',colmap(1,:))
ylabel('$v$','Interpreter','Latex','FontSize',FontSize)
set(gca,'FontSize',FontSize,'FontName','Times')
box on; grid on
subplot(3,1,3); hold on
plot(xb,wr,'LineWidth',LineWidth,'Color',colmap(1,:))
ylabel('$w_{\overline{\rho}}$','Interpreter','Latex','FontSize',FontSize)
xlabel('$x$','Interpreter','Latex','FontSize',FontSize)
set(gca,'FontSize',FontSize,'FontName','Times')
box on; grid on

% print to file
eval(['print(''-depsc2'',''figures/',filename,'_fig.eps'')'])

