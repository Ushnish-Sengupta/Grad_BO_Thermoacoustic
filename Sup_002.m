function [] = Sup_002
% Sup_002
%
% Plot the base state sensitivities as a function of the arbitrary
% parameter param.cu, which is used in the imposition of the upstream
% boundary condition in the FD case. 
%
% The base state sensitivities oscillate but these oscillations do not
% depend on the arbitrary parameter param.cu

colmap = figfun_colmap;

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');

%% Calculate the reference scales and the nondimensional parameters
param = fun_nondim(param_dim);

%% Set the starting value of s, the numerical scheme and (for linear) the max number of iterations
scheme.s0 = fun_set_s0(param);
scheme.N = 101;
% scheme.itmax = 10;

% Finite Difference strong form (nonlinear method) DA
param.cu = 1000;
[emode_FDS_DA_nonlin_1000,ds_FDS_DA_nonlin_1000] = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 100;
[emode_FDS_DA_nonlin_100 ,ds_FDS_DA_nonlin_100 ] = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 10;
[emode_FDS_DA_nonlin_10  ,ds_FDS_DA_nonlin_10  ] = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 1;
[emode_FDS_DA_nonlin_1   ,ds_FDS_DA_nonlin_1   ] = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 0.1;
[emode_FDS_DA_nonlin_0p1 ,ds_FDS_DA_nonlin_0p1 ] = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 0.01;
[emode_FDS_DA_nonlin_0p01,ds_FDS_DA_nonlin_0p01] = fun_Helm('FDS','DA','nonlin',param,scheme);

%% Plot the results on top of each other
figure(2); clf
subplot(4,1,1); hold on
plot(emode_FDS_DA_nonlin_1000.x, abs(ds_FDS_DA_nonlin_1000.t/emode_FDS_DA_nonlin_1000.M),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , abs(ds_FDS_DA_nonlin_100.t/emode_FDS_DA_nonlin_100.M)  ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , abs(ds_FDS_DA_nonlin_10.t/emode_FDS_DA_nonlin_10.M)    ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , abs(ds_FDS_DA_nonlin_1.t/emode_FDS_DA_nonlin_1.M)      ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , abs(ds_FDS_DA_nonlin_0p1.t/emode_FDS_DA_nonlin_0p1.M)  ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, abs(ds_FDS_DA_nonlin_0p01.t/emode_FDS_DA_nonlin_0p01.M),'Color',colmap(6,:),'LineWidth',1)
ylabel('$\partial s / \partial \tau$','Interpreter','Latex')
set(gca,'YLim',[0, 20])
subplot(4,1,2); hold on
plot(emode_FDS_DA_nonlin_1000.x, abs(ds_FDS_DA_nonlin_1000.h/emode_FDS_DA_nonlin_1000.M),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , abs(ds_FDS_DA_nonlin_100.h/emode_FDS_DA_nonlin_100.M)  ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , abs(ds_FDS_DA_nonlin_10.h/emode_FDS_DA_nonlin_10.M)    ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , abs(ds_FDS_DA_nonlin_1.h/emode_FDS_DA_nonlin_1.M)      ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , abs(ds_FDS_DA_nonlin_0p1.h/emode_FDS_DA_nonlin_0p1.M)  ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, abs(ds_FDS_DA_nonlin_0p01.h/emode_FDS_DA_nonlin_0p01.M),'Color',colmap(6,:),'LineWidth',1)
ylabel('$\partial s / \partial v_p$','Interpreter','Latex')
set(gca,'YLim',[0, 0.2])
xlabel('$x$','Interpreter','Latex')
subplot(4,1,3); hold on
plot(emode_FDS_DA_nonlin_1000.x, abs(ds_FDS_DA_nonlin_1000.wr/emode_FDS_DA_nonlin_1000.M),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , abs(ds_FDS_DA_nonlin_100.wr/emode_FDS_DA_nonlin_100.M)  ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , abs(ds_FDS_DA_nonlin_10.wr/emode_FDS_DA_nonlin_10.M)    ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , abs(ds_FDS_DA_nonlin_1.wr/emode_FDS_DA_nonlin_1.M)      ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , abs(ds_FDS_DA_nonlin_0p1.wr/emode_FDS_DA_nonlin_0p1.M)  ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, abs(ds_FDS_DA_nonlin_0p01.wr/emode_FDS_DA_nonlin_0p01.M),'Color',colmap(6,:),'LineWidth',1)
ylabel('$\partial s / \partial w_\rho$','Interpreter','Latex')
set(gca,'YLim',[0, 0.5])
subplot(4,1,4); hold on
plot(emode_FDS_DA_nonlin_1000.x, abs(ds_FDS_DA_nonlin_1000.v/emode_FDS_DA_nonlin_1000.M),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , abs(ds_FDS_DA_nonlin_100.v/emode_FDS_DA_nonlin_100.M)  ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , abs(ds_FDS_DA_nonlin_10.v/emode_FDS_DA_nonlin_10.M)    ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , abs(ds_FDS_DA_nonlin_1.v/emode_FDS_DA_nonlin_1.M)      ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , abs(ds_FDS_DA_nonlin_0p1.v/emode_FDS_DA_nonlin_0p1.M)  ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, abs(ds_FDS_DA_nonlin_0p01.v/emode_FDS_DA_nonlin_0p01.M),'Color',colmap(6,:),'LineWidth',1)
ylabel('$\partial s / \partial \overline{r}$','Interpreter','Latex')
set(gca,'YLim',[0, 10])
xlabel('$x$','Interpreter','Latex')

end
