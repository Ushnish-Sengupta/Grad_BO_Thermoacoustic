function [] = Sup_001()
% Sup_001
%
% Plot the direct and adjoint eigenfunctions as a function of the arbitrary
% parameter param.cu, which is used in the imposition of the upstream
% boundary condition in the FD case. The direct eigenfunction is unaffected
% by cu but the boundary point of the discrete adjoint (DA) eigenfunction 
% is strongly affected by cu. 

colmap = figfun_colmap;

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');

%% Calculate the reference scales and the nondimensional parameters
param = fun_nondim(param_dim);

%% Set the starting value of s, the numerical scheme and max number of iterations
scheme.s0 = fun_set_s0(param);
scheme.N = 100;

%% Calculate the eigenvalues with different methods

% Finite Difference strong form (nonlinear method) DA
param.cu = 1000;
[emode_FDS_DA_nonlin_1000] = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 100;
[emode_FDS_DA_nonlin_100]  = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 10;
[emode_FDS_DA_nonlin_10]   = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 1;
[emode_FDS_DA_nonlin_1]    = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 0.1;
[emode_FDS_DA_nonlin_0p1]  = fun_Helm('FDS','DA','nonlin',param,scheme);
param.cu = 0.01;
[emode_FDS_DA_nonlin_0p01] = fun_Helm('FDS','DA','nonlin',param,scheme);

%% Plot the results on top of each other
figure(1); clf
subplot(6,1,1); hold on
plot(emode_FDS_DA_nonlin_1000.x, real(emode_FDS_DA_nonlin_1000.pr),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , real(emode_FDS_DA_nonlin_100.pr) ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , real(emode_FDS_DA_nonlin_10.pr)  ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , real(emode_FDS_DA_nonlin_1.pr)   ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , real(emode_FDS_DA_nonlin_0p1.pr) ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, real(emode_FDS_DA_nonlin_0p01.pr),'Color',colmap(6,:),'LineWidth',1)
ylabel('$p_r$','Interpreter','Latex')
set(gca,'YLim',[0,2])
subplot(6,1,2); hold on
plot(emode_FDS_DA_nonlin_1000.x, imag(emode_FDS_DA_nonlin_1000.pr),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , imag(emode_FDS_DA_nonlin_100.pr) ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , imag(emode_FDS_DA_nonlin_10.pr)  ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , imag(emode_FDS_DA_nonlin_1.pr)   ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , imag(emode_FDS_DA_nonlin_0p1.pr) ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, imag(emode_FDS_DA_nonlin_0p01.pr),'Color',colmap(6,:),'LineWidth',1)
ylabel('$p_i$','Interpreter','Latex')
set(gca,'YLim',[-0.2,0.2])
xlabel('$x$','Interpreter','Latex')
subplot(6,1,3); hold on
plot(emode_FDS_DA_nonlin_1000.x, real(emode_FDS_DA_nonlin_1000.pl),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , real(emode_FDS_DA_nonlin_100.pl) ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , real(emode_FDS_DA_nonlin_10.pl)  ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , real(emode_FDS_DA_nonlin_1.pl)   ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , real(emode_FDS_DA_nonlin_0p1.pl) ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, real(emode_FDS_DA_nonlin_0p01.pl),'Color',colmap(6,:),'LineWidth',1)
ylabel('$p^\dag_r$','Interpreter','Latex')
set(gca,'YLim',[0,2])
subplot(6,1,4); hold on
plot(emode_FDS_DA_nonlin_1000.x, imag(emode_FDS_DA_nonlin_1000.pl),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , imag(emode_FDS_DA_nonlin_100.pl) ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , imag(emode_FDS_DA_nonlin_10.pl)  ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , imag(emode_FDS_DA_nonlin_1.pl)   ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , imag(emode_FDS_DA_nonlin_0p1.pl) ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, imag(emode_FDS_DA_nonlin_0p01.pl),'Color',colmap(6,:),'LineWidth',1)
ylabel('$p^\dag_i$','Interpreter','Latex')
set(gca,'YLim',[-0.2,0.2])
xlabel('$x$','Interpreter','Latex')
subplot(6,1,5); hold on
plot(emode_FDS_DA_nonlin_1000.x, real(emode_FDS_DA_nonlin_1000.pl),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , real(emode_FDS_DA_nonlin_100.pl) ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , real(emode_FDS_DA_nonlin_10.pl)  ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , real(emode_FDS_DA_nonlin_1.pl)   ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , real(emode_FDS_DA_nonlin_0p1.pl) ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, real(emode_FDS_DA_nonlin_0p01.pl),'Color',colmap(6,:),'LineWidth',1)
ylabel('$p^\dag_r$','Interpreter','Latex')
set(gca,'YLim',[0,0.2])
set(gca,'XLim',[0 0.01])
xlabel('$x$ (zoom)','Interpreter','Latex')
subplot(6,1,6); hold on
plot(emode_FDS_DA_nonlin_1000.x, imag(emode_FDS_DA_nonlin_1000.pl),'Color',colmap(1,:),'LineWidth',6)
plot(emode_FDS_DA_nonlin_100.x , imag(emode_FDS_DA_nonlin_100.pl) ,'Color',colmap(2,:),'LineWidth',5)
plot(emode_FDS_DA_nonlin_10.x  , imag(emode_FDS_DA_nonlin_10.pl)  ,'Color',colmap(3,:),'LineWidth',4)
plot(emode_FDS_DA_nonlin_1.x   , imag(emode_FDS_DA_nonlin_1.pl)   ,'Color',colmap(4,:),'LineWidth',3)
plot(emode_FDS_DA_nonlin_0p1.x , imag(emode_FDS_DA_nonlin_0p1.pl) ,'Color',colmap(5,:),'LineWidth',2)
plot(emode_FDS_DA_nonlin_0p01.x, imag(emode_FDS_DA_nonlin_0p01.pl),'Color',colmap(6,:),'LineWidth',1)
ylabel('$p^\dag_i$','Interpreter','Latex')
set(gca,'YLim',[0,0.1])
set(gca,'XLim',[0 0.01])
xlabel('$x$ (zoom)','Interpreter','Latex')

end
