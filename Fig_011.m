function [] = Fig_011()
% Fig_010
%
% Calculate the influence of heat input from a hot mesh

% Influence of heat transfer from a secondary heater

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');

% Set the power of the secondary heater
Q = 200;

%% Calculate the reference scales and the nondimensional parameters
[param,ref] = fun_nondim(param_dim);

%% Set the numerical scheme, the max number of iterations, and the starting value of s
scheme.N     = 100;
scheme.s0    = fun_set_s0(param);

% Calculate the eigenmode and the internal sensitivities
[emode,ds_int,~] = fun_Helm('FEW','DA','nonlin',param,scheme);
% Calculate the dimensional s (in radians/sec)
s_dim = emode.s * ref.u_dim / ref.l_dim;

% Calculate nw by using the same FTF and U as the main heater
nw_dim = Q * param_dim.FTF/param_dim.U;
% Set the tube area 
A_dim = param_dim.S_c; % tube x-sect area
tau_dim = param_dim.tau;

% Extract the non-dimensional feedback sensitivity from velocity into heat
% input (in continuous nondim form)
dsdqpu_nondim = ds_int.qpu/emode.M00;

% Calculate ds
ds = (ref.u_dim/ref.l_dim) * nw_dim / (ref.p_dim * A_dim) * exp(-s_dim*tau_dim) * dsdqpu_nondim;

%% Plot
figure(1); clf
subplot(3,1,1); hold on
han_adj = plot(emode.x0*ref.l_dim,real(ds)/Q,'k-');
xlabel('Position of hot mesh, $X_m$ (m)','Interpreter','Latex')
ylabel('growth rate shift, $\delta s_r$ (rad s$^{-1} W^{-1}$)','Interpreter','Latex')
title(['Frequency = ',num2str(imag(s_dim),4),' rad s$^{-1}$; growth rate = ',num2str(real(s_dim),4),'  rad s$^{-1}$'],'Interpreter','Latex')
% text(0.1,-1,'Ugly units, to be consistent with Rigas and Jamieson')

x_exp = [
    0.0500
    0.1000
    0.1500
    0.2000
    0.3000
    0.3500
    0.4000
    0.4500
    0.5000
    0.5500
    0.6000
    0.6500
    0.7000
    0.7500
    0.8000
    0.8500
    0.9000
    0.9500
];

y_exp = [
   -0.0504
   -0.0248
   -0.0028
    0.0190
    0.0376
    0.0337
    0.0108
   -0.0117
   -0.0336
   -0.0597
   -0.0767
   -0.0940
   -0.0988
   -0.1032
   -0.1022
   -0.0891
   -0.0765
   -0.0471
];

han_exp = plot(x_exp,y_exp,'ko-');

han_leg = legend([han_adj,han_exp],{'prediction from feedback sensitivity','measurement from Jamieson (2017)'});
set(han_leg,'Interpreter','Latex','Location','SouthWest')
box on 

print('-depsc2','figures/Fig_011.eps')


end


