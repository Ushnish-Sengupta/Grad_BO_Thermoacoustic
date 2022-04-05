function [] = Fig_010()
% Fig_010
%
% Calculate the influence of drag from an adiabatic mesh

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');

% Set the dimensional flow speed
U0 = 0.10;   % m/s
% Set the pressure drop coefficient
K = 20; % dimensionless

%% Calculate the reference scales and the nondimensional parameters
[param,ref] = fun_nondim(param_dim);

%% Set the numerical scheme, the max number of iterations, and the starting value of s
scheme.N     = 100;
scheme.s0    = fun_set_s0(param);

%% Calculate the influence of drag from the heater
% Calculate the eigenmode and the internal sensitivities
[emode,ds_int,~] = fun_Helm('FEW','DA','nonlin',param,scheme);
% Extract the non-dimensional feedback sensitivity from velocity into force
% (in continuous form)
dsdfru_nondim = ds_int.fru/emode.M00;
% Convert to dimensional (rad/distance)
dsdfru_dim = dsdfru_nondim / ref.l_dim;
% Calculate the dimensional eigenvalue shift due to drag from the heater
% (units m/s * rad/m = rad/s)
dsdheater = - K * U0 * dsdfru_dim;
% Calculate the dimensional s (in radians/sec)
s_dim = emode.s * ref.u_dim / ref.l_dim;

%% Plot all in terms of cycles / sec
figure(1); clf
subplot(3,1,1); hold on
han_adj = plot(emode.x0*ref.l_dim,real(dsdheater),'k-');
xlabel('Position of adiabatic mesh, $X_m$ (m)','Interpreter','Latex')
ylabel('growth rate shift, $\delta s_r$ (rad s$^{-1}$)','Interpreter','Latex')
title(['Frequency = ',num2str(imag(s_dim),4),' rad s$^{-1}$; growth rate = ',num2str(real(s_dim),4),'  rad s$^{-1}$'],'Interpreter','Latex')
% text(0.1,-1,'Ugly units, to be consistent with Rigas and Jamieson')

%% Experimental results
x_exp = [    
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

ds_r_exp = [
   -0.6202
   -0.3188
   -0.1383
   -0.0606
   -0.1020
   -0.0216
   -0.1903
   -0.3987
   -0.7073
   -0.8606
   -1.1813
   -1.3943
   -1.4861
   -1.6321
];

han_exp = plot(x_exp,ds_r_exp,'ko-');

han_leg = legend([han_adj,han_exp],{'prediction from feedback sensitivity','measurement from Rigas (2016)'});
set(han_leg,'Interpreter','Latex','Location','South')
set(gca,'YLim',[-2.5,0])

box on

print('-depsc2','figures/Fig_010.eps')

end


