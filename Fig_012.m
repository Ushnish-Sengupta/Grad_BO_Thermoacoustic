function [] = Fig_012()
% Fig_012
%
% Influence of the viscous boundary layer

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');
% Set the kinematic viscosity 
nu = 1.48e-5; % m^2/s
% Set the tube diameter
D = 0.047; % m
% Calculate perimeter / area of the tube
PonA = 4/D; % 1/m

%% Calculate the reference scales and the nondimensional parameters
[param,ref] = fun_nondim(param_dim);

%% Set the numerical scheme, the max number of iterations, and the starting value of s
scheme.N     = 100;
scheme.s0    = fun_set_s0(param);

%% Calculate the influence of drag from the heater
% Calculate the eigenmode and the internal sensitivities
[emode,ds_int,~] = fun_Helm('FEW','DA','nonlin',param,scheme);

% Calculate the dimensional s (in radians/sec)
s_dim = emode.s * ref.u_dim / ref.l_dim;

% Extract the non-dimensional feedback sensitivity from velocity into force
% (in continuous form)
dsdfru_nondim = ds_int.fru/emode.M00;
% Convert to dimensional (rad/distance)
dsdfru_dim = dsdfru_nondim / ref.l_dim;

%% Calculate the thickness of the acoustic boundary layer
% (for oscillatory flow at angular frequency omega above a stationary boundary)
% delta = 2*pi*sqrt(2*nu/omega)
dbl = 2*pi*sqrt(2*nu/imag(s_dim));

% Calculate the local dimensional eigenvalue shift per unit length, due to the viscous bl
dsdX = - dsdfru_dim * nu/dbl * PonA;

% Integrate this
ds_tot = dsdX * emode.M00 * ones(length(dsdX),1);

%% Plot all in terms of cycles / sec
figure(1); clf
subplot(3,1,1); hold on
han_local = plot(emode.x0*ref.l_dim,real(dsdX),'k-');
han_global = plot([0,1],[real(ds_tot),real(ds_tot)],'k--');
text(0.8,0.9*real(ds_tot),['$\delta s_r = ',num2str(real(ds_tot),3),'$ rad s$^{-1}$'],'Interpreter','Latex')
xlabel('$X$ (m)','Interpreter','Latex')
ylabel('$\delta s_r / \delta X$ (rad s$^{-1}$ m$^{-1}$)','Interpreter','Latex')
title(['Frequency = ',num2str(imag(s_dim),4),' rad s$^{-1}$; growth rate = ',num2str(real(s_dim),4),'  rad s$^{-1}$'],'Interpreter','Latex')
han_leg = legend([han_local,han_global],{'local eigenvalue shift','total eigenvalue shift $\int (ds_r/dX) d X$'});
set(han_leg,'Interpreter','Latex','Location','South')

axis([0 1 -1.5 0])
box on

print('-depsc2','figures/Fig_012.eps')

end


