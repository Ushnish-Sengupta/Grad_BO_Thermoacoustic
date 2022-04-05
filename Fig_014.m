function [] = Fig_014()
% Fig_014
%
% Influence of a Helmholtz Resonator

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');
% Set the tube diameter
D = 0.047; % m
% Calculate the area of the tube
A_dim = pi*D^2/4; % m^2

%% Calculate the reference scales and the nondimensional parameters
[param,ref] = fun_nondim(param_dim);

%% Set the numerical scheme, the max number of iterations, and the starting value of s
scheme.N     = 100;
scheme.s0    = fun_set_s0(param);

%% Perform the eigenvalue calculation and output the sensitivities
% Calculate the eigenmode and the internal sensitivities
[emode,ds_int,~] = fun_Helm('FEW','DA','nonlin',param,scheme);
% Calculate the dimensional s (in radians/sec)
s_dim = emode.s * ref.u_dim / ref.l_dim;
% Extract the oscillation anular frequency
w_dim = imag(s_dim);
disp(['Frequency of oscillations is ',num2str(w_dim/2/pi),' Hz'])

%% Calculate the influence of the Helmholtz resonator
% 
% Mean velocity through neck of Helmholtz resonator
U_neck_dim = 10; % m/s
% Cross-sectional area of throat of Helmholtz resonator 
A_neck_dim = 0.0001; % m^2
% Length of neck of Helmholtz resonator 
L_neck_dim = 0.05; % metres
% Volume of Helmholtz resonator
V_helm_dim = 0.0001; % metres^3
% End correction factor
d_dim = 0.61*L_neck_dim;
%
% Natural frequency of resonator
w0_dim = sqrt(A_neck_dim * ref.u_dim^2 / (L_neck_dim + 2*d_dim) / V_helm_dim);
disp(['Natural frequency of Helmholtz resonator is ',num2str(w0_dim/2/pi),' Hz'])
% Natural frequency of neck of resonator
w1_dim = ref.u_dim / (L_neck_dim + 2*d_dim);
disp(['Natural frequency of neck of resonator is ',num2str(w1_dim/2),' Hz'])
%
% Mach number of flow in neck of Helmholtz resonator
M_neck = U_neck_dim / ref.u_dim;
% Calculate 1/(dimensionless) complex impedance of resonator
% Z = (1 - (w_dim/w0_dim)^2 + 1i*M_neck*w1_dim*w_dim/w0_dim^2) ...
%    /(M_neck + 1i*w1_dim*w_dim/w0_dim^2);
% disp(['Complex impedance = ',num2str(real(Z),3),' + ',num2str(imag(Z),3),'i'])
Zm1 = (M_neck + 1i*w1_dim*w_dim/w0_dim^2) ...
/ (1 - (w_dim/w0_dim)^2 + 1i*M_neck*w1_dim*w_dim/w0_dim^2);

% Extract the non-dimensional feedback sensitivity from velocity into force
% (in continuous form)
mrp1_nondim = ds_int.mrp/emode.M11;

% Note that ref.p_dim = pbar and c0 = sqrt(gamma*p/rho)

% Calculate the local dimensional eigenvalue shift per unit length, due to
% the acoustic liner
ds = - (ref.u_dim/param.gam/ref.l_dim) * A_neck_dim/A_dim * Zm1 * mrp1_nondim;

%% Plot ds as a function of X in terms of radians / sec
figure(1); clf
subplot(3,1,1); hold on
plot(emode.x1*ref.l_dim,real(ds),'k-');
xlabel('$X_n$ (m)','Interpreter','Latex')
ylabel('growth rate shift, $\delta s_r$ (rad s$^{-1}$)','Interpreter','Latex')
title(['Frequency = ',num2str(imag(s_dim),4),' rad s$^{-1}$; growth rate = ',num2str(real(s_dim),4),'  rad s$^{-1}$'],'Interpreter','Latex')

axis([0 1 -80 0])
box on

print('-depsc2','figures/Fig_014.eps')

end


