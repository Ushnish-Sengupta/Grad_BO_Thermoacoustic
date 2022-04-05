function [] = Fig_002()

%% Create the figure and table for the direct mode for the Rijke tube
figfun_adjoint_mode('Rijke',100,'nonlin','Fig_002_Rijke');

%% Create the figure and table for the direct mode for the Rocket
figfun_adjoint_mode('rocket',100,'linear','Fig_002_rocket');

%% Create the figure and table for the direct mode for the Gas Turbine
figfun_adjoint_mode('GT',100,'nonlin','Fig_002_GT');

end
