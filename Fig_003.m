function [] = Fig_003()

%% Create the figure and table for the base state sensitivities of the rocket
figfun_base_sen('rocket',100,'linear','Fig_003_rocket');

%% Create the figure and table for the base state sensitivities of the Rijke tube
figfun_base_sen('Rijke',100,'nonlin','Fig_003_Rijke');

%% Create the figure and table for the base state sensitivities of the Gas Turbine
figfun_base_sen('GT',100,'nonlin','Fig_003_GT');

end
