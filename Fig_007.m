function [] = Fig_007()

%% Create the figure for the base state sensitivity for the Rijke tube
figfun_base_sen_z('Rijke',100,'nonlin','Fig_007_Rijke');

%% Create the figure for the base state sensitivity for the Gas Turbine
figfun_base_sen_z('GT',100,'nonlin','Fig_007_GT');

%% Create the figure for the base state sensitivity for the Rocket
figfun_base_sen_z('rocket',100,'linear','Fig_007_rocket');

end
