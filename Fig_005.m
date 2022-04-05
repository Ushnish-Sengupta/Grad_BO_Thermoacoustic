function [] = Fig_005()

% %% Create the figure of the receptivities for the Rijke tube
figfun_receptivities('Rijke',100,'nonlin','Fig_005_Rijke');

% %% Create the figure of the receptivities for the Gas Turbine
figfun_receptivities('GT',100,'nonlin','Fig_005_GT');

%% Create the figure of the receptivities for the rocket
figfun_receptivities('rocket',100,'linear','Fig_005_rocket');

end
