function [] = Fig_001()

%% Create the figure and table for the direct mode for the Rijke tube
figfun_direct_mode('Rijke',100,10,'Fig_001_Rijke'); % print to file

%% Create the figure and table for the direct mode for the Gas Turbine
figfun_direct_mode('GT',100,10,'Fig_001_GT'); % print to file

%% Create the figure and table for the direct mode for the rocket
figfun_direct_mode('rocket',100,10,'Fig_001_rocket'); % print to file

end
