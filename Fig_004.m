function [] = Fig_004()

%% Create the figure for the feedback sensitivity for the Rijke tube
figfun_feedback_sen('Rijke',100,'nonlin','Fig_004_Rijke');

%% Create the figure for the feedback sensitivity for the Gas Turbine
figfun_feedback_sen('GT',100,'nonlin','Fig_004_GT');

%% Create the figure for the feedback sensitivity for the rocket
figfun_feedback_sen('rocket',100,'linear','Fig_004_rocket');

end
