function [] = Fig_006()

%% Create the figure for the feedback sensitivity for the Rocket
figfun_feedback_sen_z('rocket',100,'linear','Fig_006_rocket');

%% Create the figure for the feedback sensitivity for the Rijke tube
figfun_feedback_sen_z('Rijke',100,'nonlin','Fig_006_Rijke');

%% Create the figure for the feedback sensitivity for the Gas Turbine
figfun_feedback_sen_z('GT',100,'nonlin','Fig_006_GT');

end
