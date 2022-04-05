function [] = Tab_001()
% Tab_001
% 
% Create the data and figures for table 1

%% Create the baseflow figure and table for the Rijke tube
tabfun_baseflow('Rijke','Tab_001_Rijke'); % print to file

%% Create the baseflow figure and table for the Rijke tube
tabfun_baseflow('GT','Tab_001_GT'); % print to file

%% Create the baseflow figure and table for the rocket
tabfun_baseflow('rocket','Tab_001_rocket'); % print to file

end

