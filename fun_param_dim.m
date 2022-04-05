function [param_dim] = fun_param_dim(model)
% fun_param_dim
%
% Set the dimensional parameters
%
% INPUTS: 
% model     'Rijke' | 'GT' | 'play' | 'rocket'
%
% OUTPUTS:
% param_dim  structure containing the dimensional parameters

%% If no model specified, set to 'play'
if ~nargin
    model = 'play'; % 'Rijke' | 'GT' | 'rocket' | 'play'
end

switch model
    case 'play'
        % Gas parameters
        param_dim.gam    = 1.4;     % ratio of specific heats (non-dim)
        % Acoustic chamber parameters
        param_dim.rhobar = 1.22;    % mean density (kg/m3)
        param_dim.pbar   = 1e5;     % mean pressure (Pascals)
        param_dim.X      = 1.0;     % length of chamber (metres)
        param_dim.S_c    = pi*0.047^2/4; % Cross-sectional area of chamber (metres^2)
        % Heater parameters, which are used internally to calculate param_dim.n
        param_dim.Q      = 10000;                    % mean heat release rate (Joules per second)
        param_dim.U      = 0.10;                   % mean velocity through tube (metres per second)
        param_dim.FTF    = 0.014;                % Flame/heater Transfer Function (non-dim)
        param_dim.eta    = param_dim.FTF*param_dim.Q/param_dim.U; % 3D heat release rate fluctuation / U (Joules/metre)
        param_dim.n      = param_dim.eta/param_dim.pbar/param_dim.S_c; % 1D heat release rate fluctuation / U (no units)
        param_dim.tau    = 0.0015;  % heat release rate time delay (seconds)
        % Heater position parameters
        param_dim.X_w    = 0.250;    % position of measurement zone (metres)
        param_dim.L_w    = 0.05;   % length of measurement zone (metres)
        param_dim.X_h    = 0.25;    % position of heat release region (metres)
        param_dim.L_h    = 0.05;   % length of heat release region (metres)
        % Boundary conditions
        param_dim.Ru     = complex(-0.975,0.05);  % reflection coefficient at upstream end
        param_dim.Rd     = complex(-0.975,0.05);  % reflection coefficient at downstream end
        param_dim.rhu    = param_dim.rhobar; % upstream density (kg/m3)
        param_dim.rhd    = 0.85;     % downstream density (kg/m3)
        % Arbitrary boundary condition coefficients
        param_dim.cu     = 100;      % Arbitrary weighting for upstream boundary condition (used in FD only)
        param_dim.cd     = 100;      % Arbitrary weighting for downstream boundary condition (used in FD only)
        % Starting point for iteration
        param_dim.mode   = 1;        % (integer) acoustic mode to start from
    case 'Rijke' % Based on the Rijke tube in Rigas et al (2016)
        % Gas parameters
        param_dim.gam    = 1.4;     % ratio of specific heats (non-dim)
        % Acoustic chamber parameters
        param_dim.rhobar = 1.22;    % mean density (kg/m3)
        param_dim.pbar   = 1e5;     % mean pressure (Pascals)
        param_dim.X      = 1.0;     % length of chamber (metres)
        param_dim.S_c    = pi*0.047^2/4; % Cross-sectional area of chamber (metres^2)
        % Heater parameters, which are used internally to calculate param_dim.n
        param_dim.Q      = 200;                    % mean heat release rate (Joules per second)
        param_dim.U      = 0.10;                   % mean velocity through tube (metres per second)
        param_dim.FTF    = 0.014;                % Flame/heater Transfer Function (non-dim)
        param_dim.eta    = param_dim.FTF*param_dim.Q/param_dim.U; % 3D heat release rate fluctuation / U (Joules/metre)
        param_dim.n      = param_dim.eta/param_dim.pbar/param_dim.S_c; % 1D heat release rate fluctuation / U (no units)
        param_dim.tau    = 0.0015;  % heat release rate time delay (seconds)
        % Heater position parameters
        param_dim.X_w    = 0.20;    % position of measurement zone (metres)
        param_dim.L_w    = 0.025;   % length of measurement zone (metres)
        param_dim.X_h    = 0.25;    % position of heat release region (metres)
        param_dim.L_h    = 0.025;   % length of heat release region (metres)
        % Boundary conditions
        param_dim.Ru     = complex(-0.975,0.05);  % reflection coefficient at upstream end
        param_dim.Rd     = complex(-0.975,0.05);  % reflection coefficient at downstream end
        param_dim.rhu    = param_dim.rhobar; % upstream density (kg/m3)
        param_dim.rhd    = 0.85;     % downstream density (kg/m3)
        % Arbitrary boundary condition coefficients
        param_dim.cu     = 100;      % Arbitrary weighting for upstream boundary condition (used in FD only)
        param_dim.cd     = 100;      % Arbitrary weighting for downstream boundary condition (used in FD only)
        % Starting point for iteration
        param_dim.mode   = 1;        % (integer) acoustic mode to start from
    case 'GT'
        % Based on a typical gas turbine at take-off:
        %
        % Atmospheric pressure 1 bar, atmospheric temperature 290K
        % ==> ambient density = p/RT = 1.2015
        % isentropic compression over pressure ratio 40
        % ==> chamber pressure = 40 bar
        % ==> chamber inlet density = 1.2015 * 40^(1/gam) = 16.75 kg/m3
        % chamber length = 0.2m
        % chamber heat release = 20MW
        % mean speed around measurement point = 100 m/s
        % chamber exit temperature = 1800K
        % ==> chamber exit density = 40e5/287/1800 = 7.74 kg/m3
        % chamber exit Mach number = 0.98;
        %
        % Flame parameters, which are used internally to calculate param_dim.n
        Q = 20e6;                    % mean heat release rate (Joules per second)
        U = 100;                     % mean velocity through tube (metres per second)
        FTF = 1.;                    % Flame/heater Transfer Function (non-dim)
        % Acoustic chamber parameters
        param_dim.rhobar = 16.75;    % mean density (kg/m3)
        param_dim.pbar   = 40e5;     % mean pressure (Pascals)
        param_dim.gam    = 1.4;      % ratio of specific heats (non-dim)
        param_dim.X      = 0.2;      % length of chamber (metres)
        param_dim.S_c    = pi*(0.7^2 - 0.6^2); % Cross-sectional area of chamber (metres^2)
        param_dim.X_w    = 0.02;     % position of measurement zone (metres)
        param_dim.L_w    = 0.005;    % length of measurement zone (metres)
        param_dim.X_h    = 0.1;      % position of heat release region (metres)
        param_dim.L_h    = 0.05;     % length of heat release region (metres)
        param_dim.Ru     = -1;  % reflection coefficient at upstream end
        param_dim.tau    = 0.0003;   % heat release rate time delay (seconds)
        param_dim.rhu    = param_dim.rhobar; % upstream density (kg/m3)
        param_dim.rhd    = 7.74;     % downstream density (kg/m3)
        param_dim.eta    = FTF*Q/U; % 3D heat release rate fluctuation / U (Joules/metre)
        param_dim.n      = param_dim.eta/param_dim.pbar/param_dim.S_c; % 1D heat release rate fluctuation / U (no units)
        param_dim.cu     = .1;       % Arbitrary weighting for upstream boundary condition (used in FD only)
        param_dim.cd     = .1;       % Arbitrary weighting for downstream boundary condition (used in FD only)
        % Downstream reflection coefficient from Marble and Candel JSV 1977 55(2), 225-243
        Machd            = 0.98;     % Mach number at exit from combustor
        param_dim.Rd     = (1-0.5*(param_dim.gam-1)*Machd) ...
                         / (1+0.5*(param_dim.gam-1)*Machd);
        % Starting point for iteration
        param_dim.mode   = 0;        % (integer) acoustic mode to start from
    case 'rocket'
        % Based on a H2/LOx rocket engine (Vulcain of Ariane V)
        %
        % Chamber pressure 100 bar (10 MPa), 
        P_c = 10e6;
        % Chamber temperature 2000K
        T_c = 2000;
        % Chamber gases: 50% H2O, 50% H2 by mole
        vol_H2O = 0.50;
        vol_H2  = 1-vol_H2O;
        % RMM of chamber gases
        RMM_H2O = 18;
        RMM_H2  = 2;
        % R_g of chamber gases
        R_g = vol_H2O * 8314/RMM_H2O + vol_H2 * 8314/RMM_H2;
        % Density of chamber gases (assume uniform through chamber)
        rho_c = P_c / T_c / R_g;
        % gamma
        gam_H2O = 1.31;
        gam_H2 = 1.41;
        gam_c = vol_H2O * gam_H2O + vol_H2 * gam_H2;
        % Heat release rate in chamber (2GW)
        Q_c = 2e9; % Watts
        % weighted flow speed at injector (4m/s for LOx, 100 m/s for H2 at entry)
        U_inj = 10;
        %
        % Flame parameters, which are used internally to calculate param_dim.n
        Q = Q_c;                     % mean heat release rate (Joules per second)
        U = U_inj;                     % mean velocity through tube (metres per second)
        FTF = 1.;                    % Flame/heater Transfer Function (non-dim)
        % Acoustic chamber parameters
        param_dim.rhobar = rho_c;    % mean density (kg/m3)
        param_dim.pbar   = P_c;      % mean pressure (Pascals)
        param_dim.gam    = gam_c;    % ratio of specific heats (non-dim)
        param_dim.X      = 1.0;      % length of chamber (metres)
        param_dim.S_c    = pi*(0.4^2); % Cross-sectional area of chamber (metres^2)
        param_dim.X_w    = 0.06;     % position of measurement zone (metres)
        param_dim.L_w    = 0.02;     % length of measurement zone (metres)
        param_dim.X_h    = 0.7;      % position of heat release region (metres)
        param_dim.L_h    = 0.1;      % length of heat release region (metres)
        param_dim.Ru     = +0.999;   % reflection coefficient at upstream end
        param_dim.tau    = 0.0001;   % heat release rate time delay (seconds)
        param_dim.rhu    = param_dim.rhobar; % upstream density (kg/m3)
        param_dim.rhd    = rho_c;     % downstream density (kg/m3)
        param_dim.eta    = FTF*Q/U;  % 3D heat release rate fluctuation / U (Joules/metre)
        param_dim.n      = param_dim.eta/param_dim.pbar/param_dim.S_c; % 1D heat release rate fluctuation / U (no units)
        param_dim.cu     = .1;       % Arbitrary weighting for upstream boundary condition (used in FD only)
        param_dim.cd     = .1;       % Arbitrary weighting for downstream boundary condition (used in FD only)
        % Downstream reflection coefficient from Marble and Candel JSV 1977 55(2), 225-243
        Machd            = 1.00;     % Mach number at exit from combustor
        param_dim.Rd     = (1-0.5*(param_dim.gam-1)*Machd) ...
                         / (1+0.5*(param_dim.gam-1)*Machd);
        % Starting point for iteration
        param_dim.mode   = 1;        % (integer) acoustic mode to start from
end

