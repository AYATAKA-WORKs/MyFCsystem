clear
close all
%% Add Path
addpath("SLX\");
addpath("TEST\");
addpath("DATA\");

%% User Setting Parameter
% Operational Parameter --------------------------------------------------
% N_cmd.time = [0;10;15;20;25;120];
% N_cmd.signals.values = [10e3;10e3;70e3;70e3;100e3;100e3];
W_stack_cmd.time = [0;20;120];
W_stack_cmd.signals.values = [0;0.075;0.075];
p_stack_cmd.time = [0;20;120];
p_stack_cmd.signals.values = [101325;250000;250000];

% Boundary Condition -----------------------------------------------------
T_amb = 293.15;         % [K]
p_amb = 101325;         % [Pa]
RH_amb = 0.5;           % [-]
yO2_amb = 20.9476e-2;   % [-]
% ※ 環境条件=初期値 とする

% Component Parameter ----------------------------------------------------
% Compressor
cmp.param.dc = 0.2286;      % [m]
cmp.param.eta_is = 0.8;     % [-]
cmp.param.eta_mech = 0.98;  % [-]
cmp.map.coeff.a = [-3.69906e-5, 2.70399e-4, -5.36235e-4, -4.63685e-5, 2.21195e-3];
cmp.map.coeff.b = [1.766467,-1.34837,2.44419];
cmp.map.coeff.c = [-9.78755e-3, 0.10581, -0.42937, 0.80121, -0.68344, 0.43331];
% Supply manifold
V_sm = 1.00e-2;             % [m^3]
A_sm = 1.00e-3;             % [m^2]
% Heatexchanger
eta_hex = 0.8;              % [-]
T_coolant = 293.15;         % [K]
cp_coolant = 4185;          % [J/k/kg]
% FC stack
lambda = 1.3;               % Excess O2 ratio [-]
V_stack = 3.75e-3;          % [m^3]
A_stack = 0.001^2*3750/10;  % [m^2]
% Turbine
V_tbn_m = 0.01;             % [m^3]
A_VGS_center = 0.012^2;     % [m^2]
tbn.param.eta_is = 0.8;     % [-]
tbn.param.eta_mech = 0.98;  % [-]
% Motor
k_t_M = 0.0153;             % [N*m\A]
k_e_M = 0.0153;             % [V/(rad/s}]
R_M   = 0.82;               % [Ω]
eta_M = 0.98;               % [-]
% Rotor
J_rotor = 0.001;            % [kg*m^2]


%% Automatically calculated parameter
% Initial Condition
T_init = T_amb;             % [K]
p_init = p_amb;             % [Pa]
RH_init = RH_amb;           % [-]
yO2_init = yO2_amb;         % [-]
% Gas properties
gsprop = GasProperties;     % generate instance
gsprop = gsprop.UpdateProperties(p_init, T_init, RH_init, yO2_init);
% Molar mass
M_init = gsprop.M_MA;       % [kg/mol]
M_Nmix = gsprop.M_Nmix;     % [kg/mol]
M_O2   = gsprop.M_O*2;      % [kg/mol]
M_CO2  = gsprop.M_C + gsprop.M_O*2; % [kg/mol]
M_H2O  = gsprop.M_H*2 + gsprop.M_O; % [kg/mol]
% Gas constant
R_init = gsprop.R_MA;       % [J/K/kg]
R_Nmix = gsprop.R_Nmix;     % [J/K/kg]
R_O2   = gsprop.R0/M_O2;    % [J/K/kg]
R_vp   = gsprop.R0/M_H2O;   % [J/K/kg]
% Initial pressure
p_vp_init = sat_vp_pressure(T_init)*RH_init;    % [Pa]
p_DA_init = p_init - p_vp_init;                 % [Pa]
p_O2_init = yO2_init * p_DA_init;               % [Pa]
p_Nmix_init = (1-yO2_init) * p_DA_init;         % [Pa]
% Initial mass of each volume
m_st_Nmix_init = p_Nmix_init*V_stack/(R_Nmix*T_init);   % [kg]
m_st_O2_init = p_O2_init*V_stack/(R_O2*T_init);         % [kg]
m_st_vp_init = p_vp_init*V_stack/(R_vp*T_init);         % [kg]

clear gsprop

%% Calculate saturated vapor pressure [Pa]
function psat = sat_vp_pressure(Ts)
    P_CONVERT = 0.001;
    C1 = -5.6745359e3;
    C2 = 6.3925247;
    C3 = -9.6778430e-3;
    C4 = 6.2215701e-7;
    C5 = 2.0747825e-9;
    C6 = -9.4840240e-13;
    C7 = 4.1635019;
    N1 = 0.11670521452767e4;
    N2 = -0.72421316703206e6;
    N3 = -0.17073846940092e2;
    N4 = 0.12020824702470e5;
    N5 = -0.32325550322333e7;
    N6 = 0.14915108613530e2;
    N7 = -0.4823265731591e4;
    N8 = 0.40511340542057e6;
    N9 = -0.23855557567849e0;
    N10 = 0.65017534844798e3;
    
    if Ts < 273.15+0.01
        % -100~0.01C//三重点を計算以下は wexler-hyland のシミュレーションプログラム式
        psat = exp(C1 / Ts + C2 + C3 * Ts + C4 * Ts^2 + C5 * Ts^3 + C6 * Ts^4 + C7 * log(Ts)) * P_CONVERT;
    else
        % ~647.096K//臨界温度まで一定とするは IAPWS-IF97 実用国際状態式
        alpha = Ts + N9 / (Ts - N10);
        a2 = alpha * alpha;
        A = a2 + N1 * alpha + N2;
        B = N3 * a2 + N4 * alpha + N5;
        C = N6 * a2 + N7 * alpha + N8;
        psat = (2 * C / (-B + (B * B - 4 * A * C)^0.5))^4 / P_CONVERT;
    end 
    psat = psat*1000; % kPa to Pa
end
