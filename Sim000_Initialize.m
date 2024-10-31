clear all
close all

addpath("SLX\");
addpath("TEST\");

%% User Setting Parameter

% Operational Parameter
lambda = 1.3;

% Boundary Condition
T_amb = 293.15;    % [K]
p_amb = 101325;    % [Pa]
RH_amb = 0.5;      % [-]

% Geometoric Parameter
% FC stack
V_stack = 3.75e-3;          % [m^3]
A_stack = 0.001^2*37500;    % [m^2]

%% Automatic calculation parameter
% Physical property
R0 = 8.31446261815324; % [J/K/mol]
R = 287;    % [J/kg/K]
kap = 1.4;  % [-]

% Dry air property
% Mol fraction
molf_O2_DA = 0.2;
molf_N2_DA = 0.8;
yO2_amb = molf_O2_DA;
% Molecular weight [kg/mol]
M_O2 = 32e-3;                               % [kg/mol] 
M_N2 = 28e-3;                               % [kg/mol]
M_DA = molf_O2_DA*M_O2 + molf_N2_DA*M_N2;   % [kg/mol] - Dry air molecular weight

% Gas Constant
R_O2 = R0/M_O2;     % [J/kg/K]
R_N2 = R0/M_N2;     % [J/kg/K]
R_DA = R0/M_DA;     % [J/kg/K]

% Vapor property
% Molecular weight
M_H2O = 18e-3;      % [kg/mol]
% Gas Constant
R_vp = R0/M_H2O;    % [J/kg/K]

% Moist air property
p_vp_amb = sat_vp_pressure(T_amb)*RH_amb;   % [Pa]
p_DA_amb = p_amb - p_vp_amb;               % [Pa]
p_O2_amb = molf_O2_DA * p_DA_amb;           % [Pa]
p_N2_amb = molf_N2_DA * p_DA_amb;           % [Pa]
molf_O2_MA_amb = p_O2_amb/p_amb;
molf_N2_MA_amb = p_N2_amb/p_amb;
molf_H2O_MA_amb = p_vp_amb/p_amb;
M_MA_amb = molf_O2_MA_amb * M_O2 + molf_N2_MA_amb * M_N2 + molf_H2O_MA_amb * M_H2O;
R_MA_amb = R0/M_MA_amb;

% Initial Condition
T_init = T_amb;             % [K]
p_init = p_amb;             % [Pa]
RH_init = RH_amb;           % [-]
p_vp_init = p_vp_amb;       % [Pa]
p_DA_init = p_DA_amb;       % [Pa]
p_O2_init = p_O2_amb;       % [Pa]
p_N2_init = p_N2_amb;       % [Pa]
t); 
m_st_O2_init = p_O2_init*V_stack/(R_O2*T_init);
m_st_vp_init = p_vp_init*V_stack/(R_vp*T_init);

% Set Compressor Parameter
cmp.param.dc = 0.2286;      % [m]
cmp.param.eta_is = 0.8;     % [-]
cmp.param.eta_mech = 0.98;  % [-]
cmp.map.coeff.a = [-3.69906e-5, 2.70399e-4, -5.36235e-4, -4.63685e-5, 2.21195e-3];
cmp.map.coeff.b = [1.766467,-1.34837,2.44419];
cmp.map.coeff.c = [-9.78755e-3, 0.10581, -0.42937, 0.80121, -0.68344, 0.43331];

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
        %　~647.096K//臨界温度まで一定とするは IAPWS-IF97 実用国際状態式
        alpha = Ts + N9 / (Ts - N10);
        a2 = alpha * alpha;
        A = a2 + N1 * alpha + N2;
        B = N3 * a2 + N4 * alpha + N5;
        C = N6 * a2 + N7 * alpha + N8;
        psat = (2 * C / (-B + (B * B - 4 * A * C)^0.5))^4 / P_CONVERT;
    end 
    psat = psat*1000; %kPa to Pa
end
