clear all
close all

addpath("..\SLX\");
addpath("..\TEST\");

%% User Setting Parameter

% Operational Parameter
A_VGS = 0.015^2;    % [m^2]
N = 70000;          % [rpm]

% Boundary Condition
T_in  = 373.15;     % [K]
p_in  = 200000;     % [Pa]
p_out = 101325;     % [Pa]
RH_in = 0.5;        % [-]
cp    = 1005;       % [J/kg/K]

% Geometoric Parameter

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
p_vp_in = sat_vp_pressure(T_in)*RH_in;   % [Pa]
p_DA_in = p_in - p_vp_in;               % [Pa]
p_O2_in = molf_O2_DA * p_DA_in;           % [Pa]
p_N2_in = molf_N2_DA * p_DA_in;           % [Pa]
molf_O2_MA_in = p_O2_in/p_in;
molf_N2_MA_in = p_N2_in/p_in;
molf_H2O_MA_in = p_vp_in/p_in;
M_MA_in = molf_O2_MA_in * M_O2 + molf_N2_MA_in * M_N2 + molf_H2O_MA_in * M_H2O;
R_MA_in = R0/M_MA_in;

% Initial Condition
T_init = T_in;             % [K]
p_init = p_in;             % [Pa]
RH_init = RH_in;           % [-]
p_vp_init = p_vp_in;       % [Pa]
p_DA_init = p_DA_in;       % [Pa]
p_O2_init = p_O2_in;       % [Pa]
p_N2_init = p_N2_in;       % [Pa]

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
