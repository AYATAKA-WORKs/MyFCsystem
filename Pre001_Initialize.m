clear all
close all

% Physical property
R = 287;    % [J/kg/K]
kap = 1.4;  % [-]

% Initial Condition
T_init = 293.15;    % [K]
P_init = 101325;    % [Pa]

% Boundary Condition
T_amb = 293.15;    % [K]
P_amb = 101325;    % [Pa]

% Set Compressor Parameter
cmp.param.dc = 0.2286;      % [m]
cmp.param.eta_is = 0.8;     % [-]
cmp.param.eta_mech = 0.98;  % [-]
cmp.map.coeff.a = [-3.69906e-5, 2.70399e-4, -5.36235e-4, -4.63685e-5, 2.21195e-3];
cmp.map.coeff.b = [1.766467,-1.34837,2.44419];
cmp.map.coeff.c = [-9.78755e-3, 0.10581, -0.42937, 0.80121, -0.68344, 0.43331];