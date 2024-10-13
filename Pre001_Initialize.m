clear all
close all

addpath("SLX\");
addpath("TEST\");

% Physical property
R = 287;    % [J/kg/K]
kap = 1.4;  % [-]

% Dry air propety
% Mass fraction
wf_O2_DA = 0.2;
wf_N2_DA = 0.8;
% Molecular weight
M_O2 = 32;
M_N2 = 28;
M_DA = wf_O2_DA*M_O2 + wf_N2_DA*M_N2;
% Mol fraction
molf_O2_DA = (wf_O2_DA*M_DA/M_O2) / ((wf_O2_DA*M_DA/M_O2) + (wf_N2_DA*M_DA/M_N2));
molf_N2_DA = (wf_N2_DA*M_DA/M_N2) / ((wf_O2_DA*M_DA/M_O2) + (wf_N2_DA*M_DA/M_N2));

% Initial Condition
T_init = 293.15;    % [K]
p_init = 101325;    % [Pa]

% Boundary Condition
T_amb = 293.15;    % [K]
p_amb = 101325;    % [Pa]
RH_amb = 0.5;      % [-]

% Set Compressor Parameter
cmp.param.dc = 0.2286;      % [m]
cmp.param.eta_is = 0.8;     % [-]
cmp.param.eta_mech = 0.98;  % [-]
cmp.map.coeff.a = [-3.69906e-5, 2.70399e-4, -5.36235e-4, -4.63685e-5, 2.21195e-3];
cmp.map.coeff.b = [1.766467,-1.34837,2.44419];
cmp.map.coeff.c = [-9.78755e-3, 0.10581, -0.42937, 0.80121, -0.68344, 0.43331];