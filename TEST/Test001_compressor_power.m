clear all

%% Compressor Parameter
cmp.param.dc = 0.2286;      % [m]
cmp.param.eta_is = 0.8;     % [-]
cmp.param.eta_mech = 0.98;  % [-]
cmp.map.coeff.a = [-3.69906e-5, 2.70399e-4, -5.36235e-4, -4.63685e-5, 2.21195e-3];
cmp.map.coeff.b = [1.766467,-1.34837,2.44419];
cmp.map.coeff.c = [-9.78755e-3, 0.10581, -0.42937, 0.80121, -0.68344, 0.43331];

%% Air Property
R = 287;            % [J/kg/K]
kap = 1.4;          % [-]
cp = 1005;          % [J/kg/K]

%% Operational Condition
N = 90000;          % [rpm]
omega = N*pi/30;    % [rad/s]

T_in = 293.15;      % [K]
p_in = 101325;      % [Pa]
p_out = 101325*2.5; % [Pa]

%% Simulink Model check
out = sim("..\SLX\Simple_Compressor.slx");
tau_slx = out.torque;
W_in_slx = out.W_cmp;

%% M code check
pr = p_out/p_in;
W_in_m = W_in_slx;         % [kg/s]
tau_m = cp*T_in/omega/cmp.param.eta_is/cmp.param.eta_mech*(pr^((kap-1)/kap)-1)*W_in_m;   % [Nm]
