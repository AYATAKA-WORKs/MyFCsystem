clear all
close all

%% パラメータの初期設定
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
cmp.ref.p = 101325;         % [Pa]
cmp.ref.T = 293.15;         % [K]
cmp.ref.R = 287;            % [J/kg/K]
cmp.ref.kap = 1.4;          % [-]
cmp.ref.cp = 1005;          % [J/kg/K]

% Physical property
R = cmp.ref.R;      % [J/kg/K]
kap = cmp.ref.kap;  % [-]
cp = cmp.ref.cp;    % [J/kg/K]

%% コンプレッサマップの作成
N_table = [10000:10000:100000];
pr = [1:0.1:4];
p_in = P_amb;
T_in = T_amb;
T_ref = cmp.ref.T;
R_ref = cmp.ref.R;
p_ref = cmp.ref.p;

Wcr = zeros(size(N_table,2),size(pr,2));
T_out = zeros(size(N_table,2),size(pr,2));
torque = zeros(size(N_table,2),size(pr,2));

for i = 1:size(N_table,2)
    N = N_table(i);
    for j= 1:size(pr,2)
        p_out = p_in *pr(j);
        sim(".\SLX\Simple_Compressor.slx");
        Wcr(i,j) = ans.W_cmp(end);
        T_out(i,j) = ans.T_out(end);
        torque(i,j) = ans.torque(end);
    end
end

%% Compressor P-Q Map
figure('Name','Compressor P-Q Map')
hold on
grid on
box on

for i = 1:size(N_table,2)
    plot(Wcr(i,:),pr);
end

xlim([0 0.1]);
xlabel("Corrected mass flow rate [kg/s]");
ylabel("Pressure ratio");

%% 出口温度 vs 圧力比
figure('Name','Outlet Temperature vs Pressure Ratio')
hold on
grid on
box on

for i = 1:size(N_table,2)
    plot(T_out(i,:),pr);
end
xlabel("Temrature [K]");
ylabel("Pressure [Pa]");
%% torque vs 圧力比
figure('Name','Torque vs Pressure Ratio')
hold on
grid on
box on

for i = 1:size(N_table,2)
    plot(torque,pr);
end
xlabel("Temrature [K]");
ylabel("Pressure [Pa]");