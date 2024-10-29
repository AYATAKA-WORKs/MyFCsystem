clear all
%close all

%% パラメータの初期設定
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
N_table = linspace(10000,100000,10);
pr = linspace(1,4,100);
T_ref = cmp.ref.T;
R_ref = cmp.ref.R;
p_ref = cmp.ref.p;

Wcr = zeros(size(N_table,2),size(pr,2));

Uc = N_table * pi/30 * (cmp.param.dc/2);
MUc = Uc / sqrt(cmp.ref.kap*cmp.ref.R*cmp.ref.T);
a = cmp.map.coeff.a;
b = cmp.map.coeff.b;
c = cmp.map.coeff.c;
phi_max = a(1)*MUc.^4 + a(2)*MUc.^3 + a(3)*MUc.^2 + a(4)*MUc + a(5);
beta = b(1)*MUc.^2 + b(2)*MUc + b(3);
MUy_max = c(1)*MUc.^5 + c(2)*MUc.^4 + c(3)*MUc.^3 + c(4)*MUc.^2 + c(5)*MUc + c(6);

MUy = zeros(size(Uc,2),size(pr,2));
phi = zeros(size(Uc,2),size(pr,2));
for i = 1:size(Uc,2)
    MUy(i,:) = cp*T_ref*(pr.^((kap-1)/kap)-1) / (Uc(i)^2/2);
    phi(i,:) = phi_max(i) * (1 -  exp((MUy(i,:)/MUy_max(i) -1)*beta(i)));
    Wcr(i,:) = phi(i,:) * Uc(i) * p_ref/(R_ref*T_ref) * pi/4*cmp.param.dc^2;
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

figure('Name','Compressor non-dimensional Map')
hold on
grid on
box on

for i = 1:size(N_table,2)
    plot(phi(i,:),MUy(i,:));
end

xlim([0 max(phi,[],"all")]);
ylim([0 0.4])
xlabel("Flow coefficient [-]");
ylabel("Pressure coefficient [-]");

