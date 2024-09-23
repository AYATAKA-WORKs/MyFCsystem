clear all
%% M code check
N = 50000;          % [rpm]
omega = N*pi/30;    % [rad/s]
E = 170;            % [V]

kv = 0.0153;        % [V/(rad/s)]
kt = 0.0153;        % [Nm/A]
R_M = 0.82;         % [Ω]
eta_M = 0.98;       % [-]

tau = eta_M * kt/R_M *(E - kv*omega);   % [Nm]


