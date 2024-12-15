clear all

%% 計算条件
% DRY側条件
W_in_dry = 0.1;
p_in_dry = 101325;
T_in_dry = 373.15;
RH_in_dry = 0.1;
yO2_in_dry = 0.2;

% WET側条件
W_in_wet = 0.1;
p_in_wet = 101325;
T_in_wet = 373.15;
RH_in_wet = 0.9;
yO2_in_wet = 0.2;

RH_out_dry_tar = 0.9;

% 水注入量の計算
 W_vp_inj_ctrl = calc_W_vp_inj_ctrl(W_in_dry, p_in_dry, T_in_dry, RH_in_dry, yO2_in_dry,W_in_wet, p_in_wet, T_in_wet, RH_in_wet, yO2_in_wet,RH_out_dry_tar);

%% 湿り空気計算開始
W_vp_inj_lim = calc_W_vp_inj_lim(W_in_dry, p_in_dry, T_in_dry, RH_in_dry, yO2_in_dry,W_in_wet, p_in_wet, T_in_wet, RH_in_wet, yO2_in_wet);
W_vp_inj = min(W_vp_inj_ctrl, W_vp_inj_lim);

[W_out_dry, p_out_dry, T_out_dry, RH_out_dry, yO2_out_dry] = HM_std(W_in_dry, p_in_dry, T_in_dry, RH_in_dry, yO2_in_dry, W_vp_inj);
[W_out_wet, p_out_wet, T_out_wet, RH_out_wet, yO2_out_wet] = HM_std(W_in_wet, p_in_wet, T_in_wet, RH_in_wet, yO2_in_wet, -W_vp_inj);

%% 加湿／除湿計算
function [W_out, p_out, T_out, RH_out, yO2_out] = HM_std(W_in, p_in, T_in, RH_in, yO2_in, W_vp_inj)
    % 湿り空気の状態計算
    [W_O2_in, W_N2_in, W_vp_in, p_O2_in, p_N2_in, p_vp_in, M_in, R_in, cp_in, kap_in] = calc_MA_property(W_in, p_in, T_in, RH_in, yO2_in);

    % 水蒸気の注入
    W_vp_out = W_vp_in + W_vp_inj;  % [kg/s]

    % 出口流量計算
    W_O2_out = W_O2_in; % [kg/s]
    W_N2_out = W_N2_in; % [kg/s]
    W_out = W_O2_out + W_N2_out + W_vp_out; % [kg/s]

    % 出口圧計算
    p_O2_out = p_O2_in; % [Pa]
    p_N2_out = p_N2_in; % [Pa]
    if W_vp_out <= 0
        p_vp_out = 0;
    else
        p_vp_out = W_vp_out/W_vp_in * p_vp_in;  % [Pa]
    end
    p_out = p_O2_out + p_N2_out + p_vp_out;     % [Pa]

    % 温度は一定とする
    T_out = T_in;       % [K]

    % 相対湿度計算
    RH_out = p_vp_out/sat_vp_pressure(T_out);

    % 空気組成は一定とする
    yO2_out = yO2_in;   % [-]
end

%% 水交換量のリミット計算
function W_vp_inj_lim = calc_W_vp_inj_lim(W_in_dry, p_in_dry, T_in_dry, RH_in_dry, yO2_in_dry,W_in_wet, p_in_wet, T_in_wet, RH_in_wet, yO2_in_wet)
    % 物理定数の定義
    % 分子量定義
    % M_O2 = 32e-3;   % [kg/mol]
    % M_N2 = 28e-3;   % [kg/mol]
    M_vp = 18e-3;   % [kg/mol]

    % 湿り空気の分圧計算
    % DRY側
    [W_O2_in_dry, W_N2_in_dry, W_vp_in_dry, p_O2_in_dry, p_N2_in_dry, p_vp_in_dry, M_in_dry, R_in_dry, cp_in_dry, kap_in_dry] = calc_MA_property(W_in_dry, p_in_dry, T_in_dry, RH_in_dry, yO2_in_dry);
    % WET側
    [W_O2_in_wet, W_N2_in_wet, W_vp_in_wet, p_O2_in_wet, p_N2_in_wet, p_vp_in_wet, M_in_wet, R_in_wet, cp_in_wet, kap_in_wet] = calc_MA_property(W_in_wet, p_in_wet, T_in_wet, RH_in_wet, yO2_in_wet);

    % DRY側リミット
    % RH=1の時の水蒸気流量
    W_vp_RH1 = sat_vp_pressure(T_in_dry)/p_in_dry * M_vp/M_in_dry * W_in_dry;       % 湿り空気ベース計算
    % W_vp_RH1 = sat_vp_pressure(T_in_dry)/p_O2_in_dry * M_vp/M_O2 * W_O2_in_dry;   % 酸素ベース計算
    % W_vp_RH1 = sat_vp_pressure(T_in_dry)/(p_O2_in_dry+p_N2_in_dry) * M_vp/(yO2_in_dry*M_O2+(1-yO2_in_dry)*M_N2) * (W_O2_in_dry+W_N2_in_dry); % 乾燥空気ベース計算
    % リミット交換量
    W_vp_inj_lim_dry =  W_vp_RH1 - W_vp_in_dry;

    % WET側リミット
    W_vp_inj_lim_wet =  W_vp_in_wet - 0;

    % 小さいほうがリミットとなる
    W_vp_inj_lim = min(W_vp_inj_lim_dry, W_vp_inj_lim_wet);
end

%% 水交換制御量の計算
function W_vp_inj_ctrl = calc_W_vp_inj_ctrl(W_in_dry, p_in_dry, T_in_dry, RH_in_dry, yO2_in_dry,W_in_wet, p_in_wet, T_in_wet, RH_in_wet, yO2_in_wet,RH_out_dry_tar)
    % 物理定数の定義
    % 分子量定義
    % M_O2 = 32e-3;   % [kg/mol]
    % M_N2 = 28e-3;   % [kg/mol]
    M_vp = 18e-3;   % [kg/mol]

    % 湿り空気の分圧計算
    % DRY側
    [W_O2_in_dry, W_N2_in_dry, W_vp_in_dry, p_O2_in_dry, p_N2_in_dry, p_vp_in_dry, M_in_dry, R_in_dry, cp_in_dry, kap_in_dry] = calc_MA_property(W_in_dry, p_in_dry, T_in_dry, RH_in_dry, yO2_in_dry);
    % WET側
    [W_O2_in_wet, W_N2_in_wet, W_vp_in_wet, p_O2_in_wet, p_N2_in_wet, p_vp_in_wet, M_in_wet, R_in_wet, cp_in_wet, kap_in_wet] = calc_MA_property(W_in_wet, p_in_wet, T_in_wet, RH_in_wet, yO2_in_wet);

    % DRY側リミット
    % RH=1の時の水蒸気流量
    p_vp_tar_dry = RH_out_dry_tar * sat_vp_pressure(T_in_dry);
    W_vp_RH_tar = p_vp_tar_dry/p_in_dry * M_vp/M_in_dry * W_in_dry;       % 湿り空気ベース計算
    % W_vp_RH1 = sat_vp_pressure(T_in_dry)/p_O2_in_dry * M_vp/M_O2 * W_O2_in_dry;   % 酸素ベース計算
    % W_vp_RH1 = sat_vp_pressure(T_in_dry)/(p_O2_in_dry+p_N2_in_dry) * M_vp/(yO2_in_dry*M_O2+(1-yO2_in_dry)*M_N2) * (W_O2_in_dry+W_N2_in_dry); % 乾燥空気ベース計算
    % 制御水交換量
    W_vp_inj_ctrl_dry =  W_vp_RH_tar - W_vp_in_dry;

    % WET側リミット
    W_vp_inj_lim_wet =  W_vp_in_wet - 0;

    % 小さいほうが制御量となる
    W_vp_inj_ctrl = min(W_vp_inj_ctrl_dry, W_vp_inj_lim_wet);
end

%% 湿り空気の状態量計算
function [W_O2, W_N2, W_vp, p_O2, p_N2, p_vp, M, R, cp, kap] = calc_MA_property(W, p, T, RH, yO2)
    % 物理定数の定義
    % 分子量定義
    M_O2 = 32e-3;   % [kg/mol]
    M_N2 = 28e-3;   % [kg/mol]
    M_vp = 18e-3;   % [kg/mol]
    % 気体定数定義
    R0 = 8.31446261815324; % [J/K/mol]
    
    % 分圧計算
    p_vp_sat = sat_vp_pressure(T);
    p_vp = RH * p_vp_sat;
    p_DA = p - p_vp;
    p_O2 = yO2 * p_DA;
    p_N2 = p_DA - p_O2;
    
    % 湿り空気のモル分率計算
    x_O2 = p_O2/p;
    x_N2 = p_N2/p;
    x_vp = p_vp/p;

    % 湿り空気の分子量計算
    M = M_O2*x_O2 + M_N2*x_N2 + M_vp*x_vp;       % [kg/mol]

    % ガス定数計算
    R = R0/M;   % [J/K/kg]

    % 質量分率計算
    w_O2 = x_O2*M_O2 / M;
    w_N2 = x_N2*M_N2 / M;
    w_vp = x_vp*M_vp / M;

    % 各分子の定圧比熱の計算
    % 二原子分子
    cp_O2 = 7/2*R0 / M_O2;  % [J/K/kg]
    cp_N2 = 7/2*R0 / M_N2;  % [J/K/kg]
    % 水蒸気
    cp_vp = 4*R0 / M_vp;    % [J/K/kg]

    % 混合気体の定圧比熱
    cp = cp_O2*w_O2 + cp_N2*w_N2 + cp_vp*w_vp;  % [J/K/kg]

    % 混合気体の比熱比
    kap = cp/(cp-R);

    % 各分子の質量流量
    W_O2 = W * w_O2;    % [kg/s]
    W_N2 = W * w_N2;    % [kg/s]
    W_vp = W * w_vp;    % [kg/s]
end

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