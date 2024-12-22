classdef GasProperties
    properties
        p,      T,      RH,     yO2
        M_DA,   R_DA,   cp_DA,  kap_DA
        M_MA,   R_MA,   cp_MA,  kap_MA,
        M_Nmix, R_Nmix, cp_Nmix,kap_Nmix
    end
    properties (Constant = true)
        % 物理定数の定義
        R0 = 8.31446261815324;  % 気体定数 [J/K/mol]
        h  = 6.62607015e-34;    % プランク定数 [Js]
        kB = 1.380649e-23;      % ボルツマン定数 [J/K]
        c = 2.99792458e10;      % 光速 [cm/s]
        % 気体原子量
        M_H = 1.00794e-3;
        M_O = 15.9994e-3;
        M_N = 14.00674e-3;
        M_C = 12.011e-3;
        M_Ar = 39.948e-3;
        % 気体分子の振動数の定義
        mu_N2 = 2359;       % [s/cm]
        mu_O2 = 1580;       % [s/cm]
        % 二酸化炭素が持つ振動成分
        % 1 対称伸縮振動，2 逆対称伸縮振動，3 変角振動（紙面内と紙面垂直で2種類）
        mu_CO2_1 = 1383;    % [s/cm]
        mu_CO2_2 = 2349;    % [s/cm]
        mu_CO2_3 = 667 ;    % [s/cm]
        % 水蒸気が持つ振動成分
        % 1 対称伸縮振動，2 逆対称伸縮振動，3 変角振動（軸非対称なので1種類）
        mu_H2O_1 = 3832;    % [s/cm]
        mu_H2O_2 = 3943;    % [s/cm]
        mu_H2O_3 = 1649;    % [s/cm]
    end
    methods
        function [M_DA,R_DA,cp_DA,kap_DA] = getDAprop(obj,T,yO2)
            % 乾燥空気の物性値を返す関数
            obj.T =T;
            obj.yO2 =yO2;
            obj = obj.UpdateProperties_DA;
            M_DA = obj.M_DA;
            R_DA = obj.R_DA;
            cp_DA = obj.cp_DA;
            kap_DA = obj.kap_DA;
        end
        function [M_MA,R_MA,cp_MA,kap_MA] = getMAprop(obj,p,T,RH,yO2)
            % 湿り空気の物性値を返す関数
            obj.p =p;
            obj.T =T;
            obj.RH =RH;
            obj.yO2 =yO2;
            obj = obj.UpdateProperties_MA;
            M_MA = obj.M_MA;
            R_MA = obj.R_MA;
            cp_MA = obj.cp_MA;
            kap_MA = obj.kap_MA;
        end
        function obj = UpdateProperties(obj,p,T,RH,yO2)
            % 湿り空気の物性値を入力に応じて更新する関数
            obj.p =p;
            obj.T =T;
            obj.RH =RH;
            obj.yO2 =yO2;
            obj = obj.UpdateProperties_MA;
        end
    end
    methods (Access = protected)
        function obj = UpdateProperties_DA(obj)
            % 乾燥空気の物性値更新計算
            % 気体分子量
            M_O2 = obj.M_O*2;
            M_N2 = obj.M_N*2;
            M_CO2 = obj.M_C + obj.M_O*2;
            
            % 組成比（体積比 = 大気圧下ではモル比に等しい）
            nf_O2 = obj.yO2;
            nf_Ar = 0.934e-2;
            nf_CO2 = 0.0314e-2;
            nf_N2 = 1 - (nf_O2 + nf_Ar + nf_CO2);
            
            % 混合気体の分子量
            obj.M_Nmix = (nf_N2*M_N2 + nf_Ar*obj.M_Ar + nf_CO2*M_CO2)/(nf_N2+nf_Ar+nf_CO2);
            obj.M_DA   = nf_N2*M_N2 + nf_O2*M_O2 + nf_Ar*obj.M_Ar + nf_CO2*M_CO2;
                        
            % 乾燥空気のガス定数
            obj.R_Nmix = obj.R0/obj.M_Nmix;
            obj.R_DA = obj.R0/obj.M_DA;
                        
            % 各分子の定圧比熱の計算
            % 理想気体においてはモル比熱は並進・回転・振動の自由度で決まる
            % 温度依存性は量子力学によって導かれる振動に影響する
            % 気体の定圧モル比熱計算
            cp0_Ar  = 5/2*obj.R0;   % [J/K/mol] 20.786
            cp0_O2  = 7/2*obj.R0 + obj.calc_cv_m(obj.T,obj.mu_O2*obj.c,obj.R0,obj.h,obj.kB);  % [J/K/mol] 29.335
            cp0_N2  = 7/2*obj.R0 + obj.calc_cv_m(obj.T,obj.mu_N2*obj.c,obj.R0,obj.h,obj.kB);  % [J/K/mol] 29.124
            cp0_CO2 = 7/2*obj.R0 + obj.calc_cv_m(obj.T,obj.mu_CO2_1*obj.c,obj.R0,obj.h,obj.kB) + obj.calc_cv_m(obj.T,obj.mu_CO2_2*obj.c,obj.R0,obj.h,obj.kB) + 2*obj.calc_cv_m(obj.T,obj.mu_CO2_3*obj.c,obj.R0,obj.h,obj.kB);   % [J/K/mol] 37.14
            % cp0_vp  = 4*obj.R0   + obj.calc_cv_m(obj.T,obj.mu_H2O_1*obj.c,obj.R0,obj.h,obj.kB) + obj.calc_cv_m(obj.T,obj.mu_H2O_2*obj.c,obj.R0,obj.h,obj.kB) + obj.calc_cv_m(obj.T,obj.mu_H2O_3*obj.c,obj.R0,obj.h,obj.kB);     % [J/K/mol] 33.577
            cp0_Nmix = (nf_N2*cp0_N2 + nf_Ar*cp0_Ar + nf_CO2*cp0_CO2) / (nf_N2+nf_Ar+nf_CO2);
            cp0_DA   = nf_N2*cp0_N2 + nf_O2*cp0_O2 + nf_Ar*cp0_Ar + nf_CO2*cp0_CO2;
            % 気体の定圧比熱計算
            % cp_Ar  = cp0_Ar / obj.M_Ar; % [J/K/kg] 
            % cp_O2  = cp0_O2 / M_O2;     % [J/K/kg]  
            % cp_N2  = cp0_N2 / M_N2;     % [J/K/kg] 
            % cp_CO2 = cp0_CO2 / M_CO2;   % [J/K/kg] 
            obj.cp_Nmix = cp0_Nmix / obj.M_Nmix;
            obj.cp_DA = cp0_DA / obj.M_DA;
            % 参考文献
            % 気体の熱容量と分子量  : https://www.jstage.jst.go.jp/article/kakyoshi/71/1/71_24/_pdf/-char/ja
            % 物理化学 講義資料     : https://online.lec-jp.com/images/goods_book/KL/KL10058.pdf
            % 放送大学 量子化学     : https://info.ouj.ac.jp/~hamada/Quantumch/subject/cq/chap10/figure/cq98af06.html
            
            % 比熱比の計算
            obj.kap_Nmix = obj.cp_Nmix/(obj.cp_Nmix-obj.R_Nmix);
            obj.kap_DA = obj.cp_DA/(obj.cp_DA-obj.R_DA);
        end

        function obj = UpdateProperties_MA(obj)
            % 湿り空気の物性値更新計算
            % 乾燥空気の物性値を更新してから湿り空気の物性値計算をする

            % 乾燥空気の物性値計算
            obj = UpdateProperties_DA(obj);
            
            % 水蒸気分子量
            M_H2O = obj.M_H*2 + obj.M_O;

            % 組成比（圧力比 = モル比）
            p_vp = obj.RH*obj.sat_vp_pressure(obj.T);
            p_DA = obj.p - p_vp;
            nf_vp = p_vp/obj.p;
            nf_DA = p_DA/obj.p;
            
            % 湿り空気の分子量
            obj.M_MA = (nf_DA*obj.M_DA + nf_vp*M_H2O)/(nf_DA+nf_vp); % [kg/mol]
            
            % 湿り空気空気のガス定数
            obj.R_MA = obj.R0/obj.M_MA; % [J/K/kg]
            
            % 各分子の定圧比熱の計算
            % 理想気体においてはモル比熱は並進・回転・振動の自由度で決まる
            % 温度依存性は量子力学によって導かれる振動に影響する
            % 気体の定圧モル比熱計算
            cp0_vp  = 4*obj.R0   + obj.calc_cv_m(obj.T,obj.mu_H2O_1*obj.c,obj.R0,obj.h,obj.kB) + obj.calc_cv_m(obj.T,obj.mu_H2O_2*obj.c,obj.R0,obj.h,obj.kB) + obj.calc_cv_m(obj.T,obj.mu_H2O_3*obj.c,obj.R0,obj.h,obj.kB);  % [J/K/mol] 33.577
            cp0_MA  = nf_DA*obj.cp_DA*obj.M_DA + nf_vp*cp0_vp;
            % 気体の定圧比熱計算
            obj.cp_MA = cp0_MA / obj.M_MA;  % [J/K/kg]
            % 参考文献
            % 気体の熱容量と分子量  : https://www.jstage.jst.go.jp/article/kakyoshi/71/1/71_24/_pdf/-char/ja
            % 物理化学 講義資料     : https://online.lec-jp.com/images/goods_book/KL/KL10058.pdf
            % 放送大学 量子化学     : https://info.ouj.ac.jp/~hamada/Quantumch/subject/cq/chap10/figure/cq98af06.html
            
            % 比熱比の計算
            obj.kap_MA = obj.cp_MA/(obj.cp_MA-obj.R_MA);    % [J/K/kg]
        end

        % 気体分子の振動数による熱容量[J/K/mol]の計算
        function cv_m = calc_cv_m(obj,T,mu,R0,h,kB)
            Theta = h*mu/kB;
            cv_m = R0 * (Theta/T)^2 * exp(-Theta/T)/(1-exp(-Theta/T))^2;
        end

        % Calculate saturated vapor pressure [Pa]
        function psat = sat_vp_pressure(obj,Ts)
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
    end

end

