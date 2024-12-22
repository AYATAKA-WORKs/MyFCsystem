clear rslt
%% Get simulation results
% Get signal lists
sigList = out.logsout.getElementNames;

% Ambient
rslt.p_amb = out.logsout.getElement("p_amb");
rslt.T_amb = out.logsout.getElement("T_amb");
rslt.RH_amb = out.logsout.getElement("RH_amb");
rslt.yO2_amb = out.logsout.getElement("yO2_amb");

% Compressor
rslt.Win_cmp = out.logsout.getElement("W_cmp");
rslt.pin_cmp = out.logsout.getElement("p_amb");
rslt.Tin_cmp = out.logsout.getElement("T_amb");
rslt.Wout_cmp = out.logsout.getElement("W_cmp");
rslt.pout_cmp = out.logsout.getElement("pI_sm");
rslt.Tout_cmp = out.logsout.getElement("Tout_cmp");
rslt.N_etc = out.logsout.getElement("N_etc");

% FC stack
rslt.Win_st = out.logsout.getElement("Win_st");
rslt.pin_st = out.logsout.getElement("pin_st");
rslt.Tin_st = out.logsout.getElement("Tin_st");
rslt.RHin_st = out.logsout.getElement("RHin_st");
rslt.yO2in_st = out.logsout.getElement("yO2in_st");
rslt.Wout_st = out.logsout.getElement("Wout_st");
rslt.pout_st = out.logsout.getElement("pin_tbn");
rslt.Tout_st = out.logsout.getElement("Tout_st");
rslt.RHout_st = out.logsout.getElement("RHout_st");
rslt.yO2out_st = out.logsout.getElement("yO2out_st");

% Turbine
rslt.Win_tbn = out.logsout.getElement("Win_tbn");
rslt.pin_tbn = out.logsout.getElement("pin_tbn");
rslt.Tin_tbn = out.logsout.getElement("Tin_tbn");
rslt.Wout_tbn = out.logsout.getElement("Win_tbn");
rslt.pout_tbn = out.logsout.getElement("pout_tbn");
rslt.Tout_tbn = out.logsout.getElement("Tout_tbn");

% Power
rslt.P_cmp = out.logsout.getElement("P_cmp");
rslt.P_tbn = out.logsout.getElement("P_tbn");
rslt.P_M = out.logsout.getElement("P_M");

% Efficiency
rslt.eta_is_cmp = out.logsout.getElement("eta_is_cmp");
rslt.eta_is_tbn = out.logsout.getElement("eta_is_tbn");

%% Make table of final values
rslt.table.amb = ...
[...
    rslt.p_amb.Values.Data(end);
    rslt.T_amb.Values.Data(end);
    rslt.RH_amb.Values.Data(end);
    rslt.yO2_amb.Values.Data(end); ...
];
rslt.table.cmp = ...
[...
    rslt.Win_cmp.Values.Data(end);
    rslt.pin_cmp.Values.Data(end);
    rslt.Tin_cmp.Values.Data(end);
    rslt.Wout_cmp.Values.Data(end);
    rslt.pout_cmp.Values.Data(end);
    rslt.Tout_cmp.Values.Data(end);
    rslt.N_etc.Values.Data(end); ...
];
rslt.table.FCstack = ...
[...
    rslt.Win_st.Values.Data(end);
    rslt.pin_st.Values.Data(end);
    rslt.Tin_st.Values.Data(end);
    rslt.RHin_st.Values.Data(end);
    rslt.yO2in_st.Values.Data(end);
    rslt.Wout_st.Values.Data(end);
    rslt.pout_st.Values.Data(end);
    rslt.Tout_st.Values.Data(end);
    rslt.RHout_st.Values.Data(end);
    rslt.yO2out_st.Values.Data(end); ...
];
rslt.table.tbn = ...
[...
    rslt.Win_tbn.Values.Data(end);
    rslt.pin_tbn.Values.Data(end);
    rslt.Tin_tbn.Values.Data(end);
    rslt.Wout_tbn.Values.Data(end);
    rslt.pout_tbn.Values.Data(end);
    rslt.Tout_tbn.Values.Data(end); ...
];
rslt.table.power = ...
[...
    rslt.P_cmp.Values.Data(end);
    rslt.P_tbn.Values.Data(end);
    rslt.P_M.Values.Data(end); ...
];
rslt.table.efficiency = ...
[...
    rslt.eta_is_cmp.Values.Data(end);
    rslt.eta_is_tbn.Values.Data(end);
];

