% フォルダが存在しない場合、作成する
outputFolder = 'FIG';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
Sim002_S000_make_report

%% すべての結果を表示
f = figure; 
set(f, 'Name', 'All results', 'NumberTitle', 'on');

%% figureの設定
% フィギュアのサイズを指定
figWidth = 600;  % 幅
figHeight = 1000; % 高さ

fontSizeLabel = 8;  % 軸ラベルのフォントサイズ
fontSizeTicks = 9;  % 軸目盛りのフォントサイズ
lineWidth = 1.5;

% 各カテゴリのy軸スケール設定
yScale = struct(... 
    'Flow_Rate', [0 0.1], ...  % 空白をアンダースコアに変更
    'Pressure', [0, 300000], ...
    'Temperature', [250, 400], ...
    'Humidity', [0, 1], ...
    'Speed', [0, 100000], ...
    'Voltage', [0, 300], ...
    'VGS_Opening', [-2e-4 2e-4], ...  % ここも変更
    'Power', [0 10000] ...
);

colors = containers.Map(...
    { 'Win_st', 
      'W_stack_cmd', 
      'pin_st', 
      'p_stack_cmd', 
      'Tin_st', 
      'RHin_st', 
      'RHout_st', 
      'N_etc', 
      'N_cmd', 
      'v_M', 
      'VGS_cmd', 
      'P_cmp', 
      'P_tbn', 
      'P_M' }, ...  % 変数名のキー
    { [0.2, 0.4, 0.6], ...   % 明るい青
      [0, 0, 0], ...         % cmd系は黒
      [0.4, 0.6, 0.4], ...   % グリーン
      [0, 0, 0], ...         % cmd系は黒
      [0.5, 0.3, 0.4], ...   % ワインレッド
      [0.9, 0.7, 0.3], ...   % マスタードイエロー
      [1.0, 0.6, 0.2], ...   % マスタードオレンジ
      [0.2, 0.4, 0.6], ...   % 明るい青
      [0, 0, 0], ...         % cmd系は黒
      [0.7, 0.2, 0.3], ...   % 深い赤
      [0.3, 0.5, 0.2], ...   % オリーブグリーン
      [0.5, 0.5, 0.3], ...   % オリーブ
      [0.4, 0.2, 0.6], ...   % 紫
      [0.7, 0.5, 0.2] ...    % アーストーン
      }...
  );   

%% plot処理
categories = {'Flow Rate', 'Pressure', 'Temperature', 'Humidity', 'Speed', 'Voltage', 'VGS Opening', 'Power'};
variables = {... 
    {'Flow Rate', 'Win_st', 'スタック入口流量', 'kg/s'}, ...
    {'Flow Rate', 'W_stack_cmd', 'スタック流量目標', 'kg/s'}, ...
    {'Pressure', 'pin_st', 'スタック入口圧力', 'Pa'}, ...
    {'Pressure', 'p_stack_cmd', 'スタック圧力目標', 'Pa'}, ...
    {'Temperature', 'Tin_st', 'スタック入口温度', 'K'}, ...
    {'Humidity', 'RHin_st', 'スタック入口湿度', 'RH'}, ...
    {'Humidity', 'RHout_st', 'スタック出口湿度', 'RH'}, ...
    {'Speed', 'N_etc', '実回転数', 'rpm'}, ...
    {'Speed', 'N_cmd', '指令回転数', 'rpm'}, ...
    {'Voltage', 'v_M', 'モータ電圧', 'V'}, ...
    {'VGS Opening', 'VGS_cmd', 'VGS指令開度', 'm^2'}, ...
    {'Power', 'P_cmp', 'コンプレッサ動力', 'W'}, ...
    {'Power', 'P_tbn', 'タービン動力', 'W'}, ...
    {'Power', 'P_M', 'モータ動力', 'W'}
};

numCategories = length(categories);
% Subplot の表示
for c = 1:numCategories
    category = categories{c};
    idx = find(strcmp(category, cellfun(@(x) x{1}, variables, 'UniformOutput', false)));
    numPlots = length(idx);
    
    subplot(numCategories, 1, c);
    hold on;
    for i = 1:numPlots
        varName = variables{idx(i)}{2};
        signalName = variables{idx(i)}{3};
        unit = variables{idx(i)}{4};
        
        if contains(signalName, '目標') || contains(signalName, '指令')
            plot(rslt.(varName).Values.Time, rslt.(varName).Values.Data, '--k', 'LineWidth', lineWidth, 'DisplayName', signalName);
        else
            plot(rslt.(varName).Values.Time, rslt.(varName).Values.Data, 'Color', colors(varName), 'LineWidth', lineWidth, 'DisplayName', signalName);
        end
    end
    grid on;
    box on;
    
    % 最後のグラフ以外で横軸のラベルを表示しない
    if c == numCategories
        xlabel('時間 [s]', 'FontSize', fontSizeLabel);  % 最後だけ表示
    end
    
    % 軸ラベルのフォントサイズを設定
    ylabel([category ' [' variables{idx(1)}{4} ']'], 'FontSize', fontSizeLabel);
    
    % タイトルを表示しない
    % title(category, 'FontSize', fontSizeLabel);  % この行を削除
    
    % 軸目盛りのフォントサイズを設定
    set(gca, 'FontSize', fontSizeTicks);
    
    % y軸のスケールを指定
    if strcmp(yScale.(strrep(category, ' ', '_')), 'auto')
        ylim auto;  % 'auto' の場合、自動スケーリング
    else
        ylim(yScale.(strrep(category, ' ', '_')));  % ユーザ指定の範囲
    end
    
    % レジェンドを右下に配置
    legend('FontSize', fontSizeLabel, 'Location', 'southeast');
    hold off;
end

% フィギュアを画面中央に配置
screenSize = get(0, 'ScreenSize');      % 画面サイズを取得
set(gcf, 'Position', [... 
    (screenSize(3) - figWidth) / 2, ...  % 左からの距離 (画面中央に配置)
    (screenSize(4) - figHeight) / 2, ... % 下からの距離 (画面中央に配置)
    figWidth, figHeight ...              % 幅と高さ
]);

% 画像として保存
saveas(f, fullfile(outputFolder, '0_All_results.png'));
% FIG形式で保存
savefig(f, fullfile(outputFolder, '0_All_results.fig'));

%% 各plotの表示
for c = 1:numCategories
    % 各カテゴリごとに新しいfigureを作成
    f = figure;

    % ここでfigureの名前を設定
    set(f, 'Name', categories{c}, 'NumberTitle', 'off');

    category = categories{c};
    idx = find(strcmp(category, cellfun(@(x) x{1}, variables, 'UniformOutput', false)));
    numPlots = length(idx);

    hold on;
    for i = 1:numPlots
        varName = variables{idx(i)}{2};
        signalName = variables{idx(i)}{3};
        unit = variables{idx(i)}{4};

        if contains(signalName, '目標') || contains(signalName, '指令')
            plot(rslt.(varName).Values.Time, rslt.(varName).Values.Data, '--k', 'LineWidth', lineWidth, 'DisplayName', signalName);
        else
            plot(rslt.(varName).Values.Time, rslt.(varName).Values.Data, 'Color', colors(varName), 'LineWidth', lineWidth, 'DisplayName', signalName);
        end
    end
    grid on;
    box on;

    xlabel('時間 [s]', 'FontSize', fontSizeLabel);

    % 軸ラベルのフォントサイズを設定
    ylabel([category ' [' variables{idx(1)}{4} ']'], 'FontSize', fontSizeLabel);

    % タイトルを表示しない
    % title(category, 'FontSize', fontSizeLabel);  % この行を削除

    % 軸目盛りのフォントサイズを設定
    set(gca, 'FontSize', fontSizeTicks);

    % y軸のスケールを指定
    if strcmp(yScale.(strrep(category, ' ', '_')), 'auto')
        ylim auto;  % 'auto' の場合、自動スケーリング
    else
        ylim(yScale.(strrep(category, ' ', '_')));  % ユーザ指定の範囲
    end

    % レジェンドを右下に配置
    legend('FontSize', fontSizeLabel, 'Location', 'southeast');
    hold off;

    % 画像として保存
    saveas(f, fullfile(outputFolder, [num2str(c) '_' categories{c} '.png']));
    % FIG形式で保存
    savefig(f, fullfile(outputFolder, [num2str(c) '_' categories{c} '.fig']));

    % figureを閉じる
    close(f);
end
