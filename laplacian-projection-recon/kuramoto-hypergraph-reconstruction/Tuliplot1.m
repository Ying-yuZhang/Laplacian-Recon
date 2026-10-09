%% 0. 基础设置与数据准备
clear
clc
N = 20; 
load perfor_3_2000_1_5.mat
% 适配 irange (根据你的代码)
irange = unique(round(logspace(log10(3), log10(2000), 30)));
len = length(irange); 
% 计算 M/N 比率 (X轴数据)
ratio_range = irange / N;
x_min_r = min(ratio_range);
x_max_r = max(ratio_range);
% 计算平均性能 (对第4维 trials 求均值) -> 结果维度 [9, 5, len]
% 注意：假设 perfor 在你的工作区中已经存在
valid_perfor = mean(perfor(:, :, 1:len, :), 4); 

%% 1. 样式与配置
algos = {'OLS', 'NNLS', 'SL', 'THIS', 'Ours'};
% 定义 3 种指标
metric_names = {'AUC', 'AUPR', 'F1-Score'};
% 定义 3 种区域类型 (对应 Rows 1-3, 4-6, 7-9)
region_names = {'Overall', 'Pairwise (2-Hyperedge)', 'Triadic (3-Hyperedge)'};
% 颜色与线型
custom_colors = [
    0.45, 0.45, 0.45;  % 1. OLS  - 高级灰 (Neutral Grey)
    0.50, 0.65, 0.55;  % 2. NNLS - 鼠尾草绿 (Sage Green)
    0.60, 0.55, 0.70;  % 3. SL   - 尘埃紫 (Dusty Purple)
    0.25, 0.50, 0.70;  % 4. THIS - 钢蓝色 (Steel Blue): 强对比冷色调
    0.85, 0.25, 0.25   % 5. Ours - 朱砂红 (Vermilion): 核心高亮
];
custom_markers = {'o', 's', '^', 'v', 'd'};
custom_lineStyles = {'--', '--', '--', '-', '-'};

% 刻度设置 (原始 M 值 -> 转换为 M/N)
tick_values_M = [3, 10, 30, 100, 300, 1000]; 
tick_values_r = tick_values_M / N;

%% 2. 9幅独立绘图循环
% i 对应 perfor 的第 1 维 (1-9)
for i = 1
    % --- 2.1 计算当前图表的属性 ---
    metric_idx = mod(i-1, 3) + 1;
    metric_name = metric_names{metric_idx};
    
    region_idx = ceil(i / 3);
    region_name = region_names{region_idx};
    
    % 创建 Figure
    figure('Color', 'w', 'Position', [100, 100, 600, 500], ...
           'Name', sprintf('%s - %s', region_name, metric_name)); 
    hold on;
    
    % 提取当前数据 [5, len]
    current_data = squeeze(valid_perfor(i, :, :));
    
    % --- 2.2 坐标轴核心设置 ---
    ax = gca;
    set(ax, 'XScale', 'log');
    xlim([x_min_r, x_max_r]); 
    
    % 调整坐标轴的字体大小和颜色
    ax.Box = 'on';
    ax.LineWidth = 1.5;             % 稍微加粗一点边框，配合深色字体
    ax.FontSize = 16;               % 将字体调大 (原为 11)
    ax.FontWeight = 'bold';         % 字体加粗，显得颜色更深更实
    ax.XColor = 'k';                % X轴刻度及数字设为纯黑
    ax.YColor = 'k';                % Y轴刻度及数字设为纯黑
    ax.TickDir = 'in';
    
    % Y 轴范围自适应 
    y_min = min(current_data(:));
    y_max = max(current_data(:));
    ylim([max(0, y_min*0.9), min(1.0, y_max*1.02)]);
    
    % --- 2.4 绘制 5 种算法曲线 ---
    for j = 1:5
        plot(ratio_range, current_data(j, :), ...
            'LineStyle', custom_lineStyles{j}, ...
            'LineWidth', 1.8, ...
            'Color', custom_colors(j, :), ...
            'Marker', custom_markers{j}, ...
            'MarkerSize', 6, ...         % 略微增大了Marker，让画面更平衡
            'MarkerFaceColor', 'w', ...
            'DisplayName', algos{j});
    end
    
    % --- 2.5 细节打磨 ---
    % 刻度设置
    xticks(tick_values_r);
    xticklabels(arrayfun(@(x) sprintf('%d', x), tick_values_M, 'UniformOutput', false));
    
    % [修改点] 调整图例的字体大小和颜色
    % 将 FontSize 从 9 调大到 12，并增加粗体和纯黑色设置
    lgd = legend('Location', 'southwest', 'FontSize', 16); 
    lgd.FontWeight = 'bold';        % 图例字体加粗
    lgd.TextColor = 'k';            % 图例文字纯黑
    lgd.EdgeColor = 'none';
    
    hold off;
end