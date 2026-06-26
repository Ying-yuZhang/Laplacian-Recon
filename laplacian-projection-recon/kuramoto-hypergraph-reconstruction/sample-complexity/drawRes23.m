clear; clc; close all;
% ================= 1. 数据准备 =================
data2 = load('Results_M_min_ERHG_2.mat');
data3 = load('Results_M_min_ERHG_3.mat');
N_list = 25:25:100;

% 计算均值 (兼容 3D 和 2D 矩阵)
if ndims(data2.Results_M_min) == 3
    M_mean_2nd = mean(data2.Results_M_min, 3, 'omitnan');
else
    M_mean_2nd = data2.Results_M_min;
end
if ndims(data3.Results_M_min) == 3
    M_mean_3rd = mean(data3.Results_M_min, 3, 'omitnan');
else
    M_mean_3rd = data3.Results_M_min;
end

% ================= 2. 全局设置 =================
algo_names = {'Ours', 'THIS', 'OLS', 'NNLS', 'SL'};
% 为 5 种算法配置统一的颜色矩阵 (1:Ours, 2:THIS, 3:OLS, 4:NNLS, 5:SL)
colors = {
    [0.85, 0.25, 0.25];
    [0.25, 0.50, 0.70];
    [0.45, 0.45, 0.45]; 
    [0.50, 0.65, 0.55];
    [0.60, 0.55, 0.70]  
};
color_H = [0.85, 0.45, 0.15]; % 焦糖橙 (理论边界线)

% ================= 3. 开始绘图 =================
figure('Color', 'w', 'Position', [100, 100, 600, 450]);
hold on; box on;
plot_handles = []; % 动态收集有效图例的句柄

% 循环遍历 5 种算法，统一执行绘图逻辑 (NaN 自动留白)
for i = 1:5
    % 针对 Ours 算法稍微加粗以作强调
    lw = 1.5; ms = 7;
    if i == 1
        lw = 2.0; ms = 8;
    end
    
    % --- 绘制 2 阶网络折线 (实线 + 实心圆点) ---
    h2 = plot(N_list, M_mean_2nd(:, i)', ...
        '-o', 'Color', colors{i}, 'LineWidth', lw, ...
        'MarkerSize', ms, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', ...
        'DisplayName', sprintf('%s (2-order)', algo_names{i}));
    
    % 如果该列不全是 NaN，才将其加入图例显示列表
    if ~all(isnan(M_mean_2nd(:, i)))
        plot_handles = [plot_handles, h2];
    end
    
    % --- 绘制 3 阶网络折线 (虚线 + 空心方块) ---
    h3 = plot(N_list, M_mean_3rd(:, i)', ...
        '--s', 'Color', colors{i}, 'LineWidth', lw, ...
        'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors{i}, ...
        'DisplayName', sprintf('%s (3-order)', algo_names{i}));
        
    % 同样，如果不全是 NaN，加入图例
    if ~all(isnan(M_mean_3rd(:, i)))
        plot_handles = [plot_handles, h3];
    end
end

% --- 绘制理论边界线 ---
H_line = (N_list - 1) + (N_list - 1) .* (N_list - 2) / 2;
h_N = plot(N_list, N_list, '--', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'N');
h_H = plot(N_list, H_line, '-.', 'Color', color_H, 'LineWidth', 1.5, 'DisplayName', 'H');
% plot_handles = [plot_handles, h_N, h_H]; % 将理论线加入图例
plot_handles = [h_N, h_H]; % 将理论线加入图例
% ================= 4. 细节美化 =================
ax = gca;
% 【修改点 1：坐标轴刻度放大、加粗、颜色设为纯黑】
ax.LineWidth = 1.5;             % 边框略微加粗，压住深色字
ax.FontSize = 14;               % 字体调大 (从 12 调至 14)
ax.FontWeight = 'bold';         % 字体加粗
ax.XColor = 'k';                % X轴刻度颜色设为纯黑
ax.YColor = 'k';                % Y轴刻度颜色设为纯黑
ax.TickDir = 'in';
ax.XTick = N_list;

% 开启对数坐标系
set(gca, 'YScale', 'log');
xlim([min(N_list)-2, max(N_list)+2]);

% 【修改点 2：图例字体放大、加粗、颜色设为纯黑】
lgd = legend(plot_handles, 'Location', 'northwest', 'Box', 'off'); 
lgd.FontSize = 13;              % 图例字体调大 (从 11 调至 13)
lgd.FontWeight = 'bold';        % 图例字体加粗
lgd.TextColor = 'k';            % 图例文字设为纯黑

hold off;