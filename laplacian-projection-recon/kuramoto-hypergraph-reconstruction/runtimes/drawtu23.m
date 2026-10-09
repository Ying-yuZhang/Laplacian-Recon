clear; clc; close all;
% ================= 1. 定义实验参数与配色 =================
% --- 修改核心 1：分别定义 2 阶和 3 阶的网络规模 ---
N_list_2nd = [20, 40, 50:50:200];
N_list_3rd = [20, 40, 50:50:200];
% 严格定义所有的 5 种算法名称
algos = {'Ours', 'OLS', 'NNLS', 'SL', 'THIS'}; 
% 严格对应的 5 种配色
colors = [
    0.85, 0.25, 0.25;  % 1-Ours - 朱砂红
    0.45, 0.45, 0.45;  % 2-OLS  - 高级灰
    0.50, 0.65, 0.55;  % 3-NNLS - 鼠尾草绿
    0.60, 0.55, 0.70;  % 4-SL   - 尘埃紫
    0.25, 0.50, 0.70;  % 5-THIS - 钢蓝色
];
% 定义标记形状以区分阶数 (Order)
marker_2nd = 'o'; % 2阶 - 圆圈
marker_3rd = 's'; % 3阶 - 方块

% ================= 2. 加载与处理数据 =================
file_2nd = 'runtimes_2000_HG_2.mat';   % 包含 5 种算法
file_3rd = 'runtimes_2000_HG_3.mat';   % 包含 3 种算法
% --- 处理 2 阶数据 (5种算法) ---
if isfile(file_2nd)
    data_2nd = load(file_2nd);
    time_2nd = mean(data_2nd.runtimes(:,:,2:end), 3); %计算均值时剔除掉首次启动的预热数据。
else
    error(['找不到 2阶 数据文件: ', file_2nd]);
end
% --- 处理 3 阶数据 (3种算法) ---
if isfile(file_3rd)
    data_3rd = load(file_3rd);
    time_3rd = mean(data_3rd.runtimes(:,:,2:end),3); %计算均值时剔除掉首次启动的预热数据。
else
    error(['找不到 3阶 数据文件: ', file_3rd]);
end

% ================= 3. 绘制主图 =================
fig = figure('Color', 'w', 'Position', [100, 100, 650, 500]); 
ax_main = axes(fig, 'Position', [0.10, 0.12, 0.85, 0.72]); 
hold(ax_main, 'on'); box(ax_main, 'on');
h_2nd = gobjects(1, 5); % 预分配 2-order 句柄数组
h_3rd = gobjects(1, 5); % 预分配 3-order 句柄数组

% --- A. 绘制 2 阶网络曲线 (使用 N_list_2nd) ---
for i = 1:5
    plot(ax_main, N_list_2nd, time_2nd(:, i), '-', ...
        'Color', [colors(i, :) 0.6], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    h_2nd(i) = scatter(ax_main, N_list_2nd, time_2nd(:, i), 80, colors(i, :), ...
        marker_2nd, 'filled', ...
        'MarkerEdgeColor', 'none', ...
        'DisplayName', [algos{i} ' (2-order)']);
end

% --- B. 绘制 3 阶网络曲线 (使用 N_list_3rd) ---
for i = 1:5
    plot(ax_main, N_list_3rd, time_3rd(:, i), '--', ...
        'Color', [colors(i, :) 0.6], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    h_3rd(i) = scatter(ax_main, N_list_3rd, time_3rd(:, i), 80, colors(i, :), ...
        marker_3rd, ...
        'MarkerEdgeColor', colors(i, :), 'LineWidth', 1.5, ...
        'DisplayName', [algos{i} ' (3-order)']);
end

% ================= 修改核心：交替重组句柄 =================
h_matrix = [h_2nd; h_3rd];
h_plots = h_matrix(:)'; 

% ================= 4. 美化坐标轴 =================
set(ax_main, 'XScale', 'log', 'YScale', 'log', 'FontName', 'Arial'); 
% --- 修改核心 2：动态合并并排序所有的 N 值用于坐标轴范围与刻度 ---
all_N = sort(unique([N_list_2nd, N_list_3rd])); % [20, 40, 50, 80, 100, 150, 160, 200]
xlim(ax_main, [min(all_N)*0.85, max(all_N)*1.15]);
all_times = [time_2nd(:); time_3rd(:)];
valid_times = all_times(all_times > 0); 
ylim(ax_main, [min(valid_times)*0.5, 1e5]); 

% 显示合并后的所有刻度点
xticks(ax_main, all_N);
xticklabels(ax_main, string(all_N));
% 因为刻度点变多 (如 150 和 160 靠得很近)，为防止重叠，将标签倾斜 45 度
xtickangle(ax_main, 45); 

% ---> 【坐标轴修改点：字体变大、加粗、颜色设纯黑】 <---
ax_main.LineWidth = 1.5;             % 略微加粗边框以匹配深色字体
ax_main.FontSize = 16;               % 将刻度字体从 12 调大到 14
ax_main.FontWeight = 'bold';         % 刻度字体加粗，显得颜色更深
ax_main.XColor = 'k';                % 强制 X轴 颜色为纯黑
ax_main.YColor = 'k';                % 强制 Y轴 颜色为纯黑
ax_main.TickDir = 'in'; 
ax_main.TickLength = [0.015, 0.015];

% % ================= 5. 定制图例 =================
% leg = legend(ax_main, h_plots, 'Location', 'northoutside', 'NumColumns', 5, 'Interpreter', 'none');
% leg.Box = 'off';

% ---> 【字体变大、加粗、文字设纯黑】 <---
leg.FontSize = 16;                   % 将图例字体从 10 调大到 12
leg.FontWeight = 'bold';             % 图例字体加粗
leg.TextColor = 'k';                 % 强制图例文字为纯黑
leg.ItemTokenSize = [15, 10]; 

hold(ax_main, 'off');