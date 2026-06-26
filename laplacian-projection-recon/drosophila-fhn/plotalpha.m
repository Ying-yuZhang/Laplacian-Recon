% % 数据准备 (假设 G.mat 已加载)
load("G.mat")
% 确保提取前 5 行数据对应 1 到 5 阶
data = G.alpha(1:5, 1:6);
x = 0:size(data,2)-1;  % x 轴数据
y1 = data(1, :);
y2 = data(2, :);
y3 = data(3, :);
y4 = data(4, :);
y5 = data(5, :); % 【补充】提取第 5 阶数据

% =========================================================================
% 1. 颜色配置：平滑渐变法 (将赤陶红移至 5 阶，在 3-5 之间插入橙色调)
% =========================================================================
color_1st = [123, 164, 182] / 255; % 灰蓝色 (1阶)
color_2nd = [138, 155, 110] / 255; % 鼠尾草绿 (2阶)
color_3rd = [217, 182, 110] / 255; % 柔和沙黄 (3阶)
color_4th = [220, 135,  95] / 255; % 琥珀灰橘 (4阶，暖色过渡)
color_5th = [176,  92,  80] / 255; % 赤陶红 (5阶，绝对核心高亮)

% =========================================================================
% 2. 图形初始化与全局设置
% =========================================================================
figure('Color', 'w');
hold on;

% =========================================================================
% 3. 绘制参考线 (虚线，真实值/目标值)
% 【补充】添加第 5 阶的参考线 (假设递进值为 0.6，如果 F1 上限是 1.0 请根据实际情况修改)
% =========================================================================
ref_width = 1.2; 
yline(0.2, '--', 'LineWidth', ref_width, 'Color', color_1st);
yline(0.3, '--', 'LineWidth', ref_width, 'Color', color_2nd);
yline(0.4, '--', 'LineWidth', ref_width, 'Color', color_3rd);
yline(0.5, '--', 'LineWidth', ref_width, 'Color', color_4th);
yline(0.6, '--', 'LineWidth', ref_width, 'Color', color_5th); 

% =========================================================================
% 4. 绘制数据线 (实线，估计值)
% 【优化】引入线宽递进机制：低阶细线，高阶粗线，配合颜色进一步拉开层次
% =========================================================================
lw_thin = 1.5;  % 低阶弱化线宽
lw_mid  = 1.5;  % 中阶常规线宽
lw_bold = 1.5;  % 高阶强调线宽
m_size  = 5;   

% 标记点：1-o(圆), 2-s(方), 3-^(上三角), 4-d(菱形), 5-p(五角星)
p1 = plot(x, y1, '-o', 'LineWidth', lw_thin, 'Color', color_1st, 'MarkerFaceColor', color_1st, 'MarkerSize', m_size);
p2 = plot(x, y2, '-s', 'LineWidth', lw_thin, 'Color', color_2nd, 'MarkerFaceColor', color_2nd, 'MarkerSize', m_size);
p3 = plot(x, y3, '-^', 'LineWidth', lw_mid,  'Color', color_3rd, 'MarkerFaceColor', color_3rd, 'MarkerSize', m_size);
p4 = plot(x, y4, '-d', 'LineWidth', lw_bold, 'Color', color_4th, 'MarkerFaceColor', color_4th, 'MarkerSize', m_size);
p5 = plot(x, y5, '-v', 'LineWidth', lw_bold, 'Color', color_5th, 'MarkerFaceColor', color_5th, 'MarkerSize', m_size);

% =========================================================================
% 5. 坐标轴、图例与排版细节优化
% =========================================================================
% 【补充】将第 5 阶 (p5) 加入图例
lgd = legend([p1, p2, p3, p4, p5], '1st-Order', '2nd-Order', '3rd-Order', '4th-Order', '5th-Order', 'Location', 'best');
legend('boxoff');
lgd.FontSize = 11; 
lgd.FontName = 'Arial'; 

% % 取消注释，并设定对应 F1 的专业标签
% xlabel('Iterations', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
% ylabel('Coupling Strength', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
% xticks(x);

% 边框与刻度
box on;
grid off;
set(gca, 'TickDir', 'in', 'LineWidth', 0.75, 'FontName', 'Arial', 'FontSize', 10);

% 窗口尺寸
set(gcf, 'Position', [200, 200, 400, 320]);

hold off;