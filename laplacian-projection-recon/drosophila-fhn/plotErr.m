% % 数据准备 (假设 err.mat 已加载)
load F1_Score.mat
data = F1_Score;
x = 1:size(data,2);  % x 轴数据
y1 = data(1, :);
y2 = data(2, :);
y3 = data(3, :);
y4 = data(4, :);

% =========================================================================
% 1. 颜色配置：与图(a)保持绝对一致的莫兰迪 + 赤陶红体系
% =========================================================================
color_1st = [123, 164, 182] / 255; % 灰蓝色 (1阶)
color_2nd = [138, 155, 110] / 255; % 鼠尾草绿 (2阶)
color_3rd = [217, 182, 110] / 255; % 柔和沙黄 (3阶)
color_4th = [176,  92,  80] / 255; % 赤陶红 (4阶)

% =========================================================================
% 2. 图形初始化与全局设置
% =========================================================================
figure('Color', 'w'); % 确保纯白背景
hold on;

% =========================================================================
% 3. 绘制数据线 (误差下降折线)
% 使用与图(a)相同的线宽 1.5 和实心标记点设置
% =========================================================================
l_width = 1.5; 
m_size  = 5;   

p1 = plot(x, y1, '-o', 'LineWidth', l_width, 'Color', color_1st, 'MarkerFaceColor', color_1st, 'MarkerSize', m_size);
p2 = plot(x, y2, '-s', 'LineWidth', l_width, 'Color', color_2nd, 'MarkerFaceColor', color_2nd, 'MarkerSize', m_size);
p3 = plot(x, y3, '-^', 'LineWidth', l_width, 'Color', color_3rd, 'MarkerFaceColor', color_3rd, 'MarkerSize', m_size);
p4 = plot(x, y4, '-d', 'LineWidth', l_width, 'Color', color_4th, 'MarkerFaceColor', color_4th, 'MarkerSize', m_size);

% =========================================================================
% 4. 坐标轴、图例与排版细节优化
% =========================================================================
% 图例设置：首字母大写，无边框，并使用句柄精确绑定
lgd = legend([p1, p2, p3, p4], '1st-Order', '2nd-Order', '3rd-Order', '4th-Order', 'Location', 'best');
legend('boxoff');
lgd.FontSize = 11; 
lgd.FontName = 'Arial'; 

% 标签与坐标轴设置
% xlabel('Iterations', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
% ylabel('Reconstruction Errors', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold'); % 补全原图中 Y 轴标签

xticks(x);

% 边框与刻度：封闭边框，刻度向内，无网格，维持顶级期刊质感
box on;
grid off;
set(gca, 'TickDir', 'in', 'LineWidth', 0.75, 'FontName', 'Arial', 'FontSize', 10);

% 窗口尺寸：使用您指定的尺寸
set(gcf, 'Position', [200, 200, 400, 320]); 

hold off;