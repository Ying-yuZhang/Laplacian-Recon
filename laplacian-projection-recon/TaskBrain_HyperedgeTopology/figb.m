% =========================================================================
% 实验一：认知流的时序动态追踪 (Grand Average 组级别平均分析)
% 焦点通路版：自动化批处理 20 名受试者，绘制三条关键通路的组平均趋势与误差棒
% =========================================================================

% 强制设置字体，防止中文乱码
set(groot, 'defaultAxesFontName', 'Microsoft YaHei');
set(groot, 'defaultTextFontName', 'Microsoft YaHei');

% 1. 定义全局参数与路径
base_dir = 'C:\##########\ds003505\derivatives\vepcon_esi-v1.1';
num_subjects = 20;
num_windows = 4;
win_names = {'W1 (-500ms)', 'W2 (150ms)', ...
             'W3 (300ms)', 'W4 (500ms)'};

% 2. 预分配矩阵，用于存储所有受试者在 4 个时间窗的三条通路权重
% 大小均为 20 x 4，初始化为 NaN (方便后续过滤无效数据)
all_vis_vals = NaN(num_subjects, num_windows); % 左脑视觉流 (23-30)
all_mot_vals = NaN(num_subjects, num_windows); % 决策运动流 (19-10)
all_exe_vals = NaN(num_subjects, num_windows); % 按键执行流 (10-16)

% =========================================================================
% 开启 20 名受试者的全自动批处理循环
% =========================================================================
for sub_num = 1:num_subjects
    sub_id = sprintf('sub-%02d', sub_num);
    file_name = sprintf('%s_task-motion_label-L2018_desc-scale1_rtc.h5', sub_id);
    file_path = fullfile(base_dir, sub_id, 'eeg', file_name);
    
    % 安全检查：如果该被试数据丢失，跳过并记录，程序不会报错中断
    if ~exist(file_path, 'file')
        warning('未找到文件，跳过被试: %s', sub_id);
        continue;
    end
    
    fprintf('\n>>>>>>>>>> 正在全自动处理被试: %s <<<<<<<<<<\n', sub_id);
    
    % --- 数据加载与预处理 ---
    esi_data = h5read(file_path, '/esi');
    [num_nodes, ~, ~] = size(esi_data);
    Fs = 250; dt = 1 / Fs;
    
    y_plus2  = esi_data(:, 5:end, :);       
    y_plus1  = esi_data(:, 4:end-1, :);     
    y_minus1 = esi_data(:, 2:end-3, :);     
    y_minus2 = esi_data(:, 1:end-4, :);     
    esi_derivative = (-y_plus2 + 8 * y_plus1 - 8 * y_minus1 + y_minus2) / (12 * dt);
    
    esi_data_aligned = esi_data(:, 3:end-2, :); 
    ERP_data = mean(esi_data_aligned, 3);
    
    m = 2; tau = 1;
    X = phase_space_reconstruction(ERP_data, m, tau);
    t_start = 1 + (m - 1) * tau;
    Y = mean(esi_derivative(:, t_start:end, :), 3);
    
    upper_triangle_mask = triu(true(num_nodes, num_nodes), 1);
    
    % --- 精准时间轴对齐 ---
    stim_idx_raw = 375; 
    offset = 2 + (m - 1) * tau; 
    stim_idx_X = stim_idx_raw - offset; 
    
    tp_list = round([stim_idx_X - 125, stim_idx_X + 38, stim_idx_X + 75, stim_idx_X + 125]);
    T = size(X, 1);
    
    % --- 四个时间窗的网络重构 ---
    for w = 1:num_windows
        tp = tp_list(w);
        selected_indices = tp-70:tp+70;
        
        Xp = X(selected_indices, :) - X(tp, :);
        Yp = (Y(:, selected_indices) - Y(:, tp))';
        
        % 稀疏回归与网络推断
        X_score_tmp = sol_lasso(Xp, Yp, 1, 1e-4, 100000, 10^-5);
        X_score = aggregate_multidim(X_score_tmp, 'max');
        X_score = 0.5 * (X_score + X_score');
        
        % 归一化提取 (0 到 1)
        X_score_upper = X_score(upper_triangle_mask);
        X_score_upper = (X_score_upper - min(X_score_upper)) / (max(X_score_upper) - min(X_score_upper));
        X_score_norm = zeros(size(X_score));
        X_score_norm(upper_triangle_mask) = X_score_upper;
        X_score_norm = X_score_norm + X_score_norm';
        
        % 直接从归一化的连续网络中提取三条核心通路的权重
        all_vis_vals(sub_num, w) = X_score_norm(23, 30);%[Visual pathway] LOC -> MTG 
        all_mot_vals(sub_num, w) = X_score_norm(19, 10);%[Parieto-precentral pathway] IPL -> PreCG
        all_exe_vals(sub_num, w) = X_score_norm(10, 16);%[Sensorimotor pathway] PreCG -> PoCG
    end
end

fprintf('\n所有受试者数据跑完！开始进行大组平均统计与绘图...\n');

% =========================================================================
% 数据清洗与统计：大组平均 (Grand Average)
% =========================================================================
% 找出在所有三个变量中都没有 NaN 的“完全有效”受试者索引
valid_idx = ~any(isnan(all_vis_vals), 2) & ~any(isnan(all_mot_vals), 2) & ~any(isnan(all_exe_vals), 2);
num_valid = sum(valid_idx);
fprintf('共有 %d 名受试者的数据有效参与了平均计算。\n', num_valid);

% 提取有效数据
vis_valid = all_vis_vals(valid_idx, :);
mot_valid = all_mot_vals(valid_idx, :);
exe_valid = all_exe_vals(valid_idx, :);

% 计算平均值 (Mean)
mean_vis = mean(vis_valid, 1);
mean_mot = mean(mot_valid, 1);
mean_exe = mean(exe_valid, 1);

% 计算标准误 (Standard Error of Mean, SEM = std / sqrt(N))
sem_vis = std(vis_valid, 0, 1) / sqrt(num_valid);
sem_mot = std(mot_valid, 0, 1) / sqrt(num_valid);
sem_exe = std(exe_valid, 0, 1) / sqrt(num_valid);

% =========================================================================
% 绘制学术级 Grand Average 关键通路折线图 (带误差棒)
% =========================================================================
figure('Color', 'w', 'Position', [150, 150, 600, 450], 'Name', 'Grand Average Sensorimotor Pathways');
hold on;

% 设定高对比度学术配色
color_visual = '#E34234'; % 朱红色
color_motor  = '#4682B4'; % 钢蓝色
color_exec   = '#D95F02'; % 橙棕色

% 绘制视觉流误差棒图
errorbar(1:num_windows, mean_vis, sem_vis, '-o', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', color_visual, ...
    'Color', color_visual, ...
    'CapSize', 8, ...
    'DisplayName', 'Visual Processing: LOC \rightarrow MTG');

% 绘制决策流误差棒图
errorbar(1:num_windows, mean_mot, sem_mot, '-s', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', color_motor, ...
    'Color', color_motor, ...
    'CapSize', 8, ...
    'DisplayName', 'Visuomotor Transformation: IPL \rightarrow PreCG');

% 绘制感觉运动流误差棒图
errorbar(1:num_windows, mean_exe, sem_exe, '-^', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', color_exec, ...
    'Color', color_exec, ...
    'CapSize', 8, ...
    'DisplayName', 'Sensorimotor Integration: PreCG \leftrightarrow PoCG');

% 坐标轴设置
xtick_labels = {'W1 (-500 ms)', 'W2 (150 ms)', 'W3 (300 ms)', 'W4 (500 ms)'};

% 刻度线向内 (TickDir = 'in'), 增加线宽以匹配学术风格
set(gca, 'XTick', 1:num_windows, 'XTickLabel', xtick_labels, ...
    'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold', ...
    'TickDir', 'in', 'LineWidth', 1.2);

% 图例设置与网格
legend('Location', 'northwest', 'FontSize', 11, 'FontName', 'Arial', 'FontWeight', 'bold', 'Box', 'off');

% 封闭方格与网格样式
box on; 
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.25);

% X 轴范围调整，为两端留出 0.5 的空隙
xlim([0.7, num_windows + 0.3]);

% 动态自适应 Y 轴范围，为误差棒留出空间
all_means = [mean_vis, mean_mot, mean_exe];
all_sems  = [sem_vis, sem_mot, sem_exe];
ylim([min(all_means - all_sems) - 0.05, max(all_means + all_sems) + 0.1]);

hold off;

% ---------------- 辅助函数定义 ----------------
function A_pred = aggregate_multidim(X_score, method)
    if nargin < 2; method = 'max'; end
    N = size(X_score, 2);
    m = size(X_score, 1) / N;
    if mod(size(X_score, 1), N) ~= 0
        error('聚合失败：X_score 的行数必须是列数 N 的整数倍！');
    end
    X_3D = reshape(X_score, m, N, N);
    switch lower(method)
        case 'max'
            A_pred = squeeze(max(X_3D, [], 1));
        case 'mean'
            A_pred = squeeze(mean(X_3D, 1));
        otherwise
            error('不支持的聚合方法！请传入 ''max'' 或 ''mean''。');
    end
end