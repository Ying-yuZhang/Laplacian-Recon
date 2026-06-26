% =========================================================================
% 实验：认知流的时序动态追踪 (Grand Average 组级别平均分析)
% 自动化批处理 20 名受试者，并绘制高阶网络质心的组平均轨迹
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

% 2. 提前加载节点三维坐标 (用于计算超边质心)
% 请确保工作区或当前目录下有 coord_ROI.mat，包含 82x4 的 coord 变量
% 它在数据集的\ds003505\code\Source-Reconstruction-Matlab\目录下
load coord_ROI.mat
xyz_coords=coord(:,1:3);

% 3. 预分配矩阵，用于存储所有受试者在 4 个时间窗的 Y 轴质心
% 大小为 20 x 4，初始化为 NaN (方便后续过滤无效数据)
all_Y_cog = NaN(num_subjects, num_windows); 

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 四个时间窗的中点索引计算 (1个点 = 4ms):
    % 刺激前基线 (Baseline)               W1 (-500ms) = -125 个点
    % 早期视觉响应 (Early Visual)         W2 ( 150ms) = +38 个点 (150/4 = 37.5，四舍五入)
    % 运动特征整合 (Motion Integration)   W3 ( 300ms) = +75 个点
    % 运动决策与执行(Decision Motor)       W4 ( 500ms) = +125 个点
    tp_list = round([stim_idx_X - 125, stim_idx_X + 38, stim_idx_X + 75, stim_idx_X + 125]);
    
    % 存储该受试者的 4 个时间窗的超边数据
    tris_curr_sub = cell(1, 4);
    T = size(X, 1);
    
    % --- 四个时间窗的网络重构 ---
    for w = 1:num_windows
        tp = tp_list(w);
        selected_indices = tp-70:tp+70;
        
        Xp = X(selected_indices, :) - X(tp, :);
        Yp = (Y(:, selected_indices) - Y(:, tp))';
        
        X_score_tmp = sol_lasso(Xp, Yp, 1, 1e-4, 100000, 10^-5);
        X_score = aggregate_multidim(X_score_tmp, 'max');
        X_score = 0.5 * (X_score + X_score');
        
        X_score_upper = X_score(upper_triangle_mask);
        X_score_upper = (X_score_upper - min(X_score_upper)) / (max(X_score_upper) - min(X_score_upper));
        X_score_norm = zeros(size(X_score));
        X_score_norm(upper_triangle_mask) = X_score_upper;
        X_score_norm = X_score_norm + X_score_norm';
        
        ite_num = 3; cut = 0.3; delta = 0.02; init_alpha = [0.4 1]; order = 2;
        G = Extract_simplex(X_score_norm, ite_num, init_alpha, delta, cut ,order);
        
        % 提取当前窗口的超边矩阵
        [~, tris_curr_sub{w}] = convert_to_list(G.A_estimate, G.A2_estimate);
    end
    
    % --- 核心计算：提取该受试者在 4 个时间窗的 Y 轴质心 ---
    for w = 1:num_windows
        curr_tris = tris_curr_sub{w};
        if ~isempty(curr_tris)
            active_nodes_weighted = curr_tris(:); 
            cog_xyz = mean(xyz_coords(active_nodes_weighted, :), 1);
            all_Y_cog(sub_num, w) = cog_xyz(2); % 仅提取 Y 轴 (前后) 坐标存入矩阵
        end
    end
end
fprintf('\n所有受试者数据跑完！开始进行平均统计...\n');

% =========================================================================
% 数据清洗与统计：大组平均 (Grand Average)
% =========================================================================
% 剔除那些在某个时间窗没有算出超边（产生 NaN）的受试者
valid_data = all_Y_cog(~any(isnan(all_Y_cog), 2), :);
num_valid = size(valid_data, 1);
fprintf('共有 %d 名受试者的数据有效参与了平均计算。\n', num_valid);

% 计算平均值 (Mean) 和 标准误 (Standard Error of Mean, SEM)
mean_Y = mean(valid_data, 1);
sem_Y  = std(valid_data, 0, 1) / sqrt(num_valid);
% =========================================================================
% 附加分析：计算满足“质心前移”的个体一致性 (Responder Analysis)
% =========================================================================
% 定义质心前移：W4 (决策执行, 前侧) 的 Y 坐标 > W2 (视觉响应, 后侧) 的 Y 坐标
forward_shift_idx = valid_data(:, 4) > valid_data(:, 2);

% 计算满足条件的人数和百分比
num_forward = sum(forward_shift_idx);
percent_forward = (num_forward / num_valid) * 100;

% 额外看一眼 W4 > W3 (动作执行 > 运动整合) 的延迟跃升人数
delayed_shift_idx = valid_data(:, 4) > valid_data(:, 3);
num_delayed = sum(delayed_shift_idx);
percent_delayed = (num_delayed / num_valid) * 100;

% 打印结果到命令窗口
fprintf('\n================ 质心转移个体一致性分析 ================\n');
fprintf('总有效受试者: %d 名\n', num_valid);
fprintf('--------------------------------------------------------\n');
fprintf('【宏观转移 (W4 > W2)】\n');
fprintf('质心成功从后侧视觉区移向前侧运动区的受试者: %d 名 (占比 %.1f%%)\n', num_forward, percent_forward);
fprintf('--------------------------------------------------------\n');
fprintf('【晚期跃升 (W4 > W3)】\n');
fprintf('在动作执行前夕发生质心前向突增的受试者: %d 名 (占比 %.1f%%)\n', num_delayed, percent_delayed);
fprintf('========================================================\n');

% =========================================================================
% 绘制 Grand Average 组平均趋势图 (带误差棒)
% =========================================================================
figure('Color', 'w', 'Position', [150, 150, 600, 450], 'Name', 'Grand Average Y-Axis Shift');
hold on;

% 绘制带误差棒的折线图 (保持高对比度橙色)
errorbar(1:num_windows, mean_Y, sem_Y, '-o', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', '#D95F02', 'Color', '#D95F02', ...
    'CapSize', 10); % CapSize 控制误差棒两端横线的宽度

% 坐标轴刻度与标签设置 (封闭方格、刻度向内、Arial字体)
set(gca, 'XTick', 1:num_windows, 'XTickLabel', win_names, ...
    'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold', ...
    'TickDir', 'in', 'LineWidth', 1.2);

% ylabel('Average Centroid Y-coordinate', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');


% 开启封闭边框与网格
box on; 
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.25);

% X 轴范围调整，为左右两端留出 0.5 的视觉空隙
xlim([0.7, num_windows + 0.3]);

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