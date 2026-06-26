function metrics = clc_alg_performance(AA, A00, A2, A3, this_AA, ours_AA)
% CLC_ALG_PERFORMANCE Calculate algorithm performance (9x5 matrix)
% CLC_ALG_PERFORMANCE 计算算法性能 (9x5 矩阵)
% 
% Input/输入: 
%   AA: Ground Truth / 真实值
%   A00, A2, A3, this_AA, ours_AA: Estimated values / 各算法的估计值
%
% Output/输出 metrics (9 rows x 5 columns):
%   Rows 1-3: [AUROC; AUPR; F1] (Overall / 整体)
%   Rows 4-6: [AUROC; AUPR; F1] (Pairwise / 2-超边: Top N-1 rows)
%   Rows 7-9: [AUROC; AUPR; F1] (Triadic / 3-超边: Remaining rows)

    % --- 0. Basic Parameter Definition / 基础参数定义 ---
    N = size(AA, 2);
    metrics = zeros(9, 5);

    % =========================================================
    % 1. Overall Performance (ALL) / 整体性能
    %    Rows 1-3
    % =========================================================
    metrics(1:3, 1) = get_metrics(AA, A00,     false);
    metrics(1:3, 2) = get_metrics(AA, A2,      false);
    metrics(1:3, 3) = get_metrics(AA, A3,      false);
    metrics(1:3, 4) = get_metrics(AA, this_AA, false);
    
    % ours_AA outputs 0/1 predictions, no threshold truncation needed.
    % ours_AA 输出的是0/1的预测值不需要阈值截断。
    metrics(1:3, 5) = get_metrics(AA, ours_AA, true);  

    % =========================================================
    % 2. Pairwise Performance (2-Hyperedge) / 2-超边性能
    %    Rows 4-6 (Top N-1 rows / 前 N-1 行)
    % =========================================================
    idx_pair = 1:(N-1);
    gt_pair  = AA(idx_pair, :);
    
    metrics(4:6, 1) = get_metrics(gt_pair, A00(idx_pair, :),     false);
    metrics(4:6, 2) = get_metrics(gt_pair, A2(idx_pair, :),      false);
    metrics(4:6, 3) = get_metrics(gt_pair, A3(idx_pair, :),      false);
    metrics(4:6, 4) = get_metrics(gt_pair, this_AA(idx_pair, :), false);
    
    % ours_AA outputs 0/1 predictions, no threshold truncation needed.
    % ours_AA 输出的是0/1的预测值不需要阈值截断。
    metrics(4:6, 5) = get_metrics(gt_pair, ours_AA(idx_pair, :), true);

    % =========================================================
    % 3. Triadic Performance (3-Hyperedge) / 3-超边性能
    %    Rows 7-9 (Remaining rows / 剩余行)
    % =========================================================
    idx_tri = N:size(AA, 1);
    gt_tri  = AA(idx_tri, :);
    
    metrics(7:9, 1) = get_metrics(gt_tri, A00(idx_tri, :),     false);
    metrics(7:9, 2) = get_metrics(gt_tri, A2(idx_tri, :),      false);
    metrics(7:9, 3) = get_metrics(gt_tri, A3(idx_tri, :),      false);
    metrics(7:9, 4) = get_metrics(gt_tri, this_AA(idx_tri, :), false);
    
    % ours_AA outputs 0/1 predictions, no threshold truncation needed.
    % ours_AA 输出的是0/1的预测值不需要阈值截断。
    metrics(7:9, 5) = get_metrics(gt_tri, ours_AA(idx_tri, :), true);

end

% ---------------------------------------------------------
% Helper Function: Calculate [ROC; PR; F1] for a single group
% 内部辅助函数: 计算单组数据的 [ROC; PR; F1]
% 
% Inputs:
%   gt: Ground Truth
%   score: Estimated Score
%   is_binary: Boolean flag for ours_AA (True if input is already 0/1)
% ---------------------------------------------------------
function res = get_metrics(gt, score, is_binary)
    res = zeros(3, 1);
    
    % 1. Calculate AUROC and AUPR (Always based on raw scores)
    %    计算 AUROC 和 AUPR (始终基于原始分数)
    [res(1), res(2)] = clc_AUROC_AUPR(score(:), gt(:));
    
    % 2. Calculate F1 / 计算 F1
    if is_binary
        % ours_AA outputs 0/1 predictions, no threshold truncation needed.
        % ours_AA 输出的是0/1的预测值不需要阈值截断。
        bin_pred = score;
    else
        % Other algorithms need Top-K selection
        % 其他算法需要筛选 Top-K
        bin_pred = selectTopByMask2(gt, score);
    end
    
    res(3) = clc_F1(gt, bin_pred);
end