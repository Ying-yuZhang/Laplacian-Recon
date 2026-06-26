% =========================================================================
% 通用团重构 F1-Score 计算函数
% 输入:
%   pred_cliques - 算法重构出的团矩阵 (N1 x d)
%   true_cliques - Ground Truth 真实的团矩阵 (N2 x d)
% 输出:
%   f1_score     - F1 分数 (0 到 1 之间)
%   precision    - 准确率 (查准率)
%   recall       - 召回率 (查全率)
%   TP, FP, FN   - 真正例, 假正例, 假负例的数量
% =========================================================================
function [f1_score, precision, recall, TP, FP, FN] = calc_f1(pred_cliques, true_cliques)
    
    % 1. 边界情况处理 (其中一个为空或都为空)
    if isempty(pred_cliques) && isempty(true_cliques)
        f1_score = 1; precision = 1; recall = 1;
        TP = 0; FP = 0; FN = 0;
        return;
    elseif isempty(pred_cliques)
        f1_score = 0; precision = 0; recall = 0;
        TP = 0; FP = 0; FN = size(true_cliques, 1);
        return;
    elseif isempty(true_cliques)
        f1_score = 0; precision = 0; recall = 0;
        TP = 0; FP = size(pred_cliques, 1); FN = 0;
        return;
    end

    % 2. 内部排序以保证无序集合的可比性
    % 将每一行(每个团)内部的节点编号从小到大排序，杜绝排列带来的误判
    P = sort(pred_cliques, 2);
    T = sort(true_cliques, 2);

    % 3. 寻找交集 (True Positives)
    % intersect(..., 'rows') 会找出两个矩阵中完全相同的行
    TP_matrix = intersect(P, T, 'rows');
    
    % 4. 统计数量
    TP = size(TP_matrix, 1);       % 真正例：成功重构出的真实团
    num_P = size(P, 1);            % 重构出的总团数
    num_T = size(T, 1);            % 真实的总体团数
    
    FP = num_P - TP;               % 假正例：重构出但实际不存在的团 (误报)
    FN = num_T - TP;               % 假负例：真实存在但未重构出的团 (漏报)

    % 5. 计算指标
    precision = TP / num_P;
    recall = TP / num_T;

    % 避免分母为0的情况
    if (precision + recall) > 0
        f1_score = 2 * (precision * recall) / (precision + recall);
    else
        f1_score = 0;
    end
end