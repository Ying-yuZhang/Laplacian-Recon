function F1 = clc_F1(AA, AA2)
% 计算二值矩阵预测结果的 F1 分数
%
% 输入:
%   AA   - H×N 二值矩阵 (真实标签, 仅包含 0/1)
%   AA2  - H×N 二值矩阵 (预测结果, 仅包含 0/1)
%
% 输出:
%   F1   - F1 分数 (范围 [0,1])
%
% 计算公式:
%   Precision = TP / (TP + FP)
%   Recall    = TP / (TP + FN)
%   F1        = 2 * (Precision * Recall) / (Precision + Recall)

    % 检查输入大小是否一致
    if ~isequal(size(AA), size(AA2))
        error('AA 和 AA2 必须具有相同的大小。');
    end

    % 转为逻辑类型以确保计算正确
    AA  = logical(AA);
    AA2 = logical(AA2);

    % 计算混淆矩阵的四个量
    TP = sum(AA(:) & AA2(:));          % 真正例
    FP = sum(~AA(:) & AA2(:));         % 假正例
    FN = sum(AA(:) & ~AA2(:));         % 假负例

    % 精确率与召回率
    if TP + FP == 0
        Precision = 0;
    else
        Precision = TP / (TP + FP);
    end

    if TP + FN == 0
        Recall = 0;
    else
        Recall = TP / (TP + FN);
    end

    % F1 分数
    if Precision + Recall == 0
        F1 = 0;
    else
        F1 = 2 * Precision * Recall / (Precision + Recall);
    end
end
