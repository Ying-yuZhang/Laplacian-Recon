function out = sol_int(alpha, X, int_max)
% 广义高阶结构整数推断 (General-Order Integer Inference)
% ==========================================================
% 输入:
%   alpha   : 1 x D 耦合强度向量 [alpha_1, alpha_2, ..., alpha_D]
%   X       : 1 x E 观测到的多阶耦合边权重 (只传入非零边)
%   int_max : 各阶最大结构数。可以是标量(全阶通用)，也可以是 1 x (D-1) 的向量
% 输出:
%   out     : D x E 的推断矩阵。第 d 行表示第 d 阶结构的推断整数值
% ==========================================================

    alpha = alpha(:)'; % 确保是 1 x D 的行向量
    D = length(alpha);
    % E = length(X);
    
    % 1. 动态生成各阶的最大上限
    % 如果传入的是标量，则所有高阶结构共用同一个上限；
    % 如果传入的是向量，则分别约束
    if isscalar(int_max)
        int_max_vec = repmat(int_max, 1, D-1);
    else
        int_max_vec = int_max(:)';
        if length(int_max_vec) ~= D-1
            error('int_max 必须是标量，或者长度为 D-1 的向量');
        end
    end
    
    % 2. 利用 ndgrid 动态生成全阶笛卡尔积
    grids = cell(1, D);
    grids{1} = [1 0]; % 第一阶 (基础连边) 永远是二元的：1 或 0
    
    % 第 2 到 D 阶是整数的：0 到 int_max_vec(d-1)
    for d = 2:D
        grids{d} = 0:int_max_vec(d-1);
    end
    
    % 动态生成高维网格
    [grids{:}]=ndgrid(grids{:}); 
    
    % 展平为 D x K 的组合矩阵 (K 为所有可能的组合总数)
    K = numel(grids{1});
    int_num = zeros(D, K);
    for d = 1:D
        gd = permute(grids{d}, D:-1:1); 
        int_num(d, :) = gd(:);
    end
    
    % 剔除全 0 组合 (因为传入的 X 都是有效边，必然包含至少一种结构)
    zero_idx = sum(int_num, 1) == 0;
    int_num(:, zero_idx) = [];
    
    % 3. 计算所有可能组合的理论期望值
    int_s = alpha * int_num; % 结果为 1 x K 的行向量
    
    % 4. 求最小值的索引
    [~, min_idx] = min(abs(X(:) - int_s), [], 2); 
    
    % 5. 映射回对应的整数组合
    out = int_num(:, min_idx);
end