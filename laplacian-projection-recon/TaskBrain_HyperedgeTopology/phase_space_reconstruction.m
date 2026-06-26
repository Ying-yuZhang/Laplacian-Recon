function X = phase_space_reconstruction(data, m, tau)
% PHASE_SPACE_RECONSTRUCTION 将多维时间序列重构为相空间嵌入矩阵
%
% 根据时间延迟嵌入定理 (Takens' Embedding Theorem)，将 N 个节点的
% 1 维时间序列，扩展为 N*m 维的相空间状态矩阵。
%
% 输入参数:
%   data : N x L 的矩阵，N 为节点数(例如 82)，L 为时间序列长度(例如 621)
%   m    : 嵌入维度 (Embedding dimension)，例如 3
%   tau  : 延迟步长 (Time delay step)，例如 2
%
% 输出参数:
%   X    : T x (N*m) 的相空间重构矩阵，其中 T = L - (m-1)*tau
%
% 调用示例:
%   X = phase_space_reconstruction(ERP_data, 3, 2);

    % 获取输入数据的维度
    [N, L] = size(data);

    % 严谨的边界条件检查：确保延迟长度不会超过时间序列总长度
    max_delay = (m - 1) * tau;
    if max_delay >= L
        error('重构失败：最大延迟步数 (m-1)*tau (%d) 超过了时间序列的总长度 L (%d)！请减小参数 m 或 tau。', max_delay, L);
    end

    % 1. 确定有效的时间起点和长度
    t_start = 1 + max_delay;  
    t_indices = t_start:L;        % 有效的时间点索引
    T = length(t_indices);        % 最终矩阵的时间点数 T

    % 2. 预分配最终的 X 矩阵 [T x (N*m)]，提升内存分配效率
    X = zeros(T, N * m);

    % 3. 按照公式构造矩阵
    for i = 1:N
        % 取出第 i 个节点的完整时间序列
        V_i = data(i, :);
        
        % 构造单个节点 i 的延迟嵌入子矩阵 [T x m]
        x_i = zeros(T, m);
        for k = 1:m
            delay = (k - 1) * tau; % 第 k 维度的延迟量
            x_i(:, k) = V_i(t_indices - delay)'; % 截取历史状态并转置为列
        end
        
        % 计算该节点在总矩阵 X 中的列索引范围，并填入数据
        col_start = (i - 1) * m + 1;
        col_end   = i * m;
        X(:, col_start:col_end) = x_i;
    end
end