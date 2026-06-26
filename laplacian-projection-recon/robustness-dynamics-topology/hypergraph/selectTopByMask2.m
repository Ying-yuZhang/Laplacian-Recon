function AA2_bin = selectTopByMask2(AA, AA2)
% selectTopByMask 保留 AA2 中权重最大的 top-e 个元素为 1，其余置 0
%
% 输入：
%   AA   - H×N 二值矩阵 (真实结构，只包含 0/1)
%   AA2  - H×N 实数矩阵 (推断或权重矩阵)
%
% 输出：
%   AA2_bin - H×N 二值矩阵 (仅保留 top-e 最大位置)
%
% 示例：
%   AA  = [1 0 1; 0 0 1];
%   AA2 = [0.2 0.8 0.5; 0.9 0.1 0.3];
%   AA2_bin = selectTopByMask(AA, AA2);



% 1. 计算AA中1的总数 e
e = nnz(AA(:));

% 2. 将AA2展平为向量
vals = AA2(:);

% 3. 找出 top-e 最大值的索引
[~, idx] = maxk(vals, e);

% 4. 构造结果矩阵
AA2_bin = zeros(size(AA2));
AA2_bin(idx) = 1;

end