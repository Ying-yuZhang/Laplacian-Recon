function alpha = Init_alpha_SC(X_score,order,cut)
%对耦合强度初始的估计
%输入：X_score多阶拉普拉斯矩阵
%order 结构最高阶数
%cut X_score中小于cut的值被视为0
%输出：初始化估计的耦合强度[0.4 0.6 ...]一阶二阶...
% % 示例数据
% X = randn(1, 1000); % 随机生成1000个数据点
X=X_score(find(X_score>cut));

% 定义直方图的边界和获取频率
edges = 0:0.05:20; % 自定义边界，细粒度
[counts, edges] = histcounts(X, 'BinEdges', edges);

% 找到峰值（局部最大值）
[peak_values , peak_indices] = findpeaks(counts);
[~, sort_idx] = sort(peak_values, 'descend'); % 按峰值大小排序

% 提取前 n 个峰的位置
n = order; % 设置提取的峰值数量
if length(peak_indices) >= n
    % 获取前 n 个峰的位置
    peak_positions = edges(peak_indices(sort_idx(1:n))+1);
    for i=n:-1:2
        if peak_positions(i)>peak_positions(i-1)
            peak_positions(i)=peak_positions(i)-peak_positions(i-1);
        end
    end
else
    % 如果峰数量少于 n，返回全部峰的位置
    peak_positions = edges(peak_indices(sort_idx)+1);
end

alpha=peak_positions;
% disp(alpha);
end


