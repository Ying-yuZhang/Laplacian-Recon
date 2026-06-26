function [EdgeList, TriangleList] = convert_to_list(A, B)
    % A: N x N 邻接矩阵
    % B: N x N x N 三超边张量
    % 返回:
    % EdgeList: 边列表，每行是一个边 (i, j)
    % TriangleList: 三角形列表，每行是一个三角形 (i, j, k)，且满足 i < j < k

    % 从邻接矩阵 A 构造 EdgeList
    [row, col] = find(tril(A, -1));  % 获取下三角部分的边
    EdgeList = [col, row];  % 边列表
    N = length(B);
    TriangleList = [];
    % 从三超边张量 B 构造 TriangleList， i < j < k
    for i = 1:N
        for j = i+1:N % 保证 j > i
            for k = j+1:N % 保证 k > j
                % 检查 B(i, j, k) 是否为 1（表示存在一个三超边）
                if B(i, j, k) == 1
                    % 对三角形进行排序，以保证顺序不影响最终结果
                    triangle = [i, j, k];
                    TriangleList = [TriangleList; triangle];
                end
            end
        end
    end

end
