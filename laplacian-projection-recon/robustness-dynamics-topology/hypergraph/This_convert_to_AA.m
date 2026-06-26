function AA = This_convert_to_AA(EdgeList, TriangleList,N)
%number of unknowns for each node
H=(N-1)+(N-1)*(N-2)/2;
AA=zeros(H,N); %"true" matrix

% 初始化邻接矩阵
A = zeros(N, N);

% 填充邻接矩阵
for e = 1:size(EdgeList, 1)
    i = EdgeList(e, 1);
    j = EdgeList(e, 2);
    A(i, j) = EdgeList(e, 3);
end

% === 预构建三元组字符串索引，用于快速匹配 ===
triStr = strcat(string(TriangleList(:,1)), '_', ...
                string(TriangleList(:,2)), '_', ...
                string(TriangleList(:,3)));

for i = 1:N
    % ---- 2-超边部分 ----
    AA(1:N-1, i) = A(i, [1:i-1 i+1:N]);

    % ---- 3-超边部分 ----
    itemp = 1;
    for ii1 = [1:i-1 i+1:N]
        for jj1 = ii1+1:N
            if jj1 == i
                continue;
            end

            % 生成当前三元组字符串（排序以避免顺序影响）
            thisStr = sprintf('%d_%d_%d', [i, ii1, jj1]);

            % 判断是否在 triStr 内
            th = find(strcmp(triStr, thisStr), 1);  % 返回匹配的行索引

            if ~isempty(th)
                % 若找到匹配，赋值为该三元组的权重
                AA(N-1+itemp, i) = TriangleList(th, 4);
            end

            itemp = itemp + 1;
        end
    end
end

end

