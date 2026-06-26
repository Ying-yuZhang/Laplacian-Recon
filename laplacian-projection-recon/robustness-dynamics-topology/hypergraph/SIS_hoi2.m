function f = SIS_hoi2(t, x, EdgeList, TriangleList, gamma, beta, I)
% x(i): 节点 i 的感染概率/状态
% gamma(i): 恢复率 γ_i
% beta(1): pairwise 感染率 β_1
% beta(2): triangle 感染率 β_2
% I: 外部输入（可为 0）

N = length(x);

coup_pair   = zeros(N,1);
coup_tri    = zeros(N,1);

%% -------- pairwise infection term: sum_j a_ij x_j ----------
for ii = 1:size(EdgeList,1)
    i1 = EdgeList(ii,1);
    i2 = EdgeList(ii,2);

    coup_pair(i1) = coup_pair(i1) + x(i2);
    coup_pair(i2) = coup_pair(i2) + x(i1);
end

%% -------- triangle infection term: sum_{j,k} b_ijk x_j x_k ----------
for ii = 1:size(TriangleList,1)
    i1 = TriangleList(ii,1);
    i2 = TriangleList(ii,2);
    i3 = TriangleList(ii,3);

    % 每个三超边对三个节点都贡献 x_j x_k 项
    coup_tri(i1) = coup_tri(i1) + x(i2) * x(i3);
    coup_tri(i2) = coup_tri(i2) + x(i1) * x(i3);
    coup_tri(i3) = coup_tri(i3) + x(i1) * x(i2);
end

%% -------- Final SIS-HOI dynamics ----------
f = -gamma .* x ...
    + beta(1) * (1 - x) .* coup_pair ...
    + beta(2) * (1 - x) .* coup_tri ...
    + I;

end
