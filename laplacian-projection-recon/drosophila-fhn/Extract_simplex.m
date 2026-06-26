function out = Extract_simplex(X_score, ite_num, init_alpha, delta, cut, order)
% 解耦出邻接矩阵 A_estimate 和 任意阶单纯形列表
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入： X_score     多阶Laplacian矩阵
%        ite_num     迭代次数
%        init_alpha  估计的初始耦合强度 (长度为 order 的向量)
%        delta       = 10^-3 误差
%        cut         判断邻接矩阵有没有边的阈值
%        order       最高重构阶数 (d = 2, 3, ..., order)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = zeros(order, ite_num + 1);
alpha(:, 1) = init_alpha(:);
N = size(X_score, 1);
X_score = (X_score + X_score') * 0.5;
A_estimate = X_score > cut;

for ite = 1:ite_num
    XX = X_score;
    extracted_simplices = cell(1, order); % 用元胞数组存储每一阶的单纯形列表
    K_tmp_all = cell(1, order);           % 存储每一阶向底层图投影的边权重贡献

    for d = 2:order
        XX = XX - alpha(d-1, ite);
        X_num = XX;

        % ---------------------------------------------------------
        % 仅在底层有效边构成的稀疏图中寻找大小为 d+1 的完全子图
        Adj_init = sparse(X_num > delta);
        Adj_init = Adj_init | Adj_init';
        Adj_init = Adj_init - diag(diag(Adj_init));

        % 调用专门的极速找团辅助函数
        candidates = get_cliques(Adj_init, d + 1);
        num_cands = size(candidates, 1);

        picked_list = zeros(0, d + 1); % 记录成功剥离的 d-单纯形
        valid_cands = true(num_cands, 1);

        if num_cands == 0
            extracted_simplices{d} = picked_list;
            K_tmp_all{d} = sparse(N, N);
            break; % 该阶已无有效结构
        end

        % 预先计算 d+1 个节点中包含的所有边对组合索引
        edge_subs = nchoosek(1:(d+1), 2);
        num_edges_per_simplex = size(edge_subs, 1);

        % ---------------------------------------------------------
        % 贪心剥离
        while true
            valid_edges = X_num > delta;

            % 快速验证候选单纯形是否仍然有效 (其包含的所有边必须都 > delta)
            if any(valid_cands)
                is_active = valid_cands;
                for e_idx = 1:num_edges_per_simplex
                    n1 = candidates(:, edge_subs(e_idx, 1));
                    n2 = candidates(:, edge_subs(e_idx, 2));
                    % 利用线性索引极速判断
                    lin_idx = sub2ind([N, N], n1, n2);
                    is_active = is_active & valid_edges(lin_idx);
                end
            else
                is_active = false;
            end

            if ~any(is_active)
                break; % 无有效候选，结束当前阶剥离
            end

            active_idx = find(is_active);
            act_cands = candidates(active_idx, :);
            num_act = size(act_cands, 1);

            % 极速构建当前活跃单纯形的边集，用于计算 K_S
            all_i = zeros(num_act * num_edges_per_simplex * 2, 1);
            all_j = zeros(num_act * num_edges_per_simplex * 2, 1);
            ptr = 1;
            for e_idx = 1:num_edges_per_simplex
                n1 = act_cands(:, edge_subs(e_idx, 1));
                n2 = act_cands(:, edge_subs(e_idx, 2));
                len = length(n1);
                all_i(ptr : ptr+len-1) = n1;
                all_j(ptr : ptr+len-1) = n2;
                ptr = ptr + len;
                all_i(ptr : ptr+len-1) = n2;
                all_j(ptr : ptr+len-1) = n1; % 保证对称
                ptr = ptr + len;
            end

            % K_S：计算当前网络中每条边被多少个活跃 d-单纯形所共享
            K_S = sparse(all_i, all_j, 1, N, N);

            % 边概率 pp = 权重 / 共享次数
            [bs_i, bs_j, bs_v] = find(K_S);
            bs_idx = sub2ind([N, N], bs_i, bs_j);
            pp = sparse(N, N);
            pp(bs_idx) = X_num(bs_idx) ./ bs_v;

            % 计算每个活跃单纯形的总概率得分 (内部边概率求和)
            scores = zeros(num_act, 1);
            for e_idx = 1:num_edges_per_simplex
                n1 = act_cands(:, edge_subs(e_idx, 1));
                n2 = act_cands(:, edge_subs(e_idx, 2));
                % 取出该边对应的概率进累加行
                scores = scores + full(pp(sub2ind([N, N], n1, n2)));
            end

            % 贪心选择：剥离得分最大的单纯形
            [~, best_loc] = max(scores);
            best_cand_idx = active_idx(best_loc);
            best_cand = candidates(best_cand_idx, :);

            picked_list = [picked_list; best_cand];
            valid_cands(best_cand_idx) = false; % 移出候选池

            % 更新底层边权重 X_num
            subtract_val = alpha(d, ite);
            for e_idx = 1:num_edges_per_simplex
                n1 = best_cand(edge_subs(e_idx, 1));
                n2 = best_cand(edge_subs(e_idx, 2));
                X_num(n1, n2) = X_num(n1, n2) - subtract_val;
                X_num(n2, n1) = X_num(n2, n1) - subtract_val;
            end
        end

        % =========================================================
        % 最小二乘法更新耦合强度
        % =========================================================
        % 根据剥离列表快速重建 K_tmp 矩阵，直接统计每条边被选中的单纯形覆盖了多少次！
        K_tmp = sparse(N, N);
        for p = 1:size(picked_list, 1)
            cand = picked_list(p, :);
            for e_idx = 1:num_edges_per_simplex
                n1 = cand(edge_subs(e_idx, 1));
                n2 = cand(edge_subs(e_idx, 2));
                K_tmp(n1, n2) = K_tmp(n1, n2) + 1;
                K_tmp(n2, n1) = K_tmp(n2, n1) + 1;
            end
        end

        extracted_simplices{d} = picked_list;
        K_tmp_all{d} = K_tmp;
    end

    % ---------------------------------------------------------
    % 步骤 3: 最小二乘法更新耦合强度
    % ---------------------------------------------------------
    C = zeros(N * N, order);
    C(:, 1) = A_estimate(:);
    for d = 2:order
        if ~isempty(K_tmp_all{d})
            Kd_tmp = K_tmp_all{d};
            C(:, d) = Kd_tmp(:);
        end
    end

    x = X_score(:);
    C_id = find(abs(sum(C, 2)) > 0);
    C_filt = C(C_id, :);
    x_filt = x(C_id);

    %al = (C_filt' * C_filt + 1e-6 * eye(size(C_filt, 2))) \ (C_filt' * x_filt);
    al = (C_filt' * C_filt ) \ (C_filt' * x_filt);
    if any(al <= 0)
        disp('alpha为负');
        disp(al);
        break;
    end
    alpha(:, ite + 1) = al;
end

out.alpha = alpha;
out.A_estimate = A_estimate;
% =========================================================
% 步骤 4: 构造兼容输出
% =========================================================
A2_estimate = zeros(N, N, N);
picked_list=extracted_simplices{2};
for p = 1:size(picked_list, 1)
    i = picked_list(p, 1); j = picked_list(p, 2); k = picked_list(p, 3);
    % 对称写入 6 个排列
    A2_estimate(i, j, k) = 1; A2_estimate(i, k, j) = 1;
    A2_estimate(k, i, j) = 1; A2_estimate(k, j, i) = 1;
    A2_estimate(j, k, i) = 1; A2_estimate(j, i, k) = 1;
end
out.A2_estimate = A2_estimate;

% 使用动态字段名输出，输出标准的 S x (d+1) 节点表
for d = 2:order
    field_name = sprintf('A%d_simplex_list', d);
    out.(field_name) = extracted_simplices{d};
end
disp('Finish the calculation');
end

% =========================================================================
% 极速找团辅助函数：在稀疏图中寻找所有的 k-完全子图 (k-cliques)
% =========================================================================
function cliques = get_cliques(Adj, k)
[r, c] = find(triu(Adj, 1));
cliques = [r, c]; % 初始化 2-cliques (边)

% 逐阶升维：用 (d-1)-单纯形 扩展寻找 d-单纯形
for step = 3:k
    num_prev = size(cliques, 1);
    if num_prev == 0
        cliques = zeros(0, step);
        return;
    end

    new_cliques = zeros(0, step);
    for i = 1:num_prev
        curr_clique = cliques(i, :);
        last_node = curr_clique(end);

        % 候选节点必须同时是当前团内所有节点的邻居
        common_neighbors = Adj(curr_clique(1), :);
        for n_idx = 2:length(curr_clique)
            common_neighbors = common_neighbors & Adj(curr_clique(n_idx), :);
        end

        candidates = find(common_neighbors);
        % 强制约束 k_node > j_node，从物理上杜绝排列组合产生的冗余！
        candidates = candidates(candidates > last_node);

        if ~isempty(candidates)
            new_cliques = [new_cliques; repmat(curr_clique, length(candidates), 1), candidates(:)];
        end
    end
    cliques = new_cliques;
end
end