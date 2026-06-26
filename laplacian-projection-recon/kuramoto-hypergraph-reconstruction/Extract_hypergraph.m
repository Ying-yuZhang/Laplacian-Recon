function out = Extract_hypergraph(X_score, ite_num, init_alpha, delta, cut, order)
% 解耦出邻接矩阵 A_estimate 和 一般阶超边列表 (任意 D 阶)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入： X_score     多阶Laplacian矩阵
%        ite_num     迭代次数
%        init_alpha  初始耦合强度向量 (长度为 order)
%        delta       = 10^-3 误差
%        cut         判断邻接矩阵有没有边的阈值
%        order       要重构的最高阶数 (例如 order=2 对应 3-超边, order=3 对应 4-超边)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = size(X_score, 1);
    X_score = (X_score + X_score') * 0.5;
    
    alpha = zeros(order, ite_num + 1);
    alpha(:, 1) = init_alpha(:); 
    
    % 提取下三角部分的索引，按列优先顺序展平为列向量
    low_tri = tril(true(N, N), -1);
    xi = X_score(low_tri)';
    int_max = N - 2; % 超边数量的最大搜索范围
    X_nzero = xi > cut;
    
    Ai_hat = zeros(order, length(xi));
    
    for ite = 1:ite_num
        % ---------------------------------------------------------
        % 步骤 0: 整数推断 
        % ---------------------------------------------------------
        Ai_hat(:, X_nzero) = sol_int(alpha(:, ite), xi(X_nzero), int_max);
        
        A_estimate = zeros(N);
        A_estimate(low_tri) = Ai_hat(1, :)';
        A_estimate = A_estimate + A_estimate';
        
        extracted_hyperedges = cell(1, order); % 存储各阶剥离出的超边列表
        K_tmp_all = cell(1, order);            % 存储各阶对边权重的贡献覆盖
        
        % ---------------------------------------------------------
        % 步骤 1 & 2: 逐阶贪心剥离 (d = 2 到 order)
        % ---------------------------------------------------------
        for d = 2:order
            X_d = zeros(N);
            X_d(low_tri) = Ai_hat(d, :)';
            X_d = X_d + X_d';
            
            X_num = X_d; % 当前 d 阶网络中每条边上期望的交互数
            
            % 基于 X_num 快速构建底层骨架的稀疏矩阵
            Adj_init = sparse(X_num > delta);
            Adj_init = Adj_init | Adj_init';
            Adj_init = Adj_init - diag(diag(Adj_init));
            
            % 极速提取所有候选的 (d+1)-节点 团
            candidates = get_cliques(Adj_init, d + 1);
            num_cands = size(candidates, 1);
            
            picked_list = zeros(0, d + 1);
            valid_cands = true(num_cands, 1);
            
            % 预生成超边内部的所有边对索引组合
            edge_subs = nchoosek(1:(d+1), 2);
            num_edges_per_he = size(edge_subs, 1);
            
            if num_cands == 0
                extracted_hyperedges{d} = picked_list;
                K_tmp_all{d} = sparse(N, N);
                continue; % 当前阶数无结构，直接跳入下一阶
            end
            
            while true
                valid_edges = X_num > delta;
                
                % 快速验证候选池中的超边是否依然成立 (所有内部边均 > delta)
                if any(valid_cands)
                    is_active = valid_cands;
                    for e_idx = 1:num_edges_per_he
                        lin_idx = sub2ind([N, N], candidates(:, edge_subs(e_idx, 1)), candidates(:, edge_subs(e_idx, 2)));
                        is_active = is_active & valid_edges(lin_idx);
                    end
                else
                    is_active = false;
                end
                
                if ~any(is_active)
                    break; % 找不到满足条件的超边，结束当前阶剥离
                end
                
                active_idx = find(is_active);
                act_cands = candidates(active_idx, :);
                num_act = size(act_cands, 1);
                
                % 构建有效边索引，用于极速统计 B_S
                all_i = zeros(num_act * num_edges_per_he * 2, 1);
                all_j = zeros(num_act * num_edges_per_he * 2, 1);
                ptr = 1;
                for e_idx = 1:num_edges_per_he
                    n1 = act_cands(:, edge_subs(e_idx, 1));
                    n2 = act_cands(:, edge_subs(e_idx, 2));
                    len = length(n1);
                    all_i(ptr : ptr+len-1) = n1;
                    all_j(ptr : ptr+len-1) = n2;
                    ptr = ptr + len;
                    all_i(ptr : ptr+len-1) = n2;
                    all_j(ptr : ptr+len-1) = n1;
                    ptr = ptr + len;
                end
                
                B_S = sparse(all_i, all_j, 1, N, N);
                
                % 计算边贡献概率 pp = X_num / B_S
                [bs_i, bs_j, bs_v] = find(B_S);
                bs_idx = sub2ind([N, N], bs_i, bs_j);
                pp = sparse(N, N);
                pp(bs_idx) = X_num(bs_idx) ./ bs_v;
                
                % 向量化计算所有活跃超边的总得分 (内部所有边概率之和)
                scores = zeros(num_act, 1);
                for e_idx = 1:num_edges_per_he
                    n1 = act_cands(:, edge_subs(e_idx, 1));
                    n2 = act_cands(:, edge_subs(e_idx, 2));
                    scores = scores + full(pp(sub2ind([N, N], n1, n2)));
                end
                
                % 选取得分最大的超边
                [~, best_loc] = max(scores);
                best_cand_idx = active_idx(best_loc);
                best_cand = candidates(best_cand_idx, :);
                
                picked_list = [picked_list; best_cand];
                valid_cands(best_cand_idx) = false;
                
                % 从网络中剥离该超边
                for e_idx = 1:num_edges_per_he
                    n1 = best_cand(edge_subs(e_idx, 1));
                    n2 = best_cand(edge_subs(e_idx, 2));
                    X_num(n1, n2) = X_num(n1, n2) - 1;
                    X_num(n2, n1) = X_num(n2, n1) - 1;
                end
            end
            
            % 利用剥离列表极速重建该阶的 K_tmp
            K_tmp = sparse(N, N);
            for p = 1:size(picked_list, 1)
                cand = picked_list(p, :);
                for e_idx = 1:num_edges_per_he
                    n1 = cand(edge_subs(e_idx, 1));
                    n2 = cand(edge_subs(e_idx, 2));
                    K_tmp(n1, n2) = K_tmp(n1, n2) + 1;
                    K_tmp(n2, n1) = K_tmp(n2, n1) + 1;
                end
            end
            
            extracted_hyperedges{d} = picked_list;
            K_tmp_all{d} = K_tmp;
        end

        % ---------------------------------------------------------
        % 步骤 3: 最小二乘法更新耦合强度 alpha
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
        
        al = (C_filt' * C_filt) \ (C_filt' * x_filt);
        
        if any(al <= 0)
            disp('alpha为负');
            disp(al);
            break;
        end
        alpha(:, ite + 1) = al;
    end
    
    % ---------------------------------------------------------
    % 步骤 4: 动态结构输出
    % ---------------------------------------------------------
    out.alpha = alpha;
    out.A_estimate = A_estimate;

    % =========================================================
    % 步骤 5: 构造兼容输出
    % =========================================================
    A2_estimate = zeros(N, N, N);
    picked_list=extracted_hyperedges{2};
    for p = 1:size(picked_list, 1)
        i = picked_list(p, 1); j = picked_list(p, 2); k = picked_list(p, 3);
        % 对称写入 6 个排列
        A2_estimate(i, j, k) = 1; A2_estimate(i, k, j) = 1;
        A2_estimate(k, i, j) = 1; A2_estimate(k, j, i) = 1;
        A2_estimate(j, k, i) = 1; A2_estimate(j, i, k) = 1;
    end
    out.A2_estimate = A2_estimate;
    
    % 输出形式为 out.A2_hyperedge_list, out.A3_hyperedge_list 等
    for d = 3:order
        field_name = sprintf('A%d_hyperedge_list', d);
        out.(field_name) = extracted_hyperedges{d};
    end
end

% =========================================================================
% 极速找团辅助函数 (与单纯复形逻辑通用)
% =========================================================================
function cliques = get_cliques(Adj, k)
    [r, c] = find(triu(Adj, 1));
    cliques = [r, c]; 
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
            common_neighbors = Adj(curr_clique(1), :);
            for n_idx = 2:length(curr_clique)
                common_neighbors = common_neighbors & Adj(curr_clique(n_idx), :);
            end
            candidates = find(common_neighbors);
            candidates = candidates(candidates > last_node);
            if ~isempty(candidates)
                new_cliques = [new_cliques; repmat(curr_clique, length(candidates), 1), candidates(:)];
            end
        end
        cliques = new_cliques;
    end
end