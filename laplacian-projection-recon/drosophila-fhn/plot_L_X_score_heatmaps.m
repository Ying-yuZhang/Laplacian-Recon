%% L 与 X_score 的热力图：隐藏对角线，共用配色和颜色条
% 使用方法：先运行 main，使工作区中存在 L 和 X_score，再运行本脚本。
% 不改动 L、X_score；不对重构结果取绝对值、截断或对称化。
% X_score 按工作区中的原始行列顺序显示。
% main 按行组织观测，回归系数对应 L.' 的索引方向。若需按接收节点/源节点
% 与 L 严格对齐，可将下方 X_heatmap 的赋值改为 full(double(X_score.'))。

%% 1. 检查数据，并仅在绘图副本中隐藏对角线
assert(exist('L','var') == 1 && exist('X_score','var') == 1, ...
    '请先运行 main.m，使工作区中存在 L 和 X_score。');
validateattributes(L, {'numeric'}, {'2d','square','real','finite'});
validateattributes(X_score, {'numeric'}, {'2d','square','real','finite'});
assert(isequal(size(L), size(X_score)), 'L 与 X_score 的尺寸必须相同。');
nNodes = size(L, 1);
assert(nNodes >= 2, '矩阵至少应包含两个节点。');

L_heatmap = full(double(L));
X_heatmap = full(double(X_score));
diagonalIndex = 1:nNodes+1:nNodes*nNodes;
L_heatmap(diagonalIndex) = NaN;
X_heatmap(diagonalIndex) = NaN;

%% 2. 统一颜色范围
% [] 表示根据两幅图所有非对角元素自动确定范围。
% 若需要固定范围，可改为例如 [0, 80]；超出该范围的数值会颜色饱和。
global_clim = [];

plotValues = [L_heatmap(~isnan(L_heatmap)); X_heatmap(~isnan(X_heatmap))];
if isempty(global_clim)
    global_clim = [min(0, min(plotValues)), max(0, max(plotValues))];
    if global_clim(1) == global_clim(2)
        global_clim = [0, 1];    % 两个矩阵非对角元素全为 0 时
    end
end
assert(numel(global_clim) == 2 && all(isfinite(global_clim)) && ...
    global_clim(1) < global_clim(2), 'global_clim 必须为递增的两个有限数。');

%% 3. 白 -> 深海军蓝 -> 赤陶色，与参考图保持相同渐变比例
c_white = [1.00, 1.00, 1.00];
c_blue  = [0.05, 0.25, 0.55];
c_harm  = [0.85, 0.40, 0.30];

% 深蓝位于色轴跨度的 1/3 处；色轴为 [0,6] 时恰好对应数值 2。
% 自动范围包含负值时，白色对应下限，而不是严格对应数值 0。
cmap = interp1([0, 1/3, 1], [c_white; c_blue; c_harm], ...
    linspace(0, 1, 256));

%% 4. 两个并排子图、一个共用颜色条
figL = figure('Color', 'w', 'Position', [100, 100,680, 310]);
layoutL = tiledlayout(figL, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

axL_true = nexttile(layoutL);
plot_matrix_heatmap(axL_true, L_heatmap, ...
    'Ground truth $L$', global_clim);

axL_recon = nexttile(layoutL);
plot_matrix_heatmap(axL_recon, X_heatmap, ...
    'Reconstructed $\widehat{L}$', global_clim);

colormap(figL, cmap);
cbL = colorbar(axL_recon);
cbL.Layout.Tile = 'east';
cbL.LineWidth = 1;
cbL.TickDirection = 'out';
cbL.FontSize = 11;
% cbL.Label.String = 'Coupling Weight';
cbL.Label.FontSize = 13;
cbL.Label.FontWeight = 'bold';

fprintf('共用色轴范围：[%.6g, %.6g]；两幅图均隐藏对角线。\n', global_clim);

%% 5. 如需保存图片，取消下面一行的注释
% exportgraphics(figL, 'L_X_score_heatmaps.png', 'Resolution', 300);

%% 辅助函数
function plot_matrix_heatmap(ax, M, titleText, colorLimits)
h = imagesc(ax, M);
h.AlphaData = ~isnan(M);    % 对角线透明，显示坐标区的白色背景
axis(ax, 'image');
clim(ax, colorLimits);
% title(ax, titleText, ...
%     'Interpreter', 'latex', ...
%     'FontSize', 15, ...
%     'FontWeight', 'bold');
xlabel(ax, 'Node Index', 'FontSize', 13, 'FontWeight', 'bold');
ylabel(ax, 'Node Index', 'FontSize', 13, 'FontWeight', 'bold');

n = size(M, 1);
nodeTicks = unique(round(linspace(1, n, min(6, n))));
set(ax, 'XTick', nodeTicks, 'YTick', nodeTicks, ...
    'TickDir', 'out', 'Box', 'on', 'LineWidth', 1, ...
    'FontSize', 11, 'Color', 'w', 'YDir', 'reverse');
end

