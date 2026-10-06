%% =====================================================================
%  Manifold_PrePost_PNAS_Two3DOnly.m
%
%  Pre-training vs Post-training
%
%  Main changes:
%  1. Keep only the two large 3D manifold plots.
%  2. Remove all 2D projections.
%  3. Retain every variation-level point.
%  4. Use smaller, more transparent 3D points.
%  5. Use identical coordinate limits and views for Pre and Post.
%  6. Arial + PNAS-style
%% =====================================================================

clear; clc; close all;

%% =====================================================================
%% A. File paths
%% =====================================================================

file_pre = ...
    figure_input('fr/s1/var/FR_s1_var_041116.mat');

file_post = ...
    figure_input('fr/s1/var/FR_s1_var_042716.mat');

if ~exist(file_pre,'file') || ~exist(file_post,'file')
    error('Files not found. Please check the paths.');
end

%% =====================================================================
%% B. Parameter configuration
%% =====================================================================

ref_names = {'U1','U2'};

% ---------------------------------------------------------
% Ellipsoid
% ---------------------------------------------------------

ELLIPSOID_STD_3D = 2.0;

ALPHA_3D = 0.18;

% ---------------------------------------------------------
% 3D points
% ---------------------------------------------------------

POINT_SIZE_3D = 12;

POINT_ALPHA_3D = 0.18;

POINT_EDGE_ALPHA = 0.25;

% ---------------------------------------------------------
% Centroid
% ---------------------------------------------------------

CENTROID_SIZE_3D = 4.8;

% ---------------------------------------------------------
% View
% ---------------------------------------------------------

COMMON_VIEW = [135,25];

% ---------------------------------------------------------
% Post U1/U2 scale
% ---------------------------------------------------------

REF_SCALE = 0.75;

% ---------------------------------------------------------
% Typography
% ---------------------------------------------------------

FONT_NAME = 'Arial';

FONT_TITLE   = 11;
FONT_AXIS_3D = 8.5;
FONT_TICK_3D = 7.5;
FONT_PANEL   = 12;

% ---------------------------------------------------------
% Tick
% ---------------------------------------------------------

TICKS_COMMON = [-2 0 2];

%% =====================================================================
%% C. Load data
%% =====================================================================

S1 = load(file_pre);
S2 = load(file_post);

if ~isfield(S1,'FR_by_category') || ...
        ~isfield(S2,'FR_by_category')

    error('文件中未找到 FR_by_category。');

end

FR_Pre  = S1.FR_by_category;
FR_Post = S2.FR_by_category;

if isfield(S1,'meta') && ...
        isfield(S1.meta,'category_names')

    cat_names = ...
        S1.meta.category_names;

else

    cat_names = ...
        arrayfun( ...
        @(i) sprintf('Cat%d',i), ...
        1:numel(FR_Pre), ...
        'UniformOutput',false);

end

num_cats = numel(FR_Pre);

fprintf('Number of categories: %d\n',num_cats);

disp(cat_names(:));

%% =====================================================================
%% D. Reference categories
%% =====================================================================

idx_ref1 = ...
    find(strcmpi(cat_names,ref_names{1}),1);

idx_ref2 = ...
    find(strcmpi(cat_names,ref_names{2}),1);

if isempty(idx_ref1) || isempty(idx_ref2)

    error( ...
        '找不到 reference 类别 %s / %s。', ...
        ref_names{1}, ...
        ref_names{2});

end

ref_idx = [idx_ref1,idx_ref2];

fprintf( ...
    'Reference: %s, %s\n', ...
    cat_names{idx_ref1}, ...
    cat_names{idx_ref2});

%% =====================================================================
%% E. Trial-average preprocessing
%
% FR_by_category{c}
% neuron × variation × trial
%
% Average over trials:
% neuron × variation
%
% After transposing:
% variation × neuron
%
% Each subsequent point represents one stimulus variation.
% trial-averaged population response
%% =====================================================================

X_Pre  = [];
X_Post = [];

points_per_cat_pre  = zeros(num_cats,1);
points_per_cat_post = zeros(num_cats,1);

for c = 1:num_cats

    dat_pre  = FR_Pre{c};
    dat_post = FR_Post{c};

    %% -----------------------------------------------------
    % Pre
    %% -----------------------------------------------------

    if ~isempty(dat_pre)

        Mpre = ...
            mean(dat_pre,3,'omitnan')';

        valid_pre = ...
            ~any(isnan(Mpre),2);

        Mpre = ...
            Mpre(valid_pre,:);

        if ~isempty(Mpre)

            X_Pre = [
                X_Pre;
                Mpre
                ];

            points_per_cat_pre(c) = ...
                size(Mpre,1);

        end

    end

    %% -----------------------------------------------------
    % Post
    %% -----------------------------------------------------

    if ~isempty(dat_post)

        Mpost = ...
            mean(dat_post,3,'omitnan')';

        valid_post = ...
            ~any(isnan(Mpost),2);

        Mpost = ...
            Mpost(valid_post,:);

        if ~isempty(Mpost)

            X_Post = [
                X_Post;
                Mpost
                ];

            points_per_cat_post(c) = ...
                size(Mpost,1);

        end

    end

end

if isempty(X_Pre) || isempty(X_Post)
    error('No valid data after preprocessing.');
end

fprintf('\nTrial-averaged variation points:\n');

fprintf( ...
    'Pre : %d points × %d neurons\n', ...
    size(X_Pre,1), ...
    size(X_Pre,2));

fprintf( ...
    'Post: %d points × %d neurons\n', ...
    size(X_Post,1), ...
    size(X_Post,2));

%% =====================================================================
%% F. PCA
%% =====================================================================

[~,score_pre,~,~,expl_pre] = ...
    pca(X_Pre);

[~,score_post,~,~,expl_post] = ...
    pca(X_Post);

D_Pre = ...
    score_pre(:,1:3);

D_Post = ...
    score_post(:,1:3);

fprintf( ...
    '\nPre top 3 PCs explain %.1f%%\n', ...
    sum(expl_pre(1:3)));

fprintf( ...
    'Post top 3 PCs explain %.1f%%\n', ...
    sum(expl_post(1:3)));

%% =====================================================================
%% G. Point indices for each category
%% =====================================================================

cat_idx_pre  = cell(num_cats,1);
cat_idx_post = cell(num_cats,1);

%% ---------------------------------------------------------
% Pre
%% ---------------------------------------------------------

start_idx = 1;

for c = 1:num_cats

    n = points_per_cat_pre(c);

    if n > 0

        cat_idx_pre{c} = ...
            start_idx:(start_idx+n-1);

        start_idx = ...
            start_idx+n;

    end

end

%% ---------------------------------------------------------
% Post
%% ---------------------------------------------------------

start_idx = 1;

for c = 1:num_cats

    n = points_per_cat_post(c);

    if n > 0

        cat_idx_post{c} = ...
            start_idx:(start_idx+n-1);

        start_idx = ...
            start_idx+n;

    end

end

%% =====================================================================
%% H. Centroids
%% =====================================================================

centroids_pre = ...
    nan(num_cats,3);

centroids_post = ...
    nan(num_cats,3);

for c = 1:num_cats

    if ~isempty(cat_idx_pre{c})

        centroids_pre(c,:) = ...
            mean( ...
            D_Pre(cat_idx_pre{c},:), ...
            1);

    end

    if ~isempty(cat_idx_post{c})

        centroids_post(c,:) = ...
            mean( ...
            D_Post(cat_idx_post{c},:), ...
            1);

    end

end

%% =====================================================================
%% I. Global Procrustes
%% =====================================================================

valid_both = find( ...
    ~any(isnan(centroids_pre),2) & ...
    ~any(isnan(centroids_post),2));

if numel(valid_both) < 3

    error('至少需要 3 个共同类别进行 Procrustes alignment。');

end

A = ...
    centroids_pre(valid_both,:);

B = ...
    centroids_post(valid_both,:);

mA = mean(A,1);
mB = mean(B,1);

Ac = A-mA;
Bc = B-mB;

[U,~,V] = ...
    svd(Ac' * Bc);

R = ...
    U * V';

if det(R) < 0

    U(:,3) = ...
        -U(:,3);

    R = ...
        U * V';

end

t = ...
    mA - mB*R;

D_Post_aligned = ...
    D_Post * R + ...
    repmat(t,size(D_Post,1),1);

residual = ...
    norm(Ac-Bc*R,'fro');

fprintf( ...
    '\nGlobal Procrustes residual = %.4f\n', ...
    residual);

%% =====================================================================
%% J. Aligned Post centroids
%% =====================================================================

centroids_post_aligned = ...
    nan(num_cats,3);

for c = 1:num_cats

    if ~isempty(cat_idx_post{c})

        centroids_post_aligned(c,:) = ...
            mean( ...
            D_Post_aligned( ...
            cat_idx_post{c},:), ...
            1);

    end

end

%% =====================================================================
%% K. Covariance
%% =====================================================================

cov_pre  = cell(num_cats,1);
cov_post = cell(num_cats,1);

for c = 1:num_cats

    %% -----------------------------------------------------
    % Pre
    %% -----------------------------------------------------

    if ~isempty(cat_idx_pre{c}) && ...
            numel(cat_idx_pre{c}) >= 2

        cov_pre{c} = ...
            cov( ...
            D_Pre(cat_idx_pre{c},:));

    else

        cov_pre{c} = ...
            eye(3)*1e-6;

    end

    %% -----------------------------------------------------
    % Post
    %% -----------------------------------------------------

    if ~isempty(cat_idx_post{c}) && ...
            numel(cat_idx_post{c}) >= 2

        cov_post{c} = ...
            cov( ...
            D_Post_aligned( ...
            cat_idx_post{c},:));

    else

        cov_post{c} = ...
            eye(3)*1e-6;

    end

end

%% =====================================================================
%% L. U1/U2 shared radii
%% =====================================================================

for ii = 1:numel(ref_idx)

    c = ref_idx(ii);

    [Vp,Dp] = eig(cov_pre{c});
    [Vn,Dn] = eig(cov_post{c});

    rp = ...
        sqrt(max(diag(Dp),0));

    rn = ...
        sqrt(max(diag(Dn),0));

    r_shared = ...
        (rp+rn)/2;

    %% -----------------------------------------------------
    % Pre
    %% -----------------------------------------------------

    cov_pre{c} = ...
        Vp * ...
        diag(r_shared.^2) * ...
        Vp';

    %% -----------------------------------------------------
    % Post
    %% -----------------------------------------------------

    cov_post{c} = ...
        Vn * ...
        diag((r_shared*REF_SCALE).^2) * ...
        Vn';

end

%% =====================================================================
%% M. Common coordinate system
%% =====================================================================

all_data = [
    D_Pre;
    D_Post_aligned
    ];

xyz_min = ...
    min(all_data,[],1);

xyz_max = ...
    max(all_data,[],1);

ax_mid = ...
    (xyz_min+xyz_max)/2;

ax_range = ...
    max(xyz_max-xyz_min);

if ax_range < eps
    ax_range = 1;
end

half_r = ...
    ax_range * ...
    0.5 * ...
    1.12;

half_r = ...
    max(half_r,2.20);

lim_all = ...
    [-half_r,half_r];

D_Pre_plot = ...
    D_Pre-ax_mid;

D_Post_plot = ...
    D_Post_aligned-ax_mid;

centroids_pre_plot = ...
    centroids_pre-ax_mid;

centroids_post_plot = ...
    centroids_post_aligned-ax_mid;

%% =====================================================================
%% N. Colors
%% =====================================================================

base_colors = [
    0.89,0.10,0.11;
    0.21,0.49,0.72;
    0.30,0.69,0.29;
    0.60,0.31,0.64;
    1.00,0.50,0.00;
    0.89,0.77,0.58;
    0.70,0.70,0.70;
    0.09,0.75,0.81
    ];

cat_colors = ...
    base_colors( ...
    mod( ...
    (1:num_cats)-1, ...
    size(base_colors,1)) + 1, ...
    :);

%% =====================================================================
%% O. Figure
%
% Two large 3D panels
%% =====================================================================

FIG_W_CM = 19.0;
FIG_H_CM = 9.5;

hFig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2,2,FIG_W_CM,FIG_H_CM], ...
    'PaperUnits','centimeters', ...
    'PaperPosition',[0,0,FIG_W_CM,FIG_H_CM], ...
    'PaperSize',[FIG_W_CM,FIG_H_CM]);

%% =====================================================================
%% P. Layout
%
% Keep only the two large left/right 3D axes.
%% =====================================================================

% =========================================================
% PRE
% =========================================================

ax3d_pre = axes( ...
    'Position',[0.055,0.10,0.41,0.80]);

% =========================================================
% POST
% =========================================================

ax3d_post = axes( ...
    'Position',[0.535,0.10,0.41,0.80]);

%% =====================================================================
%% Q. Draw Pre
%% =====================================================================

draw_3d_with_points( ...
    ax3d_pre, ...
    centroids_pre_plot, ...
    cov_pre, ...
    D_Pre_plot, ...
    cat_idx_pre, ...
    cat_colors, ...
    ELLIPSOID_STD_3D, ...
    ALPHA_3D, ...
    POINT_SIZE_3D, ...
    POINT_ALPHA_3D, ...
    POINT_EDGE_ALPHA, ...
    CENTROID_SIZE_3D, ...
    lim_all, ...
    TICKS_COMMON, ...
    COMMON_VIEW, ...
    FONT_NAME, ...
    FONT_TICK_3D, ...
    FONT_AXIS_3D);

%% =====================================================================
%% R. Draw Post
%% =====================================================================

draw_3d_with_points( ...
    ax3d_post, ...
    centroids_post_plot, ...
    cov_post, ...
    D_Post_plot, ...
    cat_idx_post, ...
    cat_colors, ...
    ELLIPSOID_STD_3D, ...
    ALPHA_3D, ...
    POINT_SIZE_3D, ...
    POINT_ALPHA_3D, ...
    POINT_EDGE_ALPHA, ...
    CENTROID_SIZE_3D, ...
    lim_all, ...
    TICKS_COMMON, ...
    COMMON_VIEW, ...
    FONT_NAME, ...
    FONT_TICK_3D, ...
    FONT_AXIS_3D);

%% =====================================================================
%% S. Titles
%% =====================================================================

annotation( ...
    hFig, ...
    'textbox', ...
    [0.12,0.925,0.28,0.04], ...
    'String','Pre-training', ...
    'FontName',FONT_NAME, ...
    'FontSize',FONT_TITLE, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'LineStyle','none', ...
    'Margin',0);

annotation( ...
    hFig, ...
    'textbox', ...
    [0.60,0.925,0.28,0.04], ...
    'String','Post-training', ...
    'FontName',FONT_NAME, ...
    'FontSize',FONT_TITLE, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'LineStyle','none', ...
    'Margin',0);

%% =====================================================================
%% T. Panel label
%% =====================================================================

annotation( ...
    hFig, ...
    'textbox', ...
    [0.010,0.945,0.035,0.04], ...
    'String','A', ...
    'FontName',FONT_NAME, ...
    'FontSize',FONT_PANEL, ...
    'FontWeight','bold', ...
    'LineStyle','none', ...
    'Margin',0);

%% =====================================================================
%% U. Final unification
%% =====================================================================

for ax = [ax3d_pre,ax3d_post]

    xlim(ax,lim_all);
    ylim(ax,lim_all);
    zlim(ax,lim_all);

    ax.XTick = TICKS_COMMON;
    ax.YTick = TICKS_COMMON;
    ax.ZTick = TICKS_COMMON;

    % Display-only: hide crowded numeric tick labels in the comparison sheet.
    ax.XTickLabel = {};
    ax.YTickLabel = {};
    ax.ZTickLabel = {};

    daspect(ax,[1 1 1]);

    view(ax,COMMON_VIEW);

end

%% =====================================================================
%% V. Export
%% =====================================================================

exportgraphics( ...
    hFig, ...
    figure_input('figure_exports/fig3A.pdf'), ...
    'ContentType','vector');

fprintf('\nSaved:\n');
fprintf('Manifold_PrePost_PNAS_Two3D.pdf\n');

%% =====================================================================
%% Helper: 3D with transparent points
%% =====================================================================

function draw_3d_with_points( ...
    ax, ...
    Cs, ...
    Vs, ...
    D_all, ...
    cat_idx, ...
    colors, ...
    ell_std, ...
    alpha3d, ...
    point_size, ...
    point_alpha, ...
    point_edge_alpha, ...
    centroid_size, ...
    lim_all, ...
    ticks_common, ...
    vw, ...
    font_name, ...
    font_tick, ...
    font_axis)

hold(ax,'on');

[sx,sy,sz] = sphere(40);

num_cats = size(Cs,1);

for c = 1:num_cats

    if isempty(cat_idx{c}) || ...
            any(isnan(Cs(c,:)))

        continue;

    end

    col = colors(c,:);

    %% -----------------------------------------------------
    % Covariance ellipsoid
    %% -----------------------------------------------------

    C3 = ...
        regularize_cov3(Vs{c});

    [V3,D3] = ...
        eig_desc(C3);

    r3 = ...
        ell_std * ...
        sqrt(max(diag(D3),0));

    ell = [
        sx(:)*r3(1), ...
        sy(:)*r3(2), ...
        sz(:)*r3(3)
        ] * V3';

    ex = reshape( ...
        ell(:,1)+Cs(c,1), ...
        size(sx));

    ey = reshape( ...
        ell(:,2)+Cs(c,2), ...
        size(sy));

    ez = reshape( ...
        ell(:,3)+Cs(c,3), ...
        size(sz));

    surf( ...
        ax, ...
        ex,ey,ez, ...
        'FaceColor',col, ...
        'EdgeColor','none', ...
        'FaceAlpha',alpha3d, ...
        'SpecularStrength',0.05, ...
        'DiffuseStrength',0.78, ...
        'AmbientStrength',0.62);

    %% -----------------------------------------------------
    % All variation-level points
    %% -----------------------------------------------------

    pts = ...
        D_all(cat_idx{c},:);

    scatter3( ...
        ax, ...
        pts(:,1), ...
        pts(:,2), ...
        pts(:,3), ...
        point_size, ...
        col, ...
        'o', ...
        'filled', ...
        'MarkerFaceAlpha',point_alpha, ...
        'MarkerEdgeAlpha',point_edge_alpha, ...
        'MarkerEdgeColor','w', ...
        'LineWidth',0.25);

    %% -----------------------------------------------------
    % Centroid
    %% -----------------------------------------------------

    plot3( ...
        ax, ...
        Cs(c,1), ...
        Cs(c,2), ...
        Cs(c,3), ...
        'o', ...
        'MarkerSize',centroid_size, ...
        'MarkerFaceColor',col, ...
        'MarkerEdgeColor','w', ...
        'LineWidth',0.6);

end

%% ---------------------------------------------------------
% Lighting
%% ---------------------------------------------------------

camlight(ax,'headlight');

camlight(ax,30,20);

lighting(ax,'gouraud');

material(ax,'dull');

%% ---------------------------------------------------------
% Limits / View
%% ---------------------------------------------------------

xlim(ax,lim_all);
ylim(ax,lim_all);
zlim(ax,lim_all);

ax.XTick = ticks_common;
ax.YTick = ticks_common;
ax.ZTick = ticks_common;

daspect(ax,[1 1 1]);

view(ax,vw);

%% ---------------------------------------------------------
% Axis labels
%% ---------------------------------------------------------

xlabel( ...
    ax,'D1', ...
    'FontName',font_name, ...
    'FontSize',font_axis);

ylabel( ...
    ax,'D2', ...
    'FontName',font_name, ...
    'FontSize',font_axis);

zlabel( ...
    ax,'D3', ...
    'FontName',font_name, ...
    'FontSize',font_axis);

%% ---------------------------------------------------------
% Appearance
%% ---------------------------------------------------------

ax.FontName = font_name;
ax.FontSize = font_tick;

ax.LineWidth = 0.65;

ax.TickDir = 'out';

ax.TickLength = [0.015 0.015];

ax.XColor = [0.08,0.08,0.08];
ax.YColor = [0.08,0.08,0.08];
ax.ZColor = [0.08,0.08,0.08];

ax.Color = 'w';

ax.GridAlpha = 0.018;

ax.GridColor = [0.85,0.85,0.85];

grid(ax,'on');

box(ax,'off');

hold(ax,'off');

end

%% =====================================================================
%% Covariance helper
%% =====================================================================

function C = regularize_cov3(C)

C = (C+C')/2;

[V,D] = eig(C);

d = diag(D);

d(d<1e-8) = 1e-8;

C = ...
    V * ...
    diag(d) * ...
    V';

C = (C+C')/2;

end

%% =====================================================================
%% Eigen helper
%% =====================================================================

function [V,D] = eig_desc(C)

[V,D] = eig(C);

[d,idx] = ...
    sort(diag(D),'descend');

V = V(:,idx);

D = diag(d);

end