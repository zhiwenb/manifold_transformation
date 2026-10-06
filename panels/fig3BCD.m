[S, fig_png, fig_pdf] = did_group_tests( ...
  figure_input('code_manifold_combine/manifold_capacity/Results_manifold/var/analysis/M1.mat'), ...
  figure_input('code_manifold_combine/manifold_capacity/Results_manifold/var/analysis/M2.mat'), ...
  'save_dir', figure_output('fig3BCD'), ...
  'n_perm', 10000, 'B_boot', 10000);



function [S, fig_png, fig_pdf] = did_group_tests(M1_file, M2_file, varargin)
% Read two DiD result files containing table T, test each group and plot three panels.
% Metrics: alphaM_DiD tests < 0; RM_DiD and DM_DiD test > 0.

% Parameters
p = inputParser;
addParameter(p,'save_dir',pwd,@ischar);
addParameter(p,'n_perm',10000,@(x)isnumeric(x)&&x>10);
addParameter(p,'B_boot',10000,@(x)isnumeric(x)&&x>100);
parse(p,varargin{:});
save_dir = p.Results.save_dir;
n_perm  = p.Results.n_perm;
B_boot  = p.Results.B_boot;
if ~exist(save_dir,'dir'), mkdir(save_dir); end

% Read tables
T1 = load_table_from_did(M1_file);
T2 = load_table_from_did(M2_file);
std_names = {'alphaM_DiD','RM_DiD','DM_DiD'};
T1 = standardize_cols(T1, std_names);
T2 = standardize_cols(T2, std_names);

% Overall statistical tests
S = struct(); S.M1 = struct(); S.M2 = struct();
S.M1.alpha = one_sample_block(T1.alphaM_DiD, 'down', n_perm, B_boot);
S.M1.R     = one_sample_block(T1.RM_DiD,     'up',   n_perm, B_boot);
S.M1.D     = one_sample_block(T1.DM_DiD,     'up',   n_perm, B_boot);

S.M2.alpha = one_sample_block(T2.alphaM_DiD, 'down', n_perm, B_boot);
S.M2.R     = one_sample_block(T2.RM_DiD,     'up',   n_perm, B_boot);
S.M2.D     = one_sample_block(T2.DM_DiD,     'up',   n_perm, B_boot);

% Print a concise summary
fprintf('\n==== OVERALL summary ====\n');
pp('M1 α', S.M1.alpha);  pp('M1 R', S.M1.R);  pp('M1 D', S.M1.D);
pp('M2 α', S.M2.alpha);  pp('M2 R', S.M2.R);  pp('M2 D', S.M2.D);

% Plot three panels comparing M1 and M2
fig = figure('Color','w','Position',[200 200 1400 450]);
tiledlayout(1,3,'TileSpacing','loose','Padding','compact');

% Define colors
colors = [0.2 0.4 0.7; 0.8 0.3 0.3];  % M1 is blue; M2 is red.

% Alpha: expected direction < 0
nexttile; 
plot_one_metric_enhanced('Capacity (α) - DiD', ...
    [S.M1.alpha.mean, S.M2.alpha.mean], ...
    [S.M1.alpha.ci(1), S.M2.alpha.ci(1)], ...
    [S.M1.alpha.ci(2), S.M2.alpha.ci(2)], ...
    [S.M1.alpha.p_dir, S.M2.alpha.p_dir], 'down', colors);

% R: expected direction > 0
nexttile; 
plot_one_metric_enhanced('Radius (R) - DiD', ...
    [S.M1.R.mean, S.M2.R.mean], ...
    [S.M1.R.ci(1), S.M2.R.ci(1)], ...
    [S.M1.R.ci(2), S.M2.R.ci(2)], ...
    [S.M1.R.p_dir, S.M2.R.p_dir], 'up', colors);

% D: expected direction > 0
nexttile; 
plot_one_metric_enhanced('Dimension (D) - DiD', ...
    [S.M1.D.mean, S.M2.D.mean], ...
    [S.M1.D.ci(1), S.M2.D.ci(1)], ...
    [S.M1.D.ci(2), S.M2.D.ci(2)], ...
    [S.M1.D.p_dir, S.M2.D.p_dir], 'up', colors);

% Save outputs
ts = datestr(now,'yyyymmdd_HHMMSS');
fig_png = fullfile(save_dir, ['DiD_overall_M1_M2_' ts '.png']);
fig_pdf = fullfile(save_dir, ['DiD_overall_M1_M2_' ts '.pdf']);
exportgraphics(fig, fig_png, 'Resolution', 300);
exportgraphics(fig, fig_pdf, 'ContentType', 'vector');

fprintf('\nSaved figure:\n  %s\n  %s\n', fig_png, fig_pdf);

% Local functions
function T = load_table_from_did(fpath)
    Sload = load(fpath);
    if isfield(Sload,'T') && istable(Sload.T)
        T = Sload.T;
    elseif isfield(Sload,'DiD')
        try
            T = struct2table(Sload.DiD);
        catch
            error('文件 %s 中的 DiD 不是结构体数组，无法转为表。', fpath);
        end
    else
        % Fallback: find the first table at the top level.
        T = [];
        fns = fieldnames(Sload);
        for ii=1:numel(fns)
            if istable(Sload.(fns{ii})), T = Sload.(fns{ii}); break; end
        end
        if isempty(T), error('未在 %s 中找到表 T。', fpath); end
    end
end

function T = standardize_cols(T, keep)
    % Rename aliases to the standard column names.
    aliases = { ...
        'alphaM_DiD', {'alphaM_DiD','alpha_DiD','alpha','alphaM'}; ...
        'RM_DiD',     {'RM_DiD','R_DiD','R','RM'}; ...
        'DM_DiD',     {'DM_DiD','D_DiD','D','DM'}; ...
    };
    for i=1:size(aliases,1)
        tgt = aliases{i,1}; cands = aliases{i,2};
        if ~ismember(tgt, T.Properties.VariableNames)
            for j=1:numel(cands)
                if ismember(cands{j}, T.Properties.VariableNames)
                    T.Properties.VariableNames{cands{j}} = tgt; break;
                end
            end
        end
    end
    for i=1:numel(keep)
        if ~ismember(keep{i}, T.Properties.VariableNames)
            error('表中缺少所需列：%s', keep{i});
        end
    end
end

function out = one_sample_block(x, dir, n_perm, B_boot)
    % One-sample mean versus zero: sign-flip permutation and bootstrap CI; return the directional p value.
    x = x(:); x = x(isfinite(x));
    N = numel(x);
    if N==0
        out = struct('N',0,'mean',NaN,'ci',[NaN NaN],'p_up',NaN,'p_down',NaN,'p_dir',NaN);
        return;
    end
    mu = mean(x);
    % Permutation test: flip signs.
    null_means = zeros(n_perm,1);
    for i=1:n_perm
        sgn = (rand(N,1)>0.5)*2 - 1;   % ±1
        null_means(i) = mean(sgn.*x);
    end
    p_up   = (sum(null_means >= mu)+1)/(n_perm+1);
    p_down = (sum(null_means <= mu)+1)/(n_perm+1);
    % Directional p value
    if strcmpi(dir,'up'),   p_dir = p_up;   else, p_dir = p_down; end
    % bootstrap CI for mean
    mn_boot = zeros(B_boot,1);
    for b=1:B_boot, mn_boot(b) = mean(x(randi(N,N,1))); end
    ci = quantile(mn_boot,[0.025 0.975]);
    out = struct('N',N,'mean',mu,'ci',ci,'p_up',p_up,'p_down',p_down,'p_dir',p_dir);
end

function pp(tag, st)
    stars = p2stars(st.p_up, st.p_down, st.mean);
    fprintf('%-6s  N=%d  mean=%.4f  CI95=[%.4f, %.4f]  sig=%s\n', ...
        tag, st.N, st.mean, st.ci(1), st.ci(2), stars);
end

function stars = p2stars(p_up, p_down, mu)
    % Use the observed direction: p_up for mu >= 0, otherwise p_down.
    if mu>=0, p = p_up; else, p = p_down; end
    if p < 1e-3, stars='***';
    elseif p < 1e-2, stars='**';
    elseif p < 0.05, stars='*';
    else, stars='n.s.';
    end
end

function add_subplot_label(label)
    % Add a panel label such as A, B or C at the upper-left corner.
    ax = gca;
    % Get the axes position.
    xl = xlim(ax);
    yl = ylim(ax);
    
    % Position the label relative to the axes limits.
    text(xl(1) - 0.15*range(xl), yl(2), label, ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'Clipping', 'off');
end

function plot_one_metric_enhanced(title_str, means, ci_lo, ci_hi, p_vec_dir, expect_dir, colors)
    % Plot colored bars, error bars, significance stars and p values.
    
    % Draw bars.
    bh = bar(1:2, means, 0.65); 
    bh.FaceColor = 'flat';
    bh.CData = colors;
    bh.EdgeColor = 'none';
    hold on;
    
    % Add error bars with thicker lines.
    yneg = means - ci_lo; 
    ypos = ci_hi - means;
    errorbar(1:2, means, yneg, ypos, 'k', 'LineWidth', 2, ...
             'LineStyle','none', 'CapSize', 10);
    
    % Zero reference line
    yline(0, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
    
    % Set adaptive y limits before drawing annotations.
    yall = [ci_lo(:); ci_hi(:); 0];
    y_range = range(yall);
    y_min = min(yall) - 0.15*y_range;
    y_max = max(yall) + 0.30*y_range;  % Leave additional space above for annotations.
    ylim([y_min, y_max]);
    
    % Add stars and p values after setting ylim.
    for i = 1:2
        p = p_vec_dir(i);
        
        % Determine significance stars.
        if p < 0.001
            star = '***';
        elseif p < 0.01
            star = '**';
        elseif p < 0.05
            star = '*';
        else
            star = 'n.s.';
        end
        
        % Place stars above the error bars.
        xpos = i;
        star_ypos = ci_hi(i) + 0.2*y_range;
        
        % Draw stars in a larger font for visibility.
        text(xpos, star_ypos, star, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold', ...
             'FontSize', 16, ...
             'Color', 'k', ...
             'Clipping', 'off');  % Disable clipping.
        
        % Draw p values above the stars.
        % if p < 0.001
        %     p_text = 'p < 0.001';
        % else
        %     p_text = sprintf('p = %.3f', p);
        % end
        % 
        % text(xpos, star_ypos + 0.08*y_range, p_text, ...
        %      'HorizontalAlignment', 'center', ...
        %      'VerticalAlignment', 'bottom', ...
        %      'FontSize', 10, ...
        %      'Color', [0.3 0.3 0.3], ...
        %      'Clipping', 'off');  % Disable clipping.
    end
    
    % Axes settings
    set(gca, 'XTick', 1:2, 'XTickLabel', {'M1','M2'}, ...
             'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
    xlabel('Group', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('DiD Effect', 'FontSize', 15, 'FontWeight', 'bold');
    title(title_str, 'FontSize', 20, 'FontWeight', 'bold');
    
    % Grid and borders
    grid on;
    set(gca, 'GridAlpha', 0.15, 'GridLineStyle', ':');
    box on;
end

end