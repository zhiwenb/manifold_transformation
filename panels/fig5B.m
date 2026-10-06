%% Final Clean Script: Session-level Orthogonality Analysis
% Keep original grouping and original angle definition
% Use more robust session-level statistics
%
% Output:
%   - one final large figure
%   - session-level summary statistics in command window
%
% Required external function:
%   - getPara(FR_by_category)

clear; clc; close all;

%% 1. Configuration
rootDir = figure_input('audit_20261005/old_global_followup/fr');

idx_G = 1;
idx_H = 2;
idx_R = 3;
idx_S = 4;
idx_T = 5;

analysis_groups = {
    struct('name', 'Polar',       'indices', [idx_R, idx_S, idx_T], 'p_x', 'ang_freq', 'p_y', 'radial_freq', 'p_type', 'linear'), ...
    struct('name', 'Grating',     'indices', [idx_G],               'p_x', 'ori_deg',  'p_y', 'freq_cpd',    'p_type', 'circular_x'), ...
    struct('name', 'Hyperbolic',  'indices', [idx_H],               'p_x', 'ori_deg',  'p_y', 'radial_freq', 'p_type', 'circular_x')
};

animal_sessions = struct();
animal_sessions.M1 = {'s1', 's4'};
animal_sessions.M2 = {'s6', 's7'};

special_session_config = struct();
special_session_config.s4 = struct('N_pre', 1, 'N_post', 2);

polar_allowed_sessions = {'s1', 's4', 's6', 's7'};
hyp_allowed_sessions   = {'s1', 's7'};

sessionDirs = dir(fullfile(rootDir, 's*'));
sessionDirs = sessionDirs([sessionDirs.isdir]);

%% 2. Collect session-level deltas
num_groups = length(analysis_groups);

session_delta_mat = [];   % [n_sessions x 3]
session_names = {};
animal_names = {};

for s = 1:length(sessionDirs)
    sessionName = sessionDirs(s).name;

    % Determine animal identity
    animal_name = '';
    animal_list = fieldnames(animal_sessions);
    for a = 1:length(animal_list)
        if any(strcmp(sessionName, animal_sessions.(animal_list{a})))
            animal_name = animal_list{a};
            break;
        end
    end

    if isempty(animal_name)
        fprintf('Skipping session %s (not in defined animal list)\n', sessionName);
        continue;
    end

    sessionPath = fullfile(rootDir, sessionName, 'global');
    matFiles = dir(fullfile(sessionPath, '*.mat'));
    if isempty(matFiles)
        fprintf('Skipping session %s (no .mat files found)\n', sessionName);
        continue;
    end

    fprintf('\nProcessing session: %s (%s)\n', sessionName, animal_name);

    Results = struct('date', {}, 'ortho', {});

    for f = 1:length(matFiles)
        current_filename = matFiles(f).name;
        filePath = fullfile(sessionPath, current_filename);

        S = load(filePath, 'FR_by_category');
        if ~isfield(S, 'FR_by_category')
            fprintf('  [SKIP] %s | FR_by_category not found\n', current_filename);
            continue;
        end
        FR_by_category = S.FR_by_category;

        [P, ~, ~] = getPara(FR_by_category);

        dateMatch = regexp(current_filename, '\d{6}', 'match');
        if isempty(dateMatch)
            fprintf('  [SKIP] %s | no date code found in filename\n', current_filename);
            continue;
        end
        current_date = dateMatch{1};

        daily_ortho = nan(1, num_groups);
        is_data_ok = true;

        for g = 1:num_groups
            group = analysis_groups{g};

            % inclusion rules
            if g == 1 && ~ismember(sessionName, polar_allowed_sessions)
                continue;
            end
            if g == 3 && ~ismember(sessionName, hyp_allowed_sessions)
                continue;
            end

            [X_combined, p_x_raw, p_y_raw] = stitch_data( ...
                FR_by_category, P, group.indices, group.p_x, group.p_y);

            if isempty(X_combined) || size(X_combined,1) < 2
                fprintf('  [SKIP] %s | %s | missing essential data\n', current_filename, group.name);
                is_data_ok = false;
                break;
            end

            [p_x_norm, p_y_norm] = normalize_params(p_x_raw, p_y_raw, group.p_type);

            if strcmp(group.p_type, 'linear')
                daily_ortho(g) = calculate_cartesian_orthogonality(X_combined, p_x_norm, p_y_norm);
            else
                daily_ortho(g) = calculate_cylinder_orthogonality_original(X_combined, p_x_raw, p_y_raw);
            end
        end

        if is_data_ok
            Results(end+1).date = current_date; %#ok<SAGROW>
            Results(end).ortho = daily_ortho;
            fprintf('  %s -> [Polar %.2f, Grating %.2f, Hyper %.2f]\n', ...
                current_date, daily_ortho(1), daily_ortho(2), daily_ortho(3));
        end
    end

    if isempty(Results)
        fprintf('  No valid files in session %s\n', sessionName);
        continue;
    end

    % Sort by date
    [~, sort_idx] = sort({Results.date});
    Results = Results(sort_idx);

    ortho_raw = [Results.ortho];
    ortho_matrix = reshape(ortho_raw, num_groups, length(Results))';
    N_total = size(ortho_matrix, 1);

    % Freeze the dates selected by the original analysis; intersect with available dates.
    frozen=struct();
    frozen.s1.early = {'041116','041416','041716'};
    frozen.s1.late = {'042116','042416','042516'};
    frozen.s4.early = {'031217'};
    frozen.s4.late = {'031517','032317'};
    frozen.s6.early = {'092217','092417','092717'};
    frozen.s6.late = {'092917','100117','100317'};
    frozen.s7.early = {'081218','081518','081818'};
    frozen.s7.late = {'082018','082218','082418'};
    early_mask=ismember({Results.date},frozen.(sessionName).early);
    late_mask=ismember({Results.date},frozen.(sessionName).late);
    assert(any(early_mask)&&any(late_mask),'Missing matched early or late stage');
    Y_pre=ortho_matrix(early_mask,:);Y_post=ortho_matrix(late_mask,:);
    audit_daily.(sessionName)=Results;
    audit_stages.(sessionName)=struct('early',{{Results(early_mask).date}},'late',{{Results(late_mask).date}});
    session_delta = mean(Y_post, 1, 'omitnan') - mean(Y_pre, 1, 'omitnan');

    % apply hyper inclusion rule
    if ~ismember(sessionName, hyp_allowed_sessions)
        session_delta(3) = NaN;
    end

    session_delta_mat = [session_delta_mat; session_delta]; %#ok<AGROW>
    session_names{end+1} = sessionName; %#ok<SAGROW>
    animal_names{end+1} = animal_name; %#ok<SAGROW>

    fprintf('  Session-level delta: [Polar %.2f, Grating %.2f, Hyper %.2f]\n', ...
        session_delta(1), session_delta(2), session_delta(3));
end

%% 3. Summary
fprintf('\n========================================\n');
fprintf('Session-level cross-session summary\n');
fprintf('========================================\n');

n_sessions = size(session_delta_mat, 1);
fprintf('Total valid sessions: %d\n\n', n_sessions);

if n_sessions < 2
    error('Not enough valid sessions for summary analysis.');
end

fprintf('%-10s %-8s %12s %12s %12s\n', 'Session', 'Animal', 'Polar', 'Grating', 'Hyper');
fprintf('%s\n', repmat('-', 1, 62));
for i = 1:n_sessions
    fprintf('%-10s %-8s %12.2f %12.2f %12.2f\n', ...
        session_names{i}, animal_names{i}, ...
        session_delta_mat(i,1), session_delta_mat(i,2), session_delta_mat(i,3));
end

%% 4. Statistics
mean_delta = nan(1, num_groups);
sem_delta  = nan(1, num_groups);
n_group    = zeros(1, num_groups);
ci_low     = nan(1, num_groups);
ci_high    = nan(1, num_groups);

for g = 1:num_groups
    d = session_delta_mat(:, g);
    d = d(~isnan(d));
    n_group(g) = numel(d);

    if ~isempty(d)
        mean_delta(g) = mean(d);
        if numel(d) > 1
            sem_delta(g) = std(d, 0, 1) / sqrt(numel(d));
        else
            sem_delta(g) = NaN;
        end
        [ci_low(g), ci_high(g)] = bootstrap_ci_session(d, 10000, 0.95);
    end
end

% paired differences vs Grating
valid_pg = ~isnan(session_delta_mat(:,1)) & ~isnan(session_delta_mat(:,2));
valid_hg = ~isnan(session_delta_mat(:,3)) & ~isnan(session_delta_mat(:,2));

diff_polar_grating = session_delta_mat(valid_pg,1) - session_delta_mat(valid_pg,2);
diff_hyper_grating = session_delta_mat(valid_hg,3) - session_delta_mat(valid_hg,2);

stats = struct();
stats.p_polar_ttest    = NaN;
stats.p_hyper_ttest    = NaN;
stats.t_polar          = NaN;
stats.t_hyper          = NaN;
stats.df_polar         = NaN;
stats.df_hyper         = NaN;
stats.ci_polar_diff    = [NaN NaN];
stats.ci_hyper_diff    = [NaN NaN];
stats.n_pos_polar      = sum(diff_polar_grating > 0);
stats.n_pos_hyper      = sum(diff_hyper_grating > 0);

if numel(diff_polar_grating) >= 2
    [~, stats.p_polar_ttest, ~, tstat_pg] = ttest(diff_polar_grating, 0, 'tail', 'right');
    stats.t_polar  = tstat_pg.tstat;
    stats.df_polar = tstat_pg.df;
    [ci_lo, ci_hi] = bootstrap_ci_session(diff_polar_grating, 10000, 0.95);
    stats.ci_polar_diff = [ci_lo, ci_hi];
end

if numel(diff_hyper_grating) >= 2
    [~, stats.p_hyper_ttest, ~, tstat_hg] = ttest(diff_hyper_grating, 0, 'tail', 'right');
    stats.t_hyper  = tstat_hg.tstat;
    stats.df_hyper = tstat_hg.df;
    [ci_lo, ci_hi] = bootstrap_ci_session(diff_hyper_grating, 10000, 0.95);
    stats.ci_hyper_diff = [ci_lo, ci_hi];
end

fprintf('\n=== Group-level session summary ===\n');
for g = 1:num_groups
    if n_group(g) > 0
        fprintf('%s: n=%d, mean=%.2f, 95%% CI=[%.2f, %.2f]\n', ...
            analysis_groups{g}.name, n_group(g), mean_delta(g), ci_low(g), ci_high(g));
    else
        fprintf('%s: no data\n', analysis_groups{g}.name);
    end
end

fprintf('\n=== Paired comparison vs Grating (one-sample t-test, one-tailed) ===\n');
if ~isempty(diff_polar_grating)
    fprintf('Polar - Grating: n=%d, mean=%.2f, positive=%d/%d, 95%% CI=[%.2f, %.2f], t(%d)=%.3f, p(one-tail)=%.4f\n', ...
        numel(diff_polar_grating), mean(diff_polar_grating), ...
        stats.n_pos_polar, numel(diff_polar_grating), ...
        stats.ci_polar_diff(1), stats.ci_polar_diff(2), ...
        stats.df_polar, stats.t_polar, stats.p_polar_ttest);
else
    fprintf('Polar - Grating: no paired data\n');
end

if ~isempty(diff_hyper_grating)
    fprintf('Hyper - Grating: n=%d, mean=%.2f, positive=%d/%d, 95%% CI=[%.2f, %.2f], t(%d)=%.3f, p(one-tail)=%.4f\n', ...
        numel(diff_hyper_grating), mean(diff_hyper_grating), ...
        stats.n_pos_hyper, numel(diff_hyper_grating), ...
        stats.ci_hyper_diff(1), stats.ci_hyper_diff(2), ...
        stats.df_hyper, stats.t_hyper, stats.p_hyper_ttest);
else
    fprintf('Hyper - Grating: no paired data\n');
end

%% 5. Final plot
AllSessionStats = struct();
AllSessionStats.delta        = session_delta_mat;
AllSessionStats.session      = session_names;
AllSessionStats.animal       = animal_names;
AllSessionStats.mean_delta   = mean_delta;
AllSessionStats.sem_delta    = sem_delta;
AllSessionStats.n_group      = n_group;
AllSessionStats.ci_low       = ci_low;
AllSessionStats.ci_high      = ci_high;
AllSessionStats.diff_pg      = diff_polar_grating;
AllSessionStats.diff_hg      = diff_hyper_grating;
AllSessionStats.stats        = stats;
audit_result=struct();
audit_result.daily=audit_daily;
audit_result.stages=audit_stages;
audit_result.AllSessionStats=AllSessionStats;
fid=fopen(figure_input('audit_20261005/old_global_followup/results/old/fig5B.json'),'w');fprintf(fid,'%s',jsonencode(audit_result,PrettyPrint=true));fclose(fid);
fprintf('AUDIT_SAVED fig5B\n');

plot_final_session_summary(AllSessionStats);

function plot_final_session_summary(AllSessionStats)

    all_deltas   = AllSessionStats.delta;
    all_animals  = AllSessionStats.animal;
    mean_delta   = AllSessionStats.mean_delta;
    sem_delta    = AllSessionStats.sem_delta;
    n_group      = AllSessionStats.n_group;
    stats        = AllSessionStats.stats;

    color_polar   = [0.16, 0.43, 0.85];
    color_grating = [0.45, 0.45, 0.45];
    color_hyper   = [0.90, 0.45, 0.15];
    colors = [color_polar; color_grating; color_hyper];

    animal_markers = struct('M1', 'o', 'M2', 's');
    animal_colors  = struct('M1', [0.20, 0.20, 0.20], 'M2', [0.70, 0.70, 0.70]);

    figure('Name', 'Final Orthogonality Summary (Session-level)', 'Color', 'w');
    set(gcf, 'Units', 'inches', 'Position', [1.5, 1.5, 8.8, 6.8]);
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 8.8, 6.8]);

    x_pos = 1:3;
    valid_bar = ~isnan(mean_delta);

    hBar = bar(x_pos(valid_bar), mean_delta(valid_bar), ...
        'EdgeColor', 'k', 'LineWidth', 1.2, 'BarWidth', 0.72);
    hBar.FaceColor = 'flat';
    hBar.CData = colors(valid_bar, :);
    hold on;

    errorbar(x_pos(valid_bar), mean_delta(valid_bar), sem_delta(valid_bar), ...
        'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 10);

    % one point per session
    h_animals = [];
    animal_names_for_legend = {};

    for i = 1:size(all_deltas,1)
        this_row = all_deltas(i, :);
        valid = ~isnan(this_row);
        if ~any(valid)
            continue;
        end

        x_jitter = x_pos(valid) + (rand(1, sum(valid)) - 0.5) * 0.14;

        if isfield(animal_colors, all_animals{i})
            pt_color = animal_colors.(all_animals{i});
            pt_marker = animal_markers.(all_animals{i});
        else
            pt_color = [0.5, 0.5, 0.5];
            pt_marker = 'o';
        end

        h_scatter = scatter(x_jitter, this_row(valid), 55, pt_color, pt_marker, ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.55, 'MarkerFaceAlpha', 0.9);

        if ~any(strcmp(animal_names_for_legend, all_animals{i}))
            h_animals = [h_animals, h_scatter]; %#ok<AGROW>
            animal_names_for_legend{end+1} = all_animals{i}; %#ok<AGROW>
        end
    end

    valid_vals = all_deltas(~isnan(all_deltas));
    y_max = max(valid_vals);
    y_min = min(valid_vals);
    y_range = max(y_max - y_min, 1);

    current_y = y_max + 0.12 * y_range;

    if ~isnan(stats.p_polar_ttest) && stats.p_polar_ttest < 0.05
        plot_sig_line(2, 1, current_y, stats.p_polar_ttest, y_range);
        current_y = current_y + 0.18 * y_range;
    end

    if ~isnan(stats.p_hyper_ttest) && stats.p_hyper_ttest < 0.05
        plot_sig_line(2, 3, current_y, stats.p_hyper_ttest, y_range);
    end

    yline(0, '-', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1);

    ax = gca;
    ax.FontSize = 16;
    ax.LineWidth = 1.2;
    ax.TickDir = 'out';
    ax.TickLength = [0.02, 0.02];
    ax.Box = 'off';

ylabel('\Delta Orthogonality (°)', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'XTick', 1:3, ...
         'XTickLabel', {'Polar', 'Grating', 'Hyperbolic'}, ...
         'FontSize', 20, ...
         'FontWeight', 'bold');

    if ~isempty(h_animals)
        legend(h_animals, animal_names_for_legend, ...
            'Location', 'northwest', 'FontSize', 13, 'Box', 'off');
    end

    ylim([min(y_min - 0.15*y_range, -5), current_y + 0.20*y_range]);

    print(gcf, 'Final_Orthogonality_SessionLevel_originalAngle', '-dpng', '-r600');
    print(gcf, 'Final_Orthogonality_SessionLevel_originalAngle', '-dpdf');
    print(gcf, 'Final_Orthogonality_SessionLevel_originalAngle', '-dsvg');

    fprintf('\nSaved figure:\n');
    fprintf('  Final_Orthogonality_SessionLevel_originalAngle.png\n');
    fprintf('  Final_Orthogonality_SessionLevel_originalAngle.pdf\n');
    fprintf('  Final_Orthogonality_SessionLevel_originalAngle.svg\n');
end

function plot_sig_line(x1, x2, y, p, y_range)
    line_h = 0.02 * y_range;
    plot([x1 x1 x2 x2], [y-line_h y y y-line_h], '-k', 'LineWidth', 1.0);

    if p < 0.001
        p_text = '***  p<0.001';
    elseif p < 0.01
        p_text = sprintf('**  p=%.3f', p);
    else
        p_text = sprintf('*');
    end

    text(mean([x1, x2]), y + 0.03*y_range, p_text, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 13, ...
        'FontWeight', 'normal');
end

%% ===========================
%% Statistical Helper Functions
%% ===========================
function [ci_lower, ci_upper] = bootstrap_ci_session(data, N_boot, confidence)
    data = data(:);
    n = length(data);

    if n == 0
        ci_lower = NaN;
        ci_upper = NaN;
        return;
    end

    boot_means = zeros(N_boot, 1);
    for i = 1:N_boot
        idx = randi(n, n, 1);
        boot_sample = data(idx);
        boot_means(i) = mean(boot_sample);
    end

    alpha = 1 - confidence;
    ci_lower = prctile(boot_means, alpha/2 * 100);
    ci_upper = prctile(boot_means, (1 - alpha/2) * 100);
end

%% ===========================
%% Geometry Helper Functions
%% ===========================
function [X_combined, p_x_raw, p_y_raw] = stitch_data(FR_by_category, P, indices, p_x_name, p_y_name)

    X_combined = [];
    p_x_raw = [];
    p_y_raw = [];

    table_names = {'G', 'H', 'R', 'S', 'T'};

    for c = indices
        if c > length(FR_by_category)
            continue;
        end

        X_c = mean(FR_by_category{c}, 3)';  % [conditions x neurons]
        if isempty(X_c) || size(X_c, 1) < 2
            X_combined = [];
            p_x_raw = [];
            p_y_raw = [];
            return;
        end

        T_name = table_names{c};
        T_params = P.(T_name);

        if ~ismember(p_x_name, T_params.Properties.VariableNames) || ...
           ~ismember(p_y_name, T_params.Properties.VariableNames)
            X_combined = [];
            p_x_raw = [];
            p_y_raw = [];
            return;
        end

        X_combined = [X_combined; X_c];
        p_x_raw = [p_x_raw; T_params.(p_x_name)];
        p_y_raw = [p_y_raw; T_params.(p_y_name)];
    end
end

function [p_x_norm, p_y_norm] = normalize_params(p_x_raw, p_y_raw, p_type) %#ok<INUSD>
    p_x_norm = (p_x_raw - min(p_x_raw)) / (range(p_x_raw) + eps);
    p_y_norm = (p_y_raw - min(p_y_raw)) / (range(p_y_raw) + eps);
end

function ortho_val = calculate_cartesian_orthogonality(X_combined, p_x_norm, p_y_norm)

    X_centered = X_combined - mean(X_combined, 1);
    [~, score, ~, explained] = pca(X_centered);

    cum_var = cumsum(explained);
    k_90 = find(cum_var >= 90, 1);
    k_final = max(k_90, 3);
    k_final = min(k_final, size(score, 2));

    X_pca_z = zscore(score(:, 1:k_final));
    B_x = pinv(X_pca_z) * p_x_norm;
    B_y = pinv(X_pca_z) * p_y_norm;

    if norm(B_x) < eps || norm(B_y) < eps
        ortho_val = NaN;
        return;
    end

    cos_theta = dot(B_x, B_y) / (norm(B_x) * norm(B_y));
    cos_theta = max(min(cos_theta, 1), -1);
    ortho_val = acosd(abs(cos_theta));
end

function ortho_val = calculate_cylinder_orthogonality_original(X_combined, p_x_raw, p_y_raw)

    ori_rad = deg2rad(p_x_raw);
    cos_comp = cos(2 * ori_rad);
    sin_comp = sin(2 * ori_rad);
    p_linear = p_y_raw;

    X_centered = X_combined - mean(X_combined, 1);
    [~, score, ~, explained] = pca(X_centered);

    cum_var = cumsum(explained);
    k_90 = find(cum_var >= 90, 1);
    k_final = max(k_90, 3);
    k_final = min(k_final, size(score, 2));

    X_pca_z = zscore(score(:, 1:k_final));
    B_F = pinv(X_pca_z) * zscore(p_linear);
    B_C = pinv(X_pca_z) * cos_comp;
    B_S = pinv(X_pca_z) * sin_comp;

    if norm(B_F) < eps || norm(B_C) < eps || norm(B_S) < eps
        ortho_val = NaN;
        return;
    end

    cos_FC = dot(B_F, B_C) / (norm(B_F) * norm(B_C));
    cos_FS = dot(B_F, B_S) / (norm(B_F) * norm(B_S));

    max_cos = max(abs(cos_FC), abs(cos_FS));
    max_cos = max(min(max_cos, 1), -1);

    ortho_val = acosd(max_cos);
end