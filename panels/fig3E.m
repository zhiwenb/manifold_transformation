%% =====================================================================
%  Batch analysis — Early vs Late, per-pair lines colored by monkey
%  M1 = s1, s4  |  M2 = s5, s6, s7
% =====================================================================

clear; clc; close all;

%% -------------------- 1. Path setup --------------------
result_root = figure_input('code_manifold_combine/Decoding_Figure2/Results_Decoding_SVM_var');

save_dir = figure_output('fig3E');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% -------------------- 2. Define stage groups --------------------
stage_groups = struct();

stage_groups.s1.early  = {'041116','041216'};
stage_groups.s1.middle = {'041816'};
stage_groups.s1.late   = {'042916'};

stage_groups.s4.early  = {'031617'};
stage_groups.s4.middle = {''};
stage_groups.s4.late   = {'031917'};

stage_groups.s5.early  = {'071217'};
stage_groups.s5.middle = {'071617'};
stage_groups.s5.late   = {'071917','072217'};

stage_groups.s6.early  = {'092817'};
stage_groups.s6.middle = {'093017'};
stage_groups.s6.late   = {'100217'};

stage_groups.s7.early  = {'081318'};
stage_groups.s7.middle = {'081718','082018'};
stage_groups.s7.late   = {'082418'};

all_sessions = fieldnames(stage_groups);

% -------- Monkey group assignment --------
M1_sessions = {'s1','s4'};
M2_sessions = {'s5','s6','s7'};

col_M1 = [0.20 0.55 0.82];   % Blue tones identify M1.
col_M2 = [0.85 0.33 0.10];   % Orange tones identify M2.

%% -------------------- 3. Find files --------------------
files_all = dir(fullfile(result_root, '**', 'Decoder_Acc_*.mat'));
if isempty(files_all)
    error('No Decoder_Acc_*.mat files found under: %s', result_root);
end
is_variation = arrayfun(@(x) contains(x.name, '_variation_'), files_all);
files = files_all(is_variation);
if isempty(files)
    error('No variation decoding files found.');
end
fprintf('Found %d variation result files.\n', numel(files));

%% -------------------- 4. Collect per-pair values --------------------
session_data = struct();
for s = 1:numel(all_sessions)
    sess = all_sessions{s};
    session_data.(sess).early_rows = {};
    session_data.(sess).late_rows  = {};
end

for f = 1:numel(files)
    fileName = files(f).name;
    filePath = fullfile(files(f).folder, files(f).name);

    tok_sess = regexp(filePath, '[\\/](s\d+)[\\/]', 'tokens', 'once');
    if isempty(tok_sess), continue; end
    sess = tok_sess{1};
    if ~isfield(stage_groups, sess), continue; end

    tok_date = regexp(fileName, '(\d{6})', 'tokens', 'once');
    if isempty(tok_date), continue; end
    date_code = tok_date{1};

    stage = '';
    if     any(strcmp(date_code, stage_groups.(sess).early)),  stage = 'early';
    elseif any(strcmp(date_code, stage_groups.(sess).middle)), stage = 'middle';
    elseif any(strcmp(date_code, stage_groups.(sess).late)),   stage = 'late';
    else, continue;
    end

    try
        S = load(filePath, 'accuracy_matrix');
    catch, continue;
    end
    if ~isfield(S,'accuracy_matrix') || isempty(S.accuracy_matrix), continue; end

    A = S.accuracy_matrix;
    mask_ut   = triu(true(size(A)), 1);
    pair_vals = A(mask_ut);
    pair_vals = pair_vals(:)';

    if strcmp(stage, 'early')
        session_data.(sess).early_rows{end+1} = pair_vals;
    elseif strcmp(stage, 'late')
        session_data.(sess).late_rows{end+1}  = pair_vals;
    end
end

%% -------------------- 5. Build paired vectors --------------------
paired_early       = [];
paired_late        = [];
pair_session_label = {};
pair_monkey_label  = {};   % 'M1' or 'M2'

for s = 1:numel(all_sessions)
    sess   = all_sessions{s};
    e_rows = session_data.(sess).early_rows;
    l_rows = session_data.(sess).late_rows;

    if isempty(e_rows) || isempty(l_rows), continue; end

    E = cell2mat(e_rows');
    L = cell2mat(l_rows');

    if size(E,2) ~= size(L,2)
        fprintf('Warning: session %s pair count mismatch (E=%d, L=%d), skipping.\n', ...
            sess, size(E,2), size(L,2));
        continue;
    end

    e_mean_per_pair = mean(E, 1)';
    l_mean_per_pair = mean(L, 1)';

    valid = ~isnan(e_mean_per_pair) & ~isnan(l_mean_per_pair);
    e_mean_per_pair = e_mean_per_pair(valid);
    l_mean_per_pair = l_mean_per_pair(valid);

    paired_early = [paired_early; e_mean_per_pair];
    paired_late  = [paired_late;  l_mean_per_pair];

    % Determine whether the session belongs to M1 or M2.
    if any(strcmp(sess, M1_sessions))
        mk = 'M1';
    else
        mk = 'M2';
    end

    for p = 1:numel(e_mean_per_pair)
        pair_session_label{end+1} = sess;
        pair_monkey_label{end+1}  = mk;
    end
end

n_pairs_total = numel(paired_early);
fprintf('Total paired (session x pair) units: %d\n', n_pairs_total);

%% -------------------- 6. Statistics (all + per-monkey) --------------------
diff_vals = paired_late - paired_early;

fprintf('\n======= All pairs: Late - Early =======\n');
fprintf('n = %d pairs\n', n_pairs_total);
fprintf('Early: mean=%.3f  SEM=%.3f\n', mean(paired_early), std(paired_early)/sqrt(n_pairs_total));
fprintf('Late:  mean=%.3f  SEM=%.3f\n', mean(paired_late),  std(paired_late)/sqrt(n_pairs_total));
fprintf('Mean Delta = %.3f +/- %.3f SEM\n', mean(diff_vals), std(diff_vals)/sqrt(n_pairs_total));

[~, p_paired_t, ~, stats_pt] = ttest(paired_late, paired_early);
p_sr = signrank(paired_late, paired_early);
fprintf('Paired t: p=%.6f, t(%d)=%.3f\n', p_paired_t, stats_pt.df, stats_pt.tstat);
fprintf('Signrank: p=%.6f\n', p_sr);

% --- Per-monkey stats ---
for mk = {'M1','M2'}
    idx = strcmp(pair_monkey_label, mk{1});
    pe  = paired_early(idx);
    pl  = paired_late(idx);
    dv  = pl - pe;
    fprintf('\n--- %s (n=%d pairs) ---\n', mk{1}, sum(idx));
    fprintf('Early: mean=%.3f  SEM=%.3f\n', mean(pe), std(pe)/sqrt(numel(pe)));
    fprintf('Late:  mean=%.3f  SEM=%.3f\n', mean(pl), std(pl)/sqrt(numel(pl)));
    fprintf('Mean Delta = %.3f\n', mean(dv));
    if numel(pe) >= 2
        [~, pm_t] = ttest(pl, pe);
        pm_sr = signrank(pl, pe);
        fprintf('Paired t p=%.6f  Signrank p=%.6f\n', pm_t, pm_sr);
    end
end

%% -------------------- 7. Figure --------------------
violin_col = [0.31 0.58 0.80;   % Blue identifies Early.
              0.93 0.40 0.36];  % Red identifies Late.

fig = figure('Color','w','Position',[200 150 560 580]);
hold on;

% Thin lines connecting each pair, colored by monkey
idx_M1 = strcmp(pair_monkey_label, 'M1');
idx_M2 = strcmp(pair_monkey_label, 'M2');

% M1 lines
for i = find(idx_M1)
    plot([1 2], [paired_early(i) paired_late(i)], '-', ...
        'Color', [col_M1 0.30], ...
        'LineWidth', 0.9);
end

% M2 lines
for i = find(idx_M2)
    plot([1 2], [paired_early(i) paired_late(i)], '-', ...
        'Color', [col_M2 0.30], ...
        'LineWidth', 0.9);
end

% Draw two violin plots above the connecting lines.
draw_single_violin(1, paired_early, 0.38, violin_col(1,:));
draw_single_violin(2, paired_late,  0.38, violin_col(2,:));

% Significance annotation
if     p_paired_t < 0.001, sig_str = '***';
elseif p_paired_t < 0.01,  sig_str = '**';
elseif p_paired_t < 0.05,  sig_str = '*';
else,                      sig_str = 'ns';
end

y_all = [paired_early; paired_late];
rng_y = range(y_all);
y_bar = max(y_all) + 0.04 * rng_y;
y_txt = y_bar + 0.018 * rng_y;

plot([1 2], [y_bar y_bar], 'k-', 'LineWidth', 1.3);
plot([1 1], [y_bar y_bar - 0.01*rng_y], 'k-', 'LineWidth', 1.3);
plot([2 2], [y_bar y_bar - 0.01*rng_y], 'k-', 'LineWidth', 1.3);
text(1.5, y_txt, sig_str, ...
    'HorizontalAlignment','center','FontSize',17,'FontWeight','bold');

% Axes
xlim([0.5 2.5]);
set(gca,'XTick',[1 2],'XTickLabel',{'Early','Late'},'FontSize',15);
ylabel('Pairwise decoding accuracy (%)','FontSize',15);
title(sprintf('Pairwise between-category Decoding: Early vs Late'), 'FontSize',20);

% Legend
h_e  = patch(NaN,NaN, violin_col(1,:), 'FaceAlpha',0.45, 'EdgeColor',violin_col(1,:)*0.75);
h_l  = patch(NaN,NaN, violin_col(2,:), 'FaceAlpha',0.45, 'EdgeColor',violin_col(2,:)*0.75);
h_m1 = plot(NaN,NaN, '-', 'Color', col_M1, 'LineWidth', 2.5);
h_m2 = plot(NaN,NaN, '-', 'Color', col_M2, 'LineWidth', 2.5);

legend([h_e, h_l, h_m1, h_m2], ...
    {'Early','Late', ...
     sprintf('M1'), ...
     sprintf('M2')}, ...
    'Location','southwest','FontSize',10,'Box','off');
box off;

%% -------------------- 8. Save --------------------
pair_table = table( ...
    string(pair_session_label'), string(pair_monkey_label'), ...
    paired_early, paired_late, diff_vals, ...
    'VariableNames',{'session','monkey','early_acc','late_acc','delta'});

save(fullfile(save_dir,'paired_per_pair_data.mat'), ...
    'paired_early','paired_late','diff_vals', ...
    'pair_session_label','pair_monkey_label');
writetable(pair_table, fullfile(save_dir,'paired_per_pair_data.csv'));

saveas(fig, fullfile(save_dir,'violin_early_vs_late_monkey.png'));
saveas(fig, fullfile(save_dir,'violin_early_vs_late_monkey.pdf'));

fprintf('\nDone. Saved to:\n%s\n', save_dir);

%% ==================== Local helper function ====================

function draw_single_violin(x_center, vals, half_width, col)
if numel(vals) < 4, return; end
[f, xi] = ksdensity(vals);
f = f / max(f) * half_width;
fill([x_center-f, fliplr(x_center+f)], [xi, fliplr(xi)], ...
     col, 'FaceAlpha',0.45, 'EdgeColor',col*0.75, 'LineWidth',1.2);
q25 = prctile(vals,25);  q75 = prctile(vals,75);
med = median(vals);       mn  = mean(vals);
line([x_center x_center],[q25 q75],'Color',col*0.55,'LineWidth',4);
line([x_center-half_width*0.35, x_center+half_width*0.35],[med med], ...
     'Color','w','LineWidth',2.5);
plot(x_center, mn,'o','MarkerFaceColor','w','MarkerEdgeColor',col*0.6, ...
     'MarkerSize',7,'LineWidth',1.8);
end