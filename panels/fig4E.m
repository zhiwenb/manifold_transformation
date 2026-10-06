%% =====================================================================
%  Batch analysis — Early vs Late, per-pair lines colored by monkey
%  M1 = s1, s2  |  M2 = s3, s4, s5
% =====================================================================

clear; clc; close all;

%% -------------------- 1. Path setup --------------------
result_root = figure_input('original_global_decoding');

save_dir = figure_output('fig4E');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% -------------------- 2. Define stage groups --------------------
stage_groups = struct();

stage_groups.s1.early  = {'040716','041116'};
stage_groups.s1.middle = {'041416','041716','041816','041916'};
stage_groups.s1.late   = {'042016','042116','042416','042516'};

stage_groups.s2.early  = {'031217'};
stage_groups.s2.middle = {'031517'};
stage_groups.s2.late   = {'032317'};

stage_groups.s3.early  = {'071117'};
stage_groups.s3.middle = {'071517'};
stage_groups.s3.late   = {'071717','072017'};

stage_groups.s4.early  = {'092417'};
stage_groups.s4.middle = {'092717','092917'};
stage_groups.s4.late   = {'100117'};

stage_groups.s5.early  = {'081218'};
stage_groups.s5.middle = {'081518','081818'};
stage_groups.s5.late   = {'082218'};

all_sessions = fieldnames(stage_groups);

% -------- Monkey group assignment --------
M1_sessions = {'s1','s2'};
M2_sessions = {'s3','s4','s5'};

col_M1 = [0.20 0.55 0.82];   % Blue tones identify M1.
col_M2 = [0.85 0.33 0.10];   % Orange tones identify M2.

%% -------------------- 3. Find files --------------------
files_all = dir(fullfile(result_root, '**', 'Decoder_Acc_*.mat'));
if isempty(files_all)
    error('No Decoder_Acc_*.mat files found under: %s', result_root);
end
is_variation = arrayfun(@(x) contains(x.name, '_prototype_'), files_all);
files = files_all(is_variation);
if isempty(files)
    error('No prototype decoding files found.');
end
fprintf('Found %d prototype result files.\n', numel(files));

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
pair_monkey_label  = {};

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

%% -------------------- 6. Statistics --------------------
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

%% -------------------- 6b. Print every pair --------------------
fprintf('\n======= Per-pair detail =======\n');
fprintf('%-8s %-6s %-10s %10s %10s %10s\n', ...
    'session','monkey','pair_idx','early','late','delta');
fprintf('%s\n', repmat('-', 1, 58));

for i = 1:n_pairs_total
    fprintf('%-8s %-6s %-10d %10.3f %10.3f %10.3f\n', ...
        pair_session_label{i}, ...
        pair_monkey_label{i}, ...
        i, ...
        paired_early(i), ...
        paired_late(i), ...
        diff_vals(i));
end

fprintf('%s\n', repmat('-', 1, 58));
fprintf('%-8s %-6s %-10s %10.3f %10.3f %10.3f\n', ...
    'MEAN','','', mean(paired_early), mean(paired_late), mean(diff_vals));
fprintf('%-8s %-6s %-10s %10.3f %10.3f %10.3f\n', ...
    'SEM','','', ...
    std(paired_early)/sqrt(n_pairs_total), ...
    std(paired_late)/sqrt(n_pairs_total), ...
    std(diff_vals)/sqrt(n_pairs_total));

%% -------------------- 7. Figure --------------------
violin_col = [0.31 0.58 0.80;   % Blue identifies Early.
              0.93 0.40 0.36];  % Red identifies Late.

fig = figure('Color','w','Position',[200 150 560 580]);
hold on;

idx_M1 = strcmp(pair_monkey_label, 'M1');
idx_M2 = strcmp(pair_monkey_label, 'M2');

% Thin connecting lines for M1
for i = find(idx_M1)
    plot([1 2], [paired_early(i) paired_late(i)], '-', ...
        'Color', [col_M1 0.30], 'LineWidth', 0.9);
end

% Thin connecting lines for M2
for i = find(idx_M2)
    plot([1 2], [paired_early(i) paired_late(i)], '-', ...
        'Color', [col_M2 0.30], 'LineWidth', 0.9);
end

% Two violin plots
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
title('Pairwise between-category decoding accuracy','FontSize',20);

% Legend
h_e  = patch(NaN,NaN, violin_col(1,:), 'FaceAlpha',0.45, 'EdgeColor',violin_col(1,:)*0.75);
h_l  = patch(NaN,NaN, violin_col(2,:), 'FaceAlpha',0.45, 'EdgeColor',violin_col(2,:)*0.75);
h_m1 = plot(NaN,NaN, '-', 'Color', col_M1, 'LineWidth', 2.5);
h_m2 = plot(NaN,NaN, '-', 'Color', col_M2, 'LineWidth', 2.5);

legend([h_e, h_l, h_m1, h_m2], ...
    {'Early','Late','M1','M2'}, ...
    'Location','southwest','FontSize',10,'Box','off');
box off;

ylim([40 120]); % Match the reference panel axis range.

%% -------------------- 8. Save --------------------
pair_index = (1:n_pairs_total)';

pair_table = table( ...
    pair_index, ...
    string(pair_session_label'), ...
    string(pair_monkey_label'), ...
    paired_early, ...
    paired_late, ...
    diff_vals, ...
    'VariableNames',{'pair_idx','session','monkey','early_acc','late_acc','delta'});

% Sort by session, then pair_idx.
pair_table = sortrows(pair_table, {'session','pair_idx'});

fprintf('\n======= Sorted pair table =======\n');
disp(pair_table);

save(fullfile(save_dir,'paired_per_pair_data.mat'), ...
    'paired_early','paired_late','diff_vals', ...
    'pair_session_label','pair_monkey_label','pair_table');
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