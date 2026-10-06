%% =====================================================================
%  SDI Early vs Late
%
%  Combined PNAS-style plot
%
%  Delta SDI = Late - Early
%
%  ALL category pairs are pooled in one plot
%
%  Blue   = M1 (s1, s2)
%  Orange = M2 (s3, s4, s5)
%
%  Small dots = individual session x category-pair observations
%  Large dot  = overall mean Delta SDI
%  Error bar  = bootstrap 95% CI
%  Dashed line = Delta SDI = 0
%
% =====================================================================

clear;
clc;
close all;


%% =====================================================================
%  1. Path setup
% =====================================================================

result_root = ...
    figure_input('var_sdi');

save_dir = figure_output('fig3F_delta');

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end


%% =====================================================================
%  2. Define training stages
% =====================================================================

stage_groups = struct();

stage_groups.s1.early  = {'041116','041216'};
stage_groups.s1.middle = {'041816'};
stage_groups.s1.late   = {'042916'};

stage_groups.s2.early  = {'031617'};
stage_groups.s2.middle = {''};
stage_groups.s2.late   = {'031917'};

stage_groups.s3.early  = {'071217'};
stage_groups.s3.middle = {'071617'};
stage_groups.s3.late   = {'071917','072217'};

stage_groups.s4.early  = {'092817'};
stage_groups.s4.middle = {'093017'};
stage_groups.s4.late   = {'100217'};

stage_groups.s5.early  = {'081318'};
stage_groups.s5.middle = {'081718','082018'};
stage_groups.s5.late   = {'082418'};


session_order = ...
    {'s1','s2','s3','s4','s5'};


%% =====================================================================
%  3. Monkey assignment
% =====================================================================

M1_sessions = {'s1','s2'};
M2_sessions = {'s3','s4','s5'};


% Publication-friendly colors
col_M1 = [0.18 0.45 0.70];
col_M2 = [0.85 0.40 0.12];


%% =====================================================================
%  4. Find SDI files
% =====================================================================

files_all = dir( ...
    fullfile( ...
        result_root, ...
        '**', ...
        'SDI_Results_*.mat'));


if isempty(files_all)

    error( ...
        'No SDI_Results_*.mat files found under:\n%s', ...
        result_root);

end


% Keep variation files only
is_var = arrayfun( ...
    @(x) contains(lower(x.folder),'var'), ...
    files_all);


files = ...
    files_all(is_var);


if isempty(files)

    error('No variation SDI files found.');

end


fprintf( ...
    '\nFound %d variation SDI files.\n', ...
    numel(files));


%% =====================================================================
%  5. Initialize session storage
% =====================================================================

session_data = struct();


for s = 1:numel(session_order)

    sess = session_order{s};

    session_data.(sess).early_rows = {};
    session_data.(sess).late_rows  = {};

end


%% =====================================================================
%  6. Load SDI data
% =====================================================================

for f = 1:numel(files)

    fileName = ...
        files(f).name;


    filePath = ...
        fullfile( ...
            files(f).folder, ...
            fileName);


    %% -----------------------------------------------------------------
    % Parse session
    % ------------------------------------------------------------------

    tok_sess = regexp( ...
        filePath, ...
        '[\\/](s\d+)[\\/]', ...
        'tokens', ...
        'once');


    if isempty(tok_sess)

        fprintf( ...
            'Skipping: cannot parse session | %s\n', ...
            fileName);

        continue;

    end


    sess = ...
        tok_sess{1};


    if ~ismember(sess,session_order)
        continue;
    end


    %% -----------------------------------------------------------------
    % Parse date
    % ------------------------------------------------------------------

    tok_date = regexp( ...
        fileName, ...
        '(\d{6})', ...
        'tokens', ...
        'once');


    if isempty(tok_date)

        fprintf( ...
            'Skipping: cannot parse date | %s\n', ...
            fileName);

        continue;

    end


    date_code = ...
        tok_date{1};


    %% -----------------------------------------------------------------
    % Determine stage
    % ------------------------------------------------------------------

    stage = '';


    if any(strcmp( ...
            date_code, ...
            stage_groups.(sess).early))

        stage = 'early';


    elseif any(strcmp( ...
            date_code, ...
            stage_groups.(sess).middle))

        stage = 'middle';


    elseif any(strcmp( ...
            date_code, ...
            stage_groups.(sess).late))

        stage = 'late';


    else

        fprintf( ...
            'Skipping: unassigned date | %s | %s\n', ...
            sess, ...
            date_code);

        continue;

    end


    % Middle sessions are not used
    if strcmp(stage,'middle')
        continue;
    end


    %% -----------------------------------------------------------------
    % Load
    % ------------------------------------------------------------------

    try

        S = load( ...
            filePath, ...
            'SDI_results');

    catch ME

        fprintf( ...
            'Skipping: load failed | %s | %s\n', ...
            fileName, ...
            ME.message);

        continue;

    end


    if ~isfield(S,'SDI_results')

        fprintf( ...
            'Skipping: SDI_results missing | %s\n', ...
            fileName);

        continue;

    end


    %% -----------------------------------------------------------------
    % Extract matrix
    % ------------------------------------------------------------------

    if isfield(S.SDI_results,'Global') && ...
       isfield(S.SDI_results.Global,'SDI_matrix')

        SDI_mat = ...
            S.SDI_results.Global.SDI_matrix;

    else

        fprintf( ...
            'Skipping: SDI_matrix missing | %s\n', ...
            fileName);

        continue;

    end


    %% -----------------------------------------------------------------
    % Validate
    % ------------------------------------------------------------------

    if isempty(SDI_mat) || ...
       ~ismatrix(SDI_mat) || ...
       size(SDI_mat,1) ~= size(SDI_mat,2)

        fprintf( ...
            'Skipping: invalid SDI_matrix | %s\n', ...
            fileName);

        continue;

    end


    %% -----------------------------------------------------------------
    % Extract upper triangle
    %
    % IMPORTANT:
    % Do NOT remove NaNs here.
    % Keeping them preserves category-pair identity across dates.
    % ------------------------------------------------------------------

    mask_ut = ...
        triu(true(size(SDI_mat)),1);


    pair_vals = ...
        SDI_mat(mask_ut);


    pair_vals = ...
        pair_vals(:)';


    %% -----------------------------------------------------------------
    % Store
    % ------------------------------------------------------------------

    if strcmp(stage,'early')

        session_data.(sess).early_rows{end+1,1} = ...
            pair_vals;


    elseif strcmp(stage,'late')

        session_data.(sess).late_rows{end+1,1} = ...
            pair_vals;

    end

end


%% =====================================================================
%  7. Build paired Delta SDI
% =====================================================================

paired_early = [];
paired_late  = [];

delta_all = [];

session_labels = {};
monkey_labels  = {};

pair_numbers = [];


fprintf('\n');
fprintf('====================================================\n');
fprintf(' Session summary\n');
fprintf('====================================================\n');


for s = 1:numel(session_order)

    sess = ...
        session_order{s};


    e_rows = ...
        session_data.(sess).early_rows;


    l_rows = ...
        session_data.(sess).late_rows;


    %% -----------------------------------------------------------------
    % Need Early and Late
    % ------------------------------------------------------------------

    if isempty(e_rows) || isempty(l_rows)

        warning( ...
            '%s does not contain both Early and Late.', ...
            sess);

        continue;

    end


    %% -----------------------------------------------------------------
    % Validate pair counts
    % ------------------------------------------------------------------

    e_lengths = ...
        cellfun(@numel,e_rows);


    l_lengths = ...
        cellfun(@numel,l_rows);


    if numel(unique(e_lengths)) ~= 1 || ...
       numel(unique(l_lengths)) ~= 1

        warning( ...
            '%s: inconsistent number of pairs across dates.', ...
            sess);

        continue;

    end


    if e_lengths(1) ~= l_lengths(1)

        warning( ...
            '%s: Early/Late pair count mismatch.', ...
            sess);

        continue;

    end


    %% -----------------------------------------------------------------
    % Rows = dates
    % Columns = category pairs
    % ------------------------------------------------------------------

    E = ...
        vertcat(e_rows{:});


    L = ...
        vertcat(l_rows{:});


    %% -----------------------------------------------------------------
    % Average dates within each stage
    % ------------------------------------------------------------------

    early_mean_per_pair = ...
        mean(E,1,'omitnan')';


    late_mean_per_pair = ...
        mean(L,1,'omitnan')';


    %% -----------------------------------------------------------------
    % Keep pair identity
    % ------------------------------------------------------------------

    valid = ...
        ~isnan(early_mean_per_pair) & ...
        ~isnan(late_mean_per_pair);


    original_pair_number = ...
        find(valid);


    early_vals = ...
        early_mean_per_pair(valid);


    late_vals = ...
        late_mean_per_pair(valid);


    delta_vals = ...
        late_vals - early_vals;


    %% -----------------------------------------------------------------
    % Monkey
    % ------------------------------------------------------------------

    if ismember(sess,M1_sessions)

        monkey = 'M1';

    elseif ismember(sess,M2_sessions)

        monkey = 'M2';

    else

        error( ...
            'Unknown monkey for session %s.', ...
            sess);

    end


    %% -----------------------------------------------------------------
    % Append
    % ------------------------------------------------------------------

    paired_early = [ ...
        paired_early;
        early_vals];


    paired_late = [ ...
        paired_late;
        late_vals];


    delta_all = [ ...
        delta_all;
        delta_vals];


    pair_numbers = [ ...
        pair_numbers;
        original_pair_number(:)];


    session_labels = [ ...
        session_labels;
        repmat({sess},numel(delta_vals),1)];


    monkey_labels = [ ...
        monkey_labels;
        repmat({monkey},numel(delta_vals),1)];


    %% -----------------------------------------------------------------
    % Session console output
    % ------------------------------------------------------------------

    fprintf( ...
        '%s (%s): n = %d | mean Delta = %.4f\n', ...
        sess, ...
        monkey, ...
        numel(delta_vals), ...
        mean(delta_vals,'omitnan'));

end


%% =====================================================================
%  8. Basic check
% =====================================================================

n_total = ...
    numel(delta_all);


if n_total == 0

    error('No valid paired Delta SDI values.');

end


fprintf('\n');
fprintf( ...
    'Total session x pair observations = %d\n', ...
    n_total);


%% =====================================================================
%  9. Overall statistics
% =====================================================================

overall_mean = ...
    mean(delta_all,'omitnan');


overall_sem = ...
    std(delta_all,'omitnan') / sqrt(n_total);


%% ---------------------------------------------------------------------
% Paired t test: mathematically equivalent to one-sample Delta vs zero
% ----------------------------------------------------------------------

[~,p_t,ci_t,stats_t] = ...
    ttest( ...
        paired_late, ...
        paired_early);


%% ---------------------------------------------------------------------
% Signed-rank test
% ----------------------------------------------------------------------

p_sr = ...
    signrank(delta_all,0);


fprintf('\n');
fprintf('====================================================\n');
fprintf(' Overall Early vs Late\n');
fprintf('====================================================\n');

fprintf( ...
    'n = %d\n', ...
    n_total);

fprintf( ...
    'Mean Delta SDI = %.4f +/- %.4f SEM\n', ...
    overall_mean, ...
    overall_sem);

fprintf( ...
    'Paired t: t(%d) = %.3f, p = %.6g\n', ...
    stats_t.df, ...
    stats_t.tstat, ...
    p_t);

fprintf( ...
    'Signed-rank: p = %.6g\n', ...
    p_sr);


%% =====================================================================
%  10. Per-monkey descriptive statistics
% =====================================================================

idx_M1 = ...
    strcmp(monkey_labels,'M1');


idx_M2 = ...
    strcmp(monkey_labels,'M2');


fprintf('\n');

fprintf( ...
    'M1: n=%d, mean Delta=%.4f\n', ...
    sum(idx_M1), ...
    mean(delta_all(idx_M1),'omitnan'));


fprintf( ...
    'M2: n=%d, mean Delta=%.4f\n', ...
    sum(idx_M2), ...
    mean(delta_all(idx_M2),'omitnan'));


%% =====================================================================
%  11. Bootstrap overall 95% CI
% =====================================================================

rng(2026);

n_boot = 10000;


boot_mean = ...
    nan(n_boot,1);


for b = 1:n_boot

    idx_boot = ...
        randi(n_total,n_total,1);


    boot_mean(b) = ...
        mean(delta_all(idx_boot),'omitnan');

end


overall_CI = ...
    prctile( ...
        boot_mean, ...
        [2.5 97.5]);


fprintf( ...
    'Bootstrap 95%% CI = [%.4f, %.4f]\n', ...
    overall_CI(1), ...
    overall_CI(2));


%% =====================================================================
%  12. Figure setup
% =====================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 7.2 7.5], ...
    'PaperPositionMode','auto');


ax = axes(fig);

hold(ax,'on');


%% =====================================================================
%  13. Determine Y limits
% =====================================================================

data_min = ...
    min(delta_all);


data_max = ...
    max(delta_all);


data_range = ...
    data_max - data_min;


if data_range <= 0

    data_range = ...
        max(abs(data_max),0.1);

end


y_lower = ...
    min(data_min,0) - 0.12*data_range;


y_upper = ...
    max(data_max,0) + 0.24*data_range;


%% =====================================================================
%  14. Zero line
% =====================================================================

yline( ...
    ax, ...
    0, ...
    '--', ...
    'Color',[0.45 0.45 0.45], ...
    'LineWidth',0.75);


%% =====================================================================
%  15. Individual data points
% =====================================================================

% M1 slightly left
% M2 slightly right

x_center = 1;

offset_M1 = -0.07;
offset_M2 =  0.07;

jitter_width = 0.15;


% Reproducible point placement
rng(7);


%% -------------------- M1 --------------------

vals_M1 = ...
    delta_all(idx_M1);


n_M1 = ...
    numel(vals_M1);


jitter_M1 = ...
    (rand(n_M1,1)-0.5) * 2*jitter_width;


x_M1 = ...
    x_center + offset_M1 + jitter_M1;


scatter( ...
    ax, ...
    x_M1, ...
    vals_M1, ...
    17, ...
    'o', ...
    'MarkerFaceColor',col_M1, ...
    'MarkerEdgeColor','none', ...
    'MarkerFaceAlpha',0.58);


%% -------------------- M2 --------------------

vals_M2 = ...
    delta_all(idx_M2);


n_M2 = ...
    numel(vals_M2);


jitter_M2 = ...
    (rand(n_M2,1)-0.5) * 2*jitter_width;


x_M2 = ...
    x_center + offset_M2 + jitter_M2;


scatter( ...
    ax, ...
    x_M2, ...
    vals_M2, ...
    17, ...
    'o', ...
    'MarkerFaceColor',col_M2, ...
    'MarkerEdgeColor','none', ...
    'MarkerFaceAlpha',0.58);


%% =====================================================================
%  16. Overall 95% CI
% =====================================================================

summary_x = ...
    1;


plot( ...
    ax, ...
    [summary_x summary_x], ...
    [overall_CI(1) overall_CI(2)], ...
    '-', ...
    'Color',[0.10 0.10 0.10], ...
    'LineWidth',1.6);


cap = 0.055;


plot( ...
    ax, ...
    [summary_x-cap summary_x+cap], ...
    [overall_CI(1) overall_CI(1)], ...
    '-', ...
    'Color',[0.10 0.10 0.10], ...
    'LineWidth',1.2);


plot( ...
    ax, ...
    [summary_x-cap summary_x+cap], ...
    [overall_CI(2) overall_CI(2)], ...
    '-', ...
    'Color',[0.10 0.10 0.10], ...
    'LineWidth',1.2);


%% =====================================================================
%  17. Overall mean
% =====================================================================

scatter( ...
    ax, ...
    summary_x, ...
    overall_mean, ...
    58, ...
    'o', ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',[0.05 0.05 0.05], ...
    'LineWidth',1.4);


% Central small dot
scatter( ...
    ax, ...
    summary_x, ...
    overall_mean, ...
    12, ...
    'o', ...
    'MarkerFaceColor',[0.05 0.05 0.05], ...
    'MarkerEdgeColor','none');


%% =====================================================================
%  18. Statistical annotation
%
%  Choose which statistic is displayed here.
%
%  Recommended for figure:
%  p_display = p_sr;
%
%  If paired t-test is your predefined primary analysis:
%  p_display = p_t;
% =====================================================================

p_display = ...
    p_sr;


if p_display < 0.001

    p_text = ...
        'p < 0.001';

elseif p_display < 0.01

    p_text = ...
        sprintf('p = %.3f',p_display);

else

    p_text = ...
        sprintf('p = %.2f',p_display);

end


text( ...
    ax, ...
    1, ...
    data_max + 0.14*data_range, ...
    p_text, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontName','Arial', ...
    'FontSize',8);


%% =====================================================================
%  19. Axis formatting
% =====================================================================

xlim( ...
    ax, ...
    [0.55 1.45]);


ylim( ...
    ax, ...
    [y_lower y_upper]);


set( ...
    ax, ...
    'XTick',1, ...
    'XTickLabel',{'All pairs'}, ...
    'FontName','Arial', ...
    'FontSize',8, ...
    'LineWidth',0.65, ...
    'TickDir','out', ...
    'TickLength',[0.025 0.025], ...
    'Box','off', ...
    'Layer','top');


ylabel( ...
    ax, ...
    '\Delta Pairwise SDI (Late - Early)', ...
    'FontName','Arial', ...
    'FontSize',8);


xlabel(ax,'');

title(ax,'');


%% =====================================================================
%  20. Legend
% =====================================================================

h_M1 = scatter( ...
    ax, ...
    NaN, ...
    NaN, ...
    18, ...
    'o', ...
    'MarkerFaceColor',col_M1, ...
    'MarkerEdgeColor','none');


h_M2 = scatter( ...
    ax, ...
    NaN, ...
    NaN, ...
    18, ...
    'o', ...
    'MarkerFaceColor',col_M2, ...
    'MarkerEdgeColor','none');


lgd = legend( ...
    ax, ...
    [h_M1 h_M2], ...
    {'M1','M2'}, ...
    'Location','northeast', ...
    'FontName','Arial', ...
    'FontSize',7, ...
    'Box','off');


try
    lgd.ItemTokenSize = [9 6];
catch
end


%% =====================================================================
%  21. Optional panel letter
% =====================================================================

panel_letter = 'A';


text( ...
    ax, ...
    -0.15, ...
    1.04, ...
    panel_letter, ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontName','Arial', ...
    'FontSize',10, ...
    'FontWeight','bold', ...
    'Clipping','off');


%% =====================================================================
%  22. Layout
% =====================================================================

ax.Position = [ ...
    0.25 ...
    0.17 ...
    0.68 ...
    0.75];


%% =====================================================================
%  23. Save table
% =====================================================================

result_table = table( ...
    string(session_labels), ...
    string(monkey_labels), ...
    pair_numbers, ...
    paired_early, ...
    paired_late, ...
    delta_all, ...
    'VariableNames',{ ...
        'session', ...
        'monkey', ...
        'pair_number', ...
        'early_SDI', ...
        'late_SDI', ...
        'delta_SDI'});


writetable( ...
    result_table, ...
    fullfile( ...
        save_dir, ...
        'SDI_combined_delta_data.csv'));


%% =====================================================================
%  24. Save MAT
% =====================================================================

save( ...
    fullfile( ...
        save_dir, ...
        'SDI_combined_delta_data.mat'), ...
    'paired_early', ...
    'paired_late', ...
    'delta_all', ...
    'session_labels', ...
    'monkey_labels', ...
    'pair_numbers', ...
    'result_table', ...
    'overall_mean', ...
    'overall_CI', ...
    'p_t', ...
    'p_sr');


%% =====================================================================
%  25. Save figure
% =====================================================================

exportgraphics( ...
    fig, ...
    fullfile( ...
        save_dir, ...
        'Fig_Delta_SDI_combined_M1_M2.pdf'), ...
    'ContentType','vector');


exportgraphics( ...
    fig, ...
    fullfile( ...
        save_dir, ...
        'Fig_Delta_SDI_combined_M1_M2.png'), ...
    'Resolution',600);


savefig( ...
    fig, ...
    fullfile( ...
        save_dir, ...
        'Fig_Delta_SDI_combined_M1_M2.fig'));


fprintf('\n');
fprintf('====================================================\n');
fprintf('Done.\n');
fprintf('Saved to:\n%s\n',save_dir);
fprintf('====================================================\n');