%% Batch Net Relative RSA Improvement - BIG CLEAN VERSION
% Logic:
%   - Treat every training day (Day 2+) as an independent sample
%   - Polar: use s1, s2, s3, s4, s5
%   - Hyperbolic: use ONLY s1 and s5

clear; clc; close all;

%% 1. Configuration
rootDir = figure_input('global_fr');
saveDir = figure_output('fig5C');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

idx_G = 1; idx_H = 2; idx_R = 3; idx_S = 4; idx_T = 5;

analysis_groups = {
    struct('name', 'Polar (S+R+T)', 'indices', [idx_R, idx_S, idx_T], 'p_x', 'ang_freq',  'p_y', 'radial_freq'), ...
    struct('name', 'Grating (G)',   'indices', [idx_G],               'p_x', 'ori_deg',   'p_y', 'freq_cpd'), ...
    struct('name', 'Hyperbolic (H)','indices', [idx_H],               'p_x', 'ori_deg',   'p_y', 'radial_freq')
};

control_group_idx = 2;

polar_M1_ids = [1, 2];
polar_M2_ids = [3, 4, 5];
hyp_M1_ids   = [1];
hyp_M2_ids   = [5];

pooled_data = struct();
pooled_data.M1.Polar = [];
pooled_data.M1.Hyp   = [];
pooled_data.M2.Polar = [];
pooled_data.M2.Hyp   = [];

%% 2. Batch Processing
sessionDirs = dir(fullfile(rootDir, 's*'));
sessionDirs = sessionDirs([sessionDirs.isdir]);
fprintf('Found %d sessions.\n', length(sessionDirs));

for s = 1:length(sessionDirs)
    sessionName = sessionDirs(s).name;
    sessionPath = fullfile(rootDir, sessionName, 'global');

    sid_tokens = regexp(sessionName, 's(\d+)', 'tokens');
    if isempty(sid_tokens)
        continue;
    end
    sessionID = str2double(sid_tokens{1}{1});

    matFiles = dir(fullfile(sessionPath, '*.mat'));
    % Match original Fig5 s1 scope: baseline 041116, excluding extra 040716/041216.
    if sessionID==1,matFiles=matFiles(~contains({matFiles.name},'040716') & ~contains({matFiles.name},'041216'));end
    if isempty(matFiles)
        continue;
    end

    dates = cell(length(matFiles), 1);
    for f = 1:length(matFiles)
        d = regexp(matFiles(f).name, '\d{6}', 'match');
        if isempty(d)
            dates{f} = '000000';
        else
            dates{f} = d{1};
        end
    end

    [~, sortIdx] = sort(dates);
    matFiles = matFiles(sortIdx);
    num_days = length(matFiles);

    if num_days < 2
        continue;
    end

    rsa_raw = nan(num_days, 3);

    for f = 1:num_days
        try
            load(fullfile(sessionPath, matFiles(f).name), 'FR_by_category');
            [P, ~, ~] = getPara(FR_by_category);

            for g = 1:3
                grp = analysis_groups{g};
                [X, px, py] = stitch_data(FR_by_category, P, grp.indices, grp.p_x, grp.p_y);

                if ~isempty(X) && size(X,1) > 1
                    rsa_raw(f,g) = calculate_rsa_mixed_model(X, px, py, grp.p_x, grp.p_y);
                end
            end
        catch ME
            warning('Skipping %s / %s: %s', sessionName, matFiles(f).name, ME.message);
        end
    end

    baseline = rsa_raw(1,:);
    if any(isnan(baseline)) || any(baseline == 0)
        fprintf('  Session s%d skipped (invalid baseline).\n', sessionID);
        continue;
    end

    rel_rate   = (rsa_raw - baseline) ./ (rsa_raw + baseline);
    net_metric = rel_rate - rel_rate(:, control_group_idx);

    days_polar = net_metric(2:end, 1);
    days_polar = days_polar(~isnan(days_polar));

    days_hyp = net_metric(2:end, 3);
    days_hyp = days_hyp(~isnan(days_hyp));

    if ismember(sessionID, polar_M1_ids)
        pooled_data.M1.Polar = [pooled_data.M1.Polar; days_polar];
    elseif ismember(sessionID, polar_M2_ids)
        pooled_data.M2.Polar = [pooled_data.M2.Polar; days_polar];
    end

    if ismember(sessionID, hyp_M1_ids)
        pooled_data.M1.Hyp = [pooled_data.M1.Hyp; days_hyp];
    elseif ismember(sessionID, hyp_M2_ids)
        pooled_data.M2.Hyp = [pooled_data.M2.Hyp; days_hyp];
    end

    fprintf('  s%d: Polar +%d days, Hyp +%d days\n', ...
        sessionID, numel(days_polar), numel(days_hyp));
end

%% 3. Statistical Analysis
fprintf('\n=== Statistical Results (H0: Mean <= 0, one-tailed t-test) ===\n');

data_cells = {pooled_data.M1.Polar, pooled_data.M1.Hyp, ...
              pooled_data.M2.Polar, pooled_data.M2.Hyp};

group_labels     = {'M1','M1','M2','M2'};
condition_labels = {'Polar','Hyp','Polar','Hyp'};

stats_res = struct('mean',[],'sem',[],'p',[],'h',[],'n',[],'t',[],'df',[]);

for i = 1:4
    d = data_cells{i};

    if isempty(d)
        stats_res.mean(i) = 0;
        stats_res.sem(i)  = 0;
        stats_res.p(i)    = 1;
        stats_res.h(i)    = 0;
        stats_res.n(i)    = 0;
        stats_res.t(i)    = 0;
        stats_res.df(i)   = 0;
        continue;
    end

    stats_res.n(i)    = numel(d);
    stats_res.mean(i) = mean(d);
    stats_res.sem(i)  = std(d) / sqrt(numel(d));
    stats_res.df(i)   = numel(d) - 1;

    [h, p, ~, stat] = ttest(d, 0, 'Tail', 'right');

    stats_res.p(i) = p;
    stats_res.h(i) = h;
    stats_res.t(i) = stat.tstat;

    sig = '';
    if p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    end

    fprintf('%-10s: n=%d, Mean=%.3f, SEM=%.3f, t(%d)=%.2f, P=%.4f %s\n', ...
        [group_labels{i} ' ' condition_labels{i}], ...
        stats_res.n(i), stats_res.mean(i), stats_res.sem(i), ...
        stats_res.df(i), stats_res.t(i), p, sig);
end

%% 4. Figure — BIG CLEAN VERSION
c_polar     = [0.22 0.45 0.72];
c_hyp       = [0.84 0.44 0.16];
c_face      = [c_polar; c_hyp; c_polar; c_hyp];
c_edge      = c_face * 0.58;
c_dot       = c_face * 0.50;

fig_w = 170 / 25.4;
fig_h = 130 / 25.4;

hFig = figure('Color','w', ...
    'Units','inches', ...
    'Position',[1 1 fig_w fig_h], ...
    'PaperUnits','inches', ...
    'PaperSize',[fig_w fig_h], ...
    'PaperPosition',[0 0 fig_w fig_h]);

set(hFig, 'DefaultAxesFontName','Arial', ...
          'DefaultTextFontName','Arial');

ax = axes('Parent', hFig, ...
    'Units','normalized', ...
    'Position',[0.12 0.16 0.84 0.74]);
hold(ax,'on');

xpos = [1.0, 2.4, 4.8, 6.2];
bw   = 0.90;

for i = 1:4
    m = stats_res.mean(i);
    patch(ax, ...
        xpos(i)+[-bw/2, bw/2, bw/2, -bw/2], ...
        [0 0 m m], ...
        c_face(i,:), ...
        'EdgeColor', c_edge(i,:), ...
        'LineWidth', 1.5, ...
        'FaceAlpha', 0.84);
end

cw = 0.12;
for i = 1:4
    m  = stats_res.mean(i);
    se = stats_res.sem(i);

    line(ax, [xpos(i) xpos(i)], [m-se m+se], ...
        'Color','k', 'LineWidth',1.6);

    line(ax, xpos(i)+[-cw cw], [m+se m+se], ...
        'Color','k', 'LineWidth',1.6);

    line(ax, xpos(i)+[-cw cw], [m-se m-se], ...
        'Color','k', 'LineWidth',1.6);
end

rng(42);
spread_range = 0.18;
for i = 1:4
    d = data_cells{i};
    if isempty(d)
        continue;
    end

    n = numel(d);
    d_sorted = sort(d);

    if n == 1
        xj = xpos(i);
    else
        xj = xpos(i) + linspace(-spread_range, spread_range, n);
    end

    scatter(ax, xj(:), d_sorted(:), 42, ...
        'MarkerFaceColor', c_dot(i,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.65);
end

yline(ax, 0, 'Color',[0.30 0.30 0.30], 'LineWidth',1.0);

y_lo = -0.72;
y_hi = 0.95;
ylim(ax, [y_lo y_hi]);
xlim(ax, [0.2 7.0]);
yticks(ax, -0.6:0.2:0.8);

ax.YGrid         = 'on';
ax.GridAlpha     = 0.12;
ax.GridColor     = [0.45 0.45 0.45];
ax.GridLineStyle = ':';

y_star = 0.82;
for i = 1:4
    p = stats_res.p(i);

    if p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end

    star_col = [0.75 0.12 0.12];
    if p >= 0.05
        star_col = [0.55 0.55 0.55];
    end

    text(ax, xpos(i), y_star, sig, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',16, ...
        'FontWeight','bold', ...
        'Color',star_col);
end

ax.XTick = xpos;
ax.XTickLabel = {'','','',''};

y_brk = y_lo - 0.07;
y_lbl = y_lo - 0.17;

brk_pairs = [xpos(1)-bw/2, xpos(2)+bw/2; ...
             xpos(3)-bw/2, xpos(4)+bw/2];

for g = 1:2
    line(ax, brk_pairs(g,:), [y_brk y_brk], ...
        'Color','k', 'LineWidth',1.3, 'Clipping','off');
end

text(ax, mean(xpos(1:2)), y_lbl, 'M1', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'Color','k', ...
    'Clipping','off');

text(ax, mean(xpos(3:4)), y_lbl, 'M2', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'Color','k', ...
    'Clipping','off');

ylabel(ax, 'Net relative RSA improvement', ...
    'FontSize',16, ...
    'FontWeight','bold');

ax.Box        = 'off';
ax.TickDir    = 'out';
ax.TickLength = [0.018 0.018];
ax.LineWidth  = 1.2;
ax.XColor     = 'k';
ax.YColor     = 'k';
ax.FontSize   = 14;
ax.FontName   = 'Arial';

h1 = patch(ax, nan, nan, c_polar, ...
    'EdgeColor',c_polar*0.58, ...
    'LineWidth',1.3, ...
    'FaceAlpha',0.84);

h2 = patch(ax, nan, nan, c_hyp, ...
    'EdgeColor',c_hyp*0.58, ...
    'LineWidth',1.3, ...
    'FaceAlpha',0.84);

leg = legend(ax, [h1 h2], {'Polar','Hyperbolic'}, ...
    'Location','northoutside', ...
    'Orientation','horizontal', ...
    'FontSize',14, ...
    'Box','off');
leg.ItemTokenSize = [18 12];


%% 5. Save
fmts = {'-dpng','-r300'; '-dpdf','-vector'; '-depsc2','-vector'; '-dsvg',''};
exts = {'png','pdf','eps','svg'};

for k = 1:4
    fname = fullfile(saveDir, sprintf('RSA_Improvement_BigClean.%s', exts{k}));
    if isempty(fmts{k,2})
        print(hFig, fname, fmts{k,1});
    else
        print(hFig, fname, fmts{k,1}, fmts{k,2});
    end
end

savefig(hFig, fullfile(saveDir, 'RSA_Improvement_BigClean.fig'));
fprintf('\nAll figures saved to: %s\n', saveDir);

%% 6. Summary
fprintf('\n=== Summary ===\n');
for i = 1:4
    lbl = [group_labels{i} ' ' condition_labels{i}];
    p   = stats_res.p(i);

    if p < 0.001
        ps = 'P<0.001***';
    elseif p < 0.01
        ps = sprintf('P=%.3f**', p);
    elseif p < 0.05
        ps = sprintf('P=%.3f*', p);
    else
        ps = sprintf('P=%.3f', p);
    end

    fprintf('%-12s n=%d  mean=%.3f  t(%d)=%.2f  %s\n', ...
        lbl, stats_res.n(i), stats_res.mean(i), stats_res.df(i), stats_res.t(i), ps);
end
fprintf('\n*P<0.05  **P<0.01  ***P<0.001\n');

%% ========================================================================
%% Helper Functions
%% ========================================================================

function rho = calculate_rsa_mixed_model(X, p_x_raw, p_y_raw, name_x, name_y)
    name_x = lower(name_x);
    name_y = lower(name_y);

    if (contains(name_x,'freq') || contains(name_x,'cpd')) && ~contains(name_x,'ori')
        p_x_model = log1p(p_x_raw);
    else
        p_x_model = p_x_raw;
    end

    if (contains(name_y,'freq') || contains(name_y,'cpd')) && ~contains(name_y,'ori')
        p_y_model = log1p(p_y_raw);
    else
        p_y_model = p_y_raw;
    end

    p_x_norm = (p_x_model - min(p_x_model)) / (range(p_x_model) + eps);
    p_y_norm = (p_y_model - min(p_y_model)) / (range(p_y_model) + eps);

    Model_RDM  = squareform(pdist([p_x_norm, p_y_norm], 'euclidean'));
    Neural_RDM = squareform(pdist(X, 'correlation'));

    mask = triu(true(size(Model_RDM)), 1);
    rho  = corr(Model_RDM(mask), Neural_RDM(mask), 'type', 'Spearman');
end

function [X_combined, p_x_raw, p_y_raw] = stitch_data(FR_by_category, P, indices, p_x_name, p_y_name)
    X_combined = [];
    p_x_raw = [];
    p_y_raw = [];

    table_names = {'G','H','R','S','T'};

    for c = indices
        if c > length(FR_by_category)
            continue;
        end

        X_c = mean(FR_by_category{c}, 3)';
        if isempty(X_c) || size(X_c,1) < 2
            continue;
        end

        T_params = P.(table_names{c});
        if ~ismember(p_x_name, T_params.Properties.VariableNames) || ...
           ~ismember(p_y_name, T_params.Properties.VariableNames)
            continue;
        end

        X_combined = [X_combined; X_c];
        p_x_raw    = [p_x_raw; T_params.(p_x_name)];
        p_y_raw    = [p_y_raw; T_params.(p_y_name)];
    end
end

function [P, scheme, cfg] = getPara(FR_by_category, scheme_in)
    if nargin < 2
        scheme_in = 'auto';
    end

    nStim = cellfun(@(X) size(X,2), FR_by_category(:)');

    if strcmpi(scheme_in,'auto')
        if isequal(nStim,[36 12 8 12 4])
            scheme = 'p72';
        elseif isequal(nStim,[48 24 6 72 6])
            scheme = 'p159';
        else
            scheme = 'p159';
            if nStim(4) < 40 && nStim(1) < 40
                scheme = 'p72';
            end
        end
    else
        scheme = scheme_in;
    end

    switch lower(scheme)
        case 'p72'
            cfg = cfg_72();
        case 'p159'
            cfg = cfg_159();
        otherwise
            error('Unknown scheme: %s', scheme);
    end

    P.G = build_table_block(nStim(1), "G", ...
        struct('freq_cpd',cfg.G.freq,'ori_deg',cfg.G.ori_deg), cfg.G.order);

    P.H = build_table_block(nStim(2), "H", ...
        struct('radial_freq',cfg.H.radial_freq,'ori_deg',cfg.H.ori_deg), cfg.H.order);

    P.R = table((1:nStim(3))', repmat("R",nStim(3),1), cfg.R.ang_freq(:), ...
        'VariableNames', {'stim_idx','class','ang_freq'});
    P.R.class_local_idx = P.R.stim_idx;
    P.R.radial_freq = zeros(nStim(3),1);

    fieldsS = struct('ang_freq',cfg.S.ang_freq,'radial_freq',cfg.S.radial_freq);
    if isfield(cfg.S,'hand')
        fieldsS.hand = cfg.S.hand;
    end

    P.S = build_table_block(nStim(4), "S", fieldsS, cfg.S.order);

    P.T = table((1:nStim(5))', repmat("T",nStim(5),1), cfg.T.radial_freq(:), ...
        'VariableNames', {'stim_idx','class','radial_freq'});
    P.T.class_local_idx = P.T.stim_idx;
    P.T.ang_freq = zeros(nStim(5),1);

    tabs = standardizeVarOrder({P.G,P.H,P.R,P.S,P.T});
    P.all = vertcat(tabs{:});
end

function cfg = cfg_159()
    six_list = [1 2 4 8 12 16];

    cfg.S.radial_freq = six_list;
    cfg.S.ang_freq    = six_list;
    cfg.S.hand        = ["L","R"];
    cfg.S.order       = 'hand->ang_freq->radial_freq';

    cfg.R.ang_freq    = six_list;
    cfg.T.radial_freq = six_list;

    cfg.H.radial_freq = six_list;
    cfg.H.ori_deg     = [0 45 90 135];
    cfg.H.order       = 'radial_freq->ori_deg';

    cfg.G.freq        = six_list;
    cfg.G.ori_deg     = (0:7)*(180/8);
    cfg.G.order       = 'freq_cpd->ori_deg';
end

function cfg = cfg_72()
    cfg.G.ori_deg = 0:15:165;
    cfg.G.freq    = [4 8 16];
    cfg.G.order   = 'ori_deg->freq_cpd';

    cfg.H.radial_freq = [4 6 8];
    cfg.H.ori_deg     = [0 45 90 135];
    cfg.H.order       = 'radial_freq->ori_deg';

    cfg.R.ang_freq = [1 1 1 1 2 4 5 8];

    cfg.S.ang_freq    = [2 4 6 8];
    cfg.S.radial_freq = [2 4 6];
    cfg.S.order       = 'ang_freq->radial_freq';

    cfg.T.radial_freq = [4 6 8 10];
end

function T = build_table_block(nStim, classTag, fieldsStruct, orderStr)
    names = fieldnames(fieldsStruct);
    order = string(split(strrep(orderStr,' ',''),'->'))';

    lin = 0:(nStim-1);
    idxMap = struct();
    stride = 1;

    for k = 1:numel(order)
        nm = order{k};
        Lk = numel(fieldsStruct.(nm));
        idxMap.(nm) = mod(floor(lin/stride), Lk) + 1;
        stride = stride * Lk;
    end

    T = table((1:nStim)', repmat(string(classTag),nStim,1), ...
        'VariableNames', {'stim_idx','class'});

    for i = 1:numel(names)
        nm = names{i};
        vec = fieldsStruct.(nm);
        ix  = idxMap.(string(nm));
        T.(nm) = vec(ix)';
    end

    T.class_local_idx = T.stim_idx;
end

function tabs = standardizeVarOrder(tablist)
    pref = {'class','stim_idx','class_local_idx','freq_cpd','ori_deg','radial_freq','ang_freq','hand'};
    hasHand = any(cellfun(@(T) any(strcmp('hand',T.Properties.VariableNames)), tablist));
    tabs = cell(size(tablist));

    for i = 1:numel(tablist)
        Ti = tablist{i};

        for j = 1:numel(pref)
            nm = pref{j};
            if ~ismember(nm, Ti.Properties.VariableNames)
                if strcmp(nm,'hand') && hasHand
                    Ti.(nm) = strings(height(Ti),1);
                else
                    Ti.(nm) = NaN(height(Ti),1);
                end
            end
        end

        keep = pref(ismember(pref, Ti.Properties.VariableNames));
        rest = setdiff(Ti.Properties.VariableNames, keep, 'stable');
        Ti   = Ti(:, [keep, rest]);
        tabs{i} = Ti;
    end
end