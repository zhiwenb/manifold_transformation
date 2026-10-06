function [DiD, table_out, out_mat_file] = calculate_manifold_did( ...
    session_files, controls_per_session, n_baseline_days, after_idx_all, out_mat_file)
    if nargin < 4, after_idx_all = []; end
    if nargin < 5, out_mat_file = ''; end
    if nargin < 3 || isempty(n_baseline_days)
        error('Provide one cell entry per session for baseline and after indices.');
    end
    if ~iscell(n_baseline_days)
        error('Provide one cell entry per session for baseline and after indices.');
    end
    if numel(n_baseline_days) ~= numel(session_files)
        error('Provide one cell entry per session for baseline and after indices.');
    end
    baseline_idx_all = n_baseline_days; 
    
    if isempty(after_idx_all)
        after_idx_all = cell(numel(session_files), 1);
    else
        if numel(after_idx_all) ~= numel(session_files)
            error('Provide one cell entry per session for baseline and after indices.');
        end
    end
    
    rows = {};
    DiD = struct('session_name', {}, 'control_set_idx', {}, 'control_set_name', {}, ...
                 'trained_cat', {}, 'trained_idx', {}, 'n_days_total', {}, ...
                 'baseline_idx', {}, 'after_idx', {}, ...
                 'alphaM_base_mean', {}, 'alphaM_after_mean', {}, 'alphaM_DiD', {}, ...
                 'RM_base_mean', {},    'RM_after_mean',    {}, 'RM_DiD',    {}, ...
                 'DM_base_mean', {},    'DM_after_mean',    {}, 'DM_DiD',    {} );
    for s = 1:numel(session_files)
        session_path = session_files{s};
        ctl_idx = controls_per_session{s};
        if ~isnumeric(ctl_idx) || isempty(ctl_idx)
            warning('Invalid or missing session data; skipping this comparison.');
            continue;
        end
        ctl_idx = unique(ctl_idx(:)');   
        S = load(session_path);
        if ~isfield(S,'session_results')
            warning('No session_results in %s. Skip.', session_path);
            continue;
        end
        SR = S.session_results;
        n_files = numel(SR.files);
        if n_files < 2
            warning('Session %s: days < 2, skip.', session_path);
            continue;
        end
        num_per_day = zeros(n_files,1);
        for d = 1:n_files
            if isfield(SR.files(d),'manifold_details') && ~isempty(SR.files(d).manifold_details)
                num_per_day(d) = numel(SR.files(d).manifold_details);
            end
        end
        NC = max(num_per_day);
        if NC < 1
            warning('Invalid or missing session data; skipping this comparison.');
            continue;
        end
        a_M_mat = nan(n_files, NC);
        R_M_mat = nan(n_files, NC);
        D_M_mat = nan(n_files, NC);
        for d = 1:n_files
            if ~isfield(SR.files(d),'manifold_details') || isempty(SR.files(d).manifold_details)
                continue;
            end
            det = SR.files(d).manifold_details;
            M = min(numel(det), NC);
            for m = 1:M
                a_M_mat(d,m) = det(m).a_M;
                R_M_mat(d,m) = det(m).R_M;
                D_M_mat(d,m) = det(m).D_M;
            end
        end
        
        base_idx = baseline_idx_all{s}; % e.g., [2] or [1, 2]
        aft_idx = after_idx_all{s};     % e.g., [3, 4]

        base_idx = base_idx(base_idx >= 1 & base_idx <= n_files);
        base_idx = unique(base_idx(:)');
        if isempty(base_idx)
            warning('Invalid or missing session data; skipping this comparison.');
            continue;
        end
        
        if isempty(aft_idx)
            aft_idx = setdiff(1:n_files, base_idx);
        else
            aft_idx = aft_idx(aft_idx >= 1 & aft_idx <= n_files);
            aft_idx = setdiff(unique(aft_idx(:)'), base_idx);
        end
        
        if isempty(aft_idx)
            warning('Invalid or missing session data; skipping this comparison.');
            continue;
        end
        
        ctl_idx = ctl_idx(ctl_idx >= 1 & ctl_idx <= NC);
        if isempty(ctl_idx)
            warning('Invalid or missing session data; skipping this comparison.');
            continue;
        end
        
        trained_idx = setdiff(1:NC, ctl_idx);
        if isempty(trained_idx)
            warning('Invalid or missing session data; skipping this comparison.');
            continue;
        end
        
        session_name_str = infer_session_name(session_path);
        is_s7_special_case = SR.stimulus_definition == 7;
        is_paired_case = is_s7_special_case && (numel(trained_idx) == numel(ctl_idx));

        if (is_paired_case)
            fprintf('Paired trained/control comparisons: %d pairs.\n', numel(trained_idx));
            
            for i = 1:numel(trained_idx)
                c_idx = trained_idx(i); % Trained idx (e.g., 5)
                u_idx = ctl_idx(i);     % Control idx (e.g., 1)
                
                a_ctl_daymean = a_M_mat(:, u_idx);
                R_ctl_daymean = R_M_mat(:, u_idx);
                D_ctl_daymean = D_M_mat(:, u_idx);
                
                [DiD, rows] = calculate_did_for_pair(DiD, rows, ...
                    a_M_mat, R_M_mat, D_M_mat, ...
                    a_ctl_daymean, R_ctl_daymean, D_ctl_daymean, ...
                    c_idx, u_idx, sprintf('paired_U%d',u_idx), session_path, ...
                    base_idx, aft_idx, n_files);
            end
            
        else
            fprintf('  -> (Standard Mode): T vs. mean(U) (T_N=%d, U_N=%d).\n', numel(trained_idx), numel(ctl_idx));
            
            a_ctl_daymean = nanmean(a_M_mat(:, ctl_idx), 2);
            R_ctl_daymean = nanmean(R_M_mat(:, ctl_idx), 2);
            D_ctl_daymean = nanmean(D_M_mat(:, ctl_idx), 2);
            
            ctl_name = 'control_mean';
            
            for c_idx = trained_idx
                [DiD, rows] = calculate_did_for_pair(DiD, rows, ...
                    a_M_mat, R_M_mat, D_M_mat, ...
                    a_ctl_daymean, R_ctl_daymean, D_ctl_daymean, ...
                    c_idx, ctl_idx, ctl_name, session_path, ...
                    base_idx, aft_idx, n_files);
            end
        end
    end 
    
    if ~isempty(rows)
        table_out = cell2table(rows, 'VariableNames', { ...
            'Session','ControlSet','TrainedCat','N_base','N_after', ...
            'alphaM_base','alphaM_after','alphaM_DiD', ...
            'RM_base','RM_after','RM_DiD', ...
            'DM_base','DM_after','DM_DiD'});
    else
        table_out = table();
    end
    if ~isempty(out_mat_file)
        summary = struct();
        summary.generated_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        summary.notes = ['diff(day)=trained - control; Δ = mean(after)-mean(base). ', ...
                         'Paired comparisons for stimulus definition 7; otherwise trained minus control mean.'];
        save(out_mat_file, 'DiD', 'table_out', 'summary', 'session_files', 'controls_per_session', ...
             'n_baseline_days', 'after_idx_all');
    else
        out_mat_file = '';
    end
end

function [DiD, rows] = calculate_did_for_pair(DiD, rows, ...
    a_M_mat, R_M_mat, D_M_mat, ...
    a_ctl_daymean, R_ctl_daymean, D_ctl_daymean, ...
    c_idx, ctl_idx_info, ctl_name, session_path, ...
    base_idx, aft_idx, n_files)

    a_diff = a_M_mat(:, c_idx) - a_ctl_daymean;
    R_diff = R_M_mat(:, c_idx) - R_ctl_daymean;
    D_diff = D_M_mat(:, c_idx) - D_ctl_daymean;
    
    base_a  = base_idx(isfinite(a_diff(base_idx)));
    after_a = aft_idx(isfinite(a_diff(aft_idx)));
    base_R  = base_idx(isfinite(R_diff(base_idx)));
    after_R = aft_idx(isfinite(R_diff(aft_idx)));
    base_D  = base_idx(isfinite(D_diff(base_idx)));
    after_D = aft_idx(isfinite(D_diff(aft_idx)));
    
    if isempty(base_a) || isempty(after_a) || isempty(base_R) || isempty(after_R) ...
                       || isempty(base_D) || isempty(after_D)
        warning('Invalid or missing session data; skipping this comparison.');
        return;
    end
    
    a_base_mean  = mean(a_diff(base_a));
    a_after_mean = mean(a_diff(after_a));
    R_base_mean  = mean(R_diff(base_R));
    R_after_mean = mean(R_diff(after_R));
    D_base_mean  = mean(D_diff(base_D));
    D_after_mean = mean(D_diff(after_D));
    
    did_alpha = a_after_mean - a_base_mean;
    did_R     = R_after_mean - R_base_mean;
    did_D     = D_after_mean - D_base_mean;
    
    rec.session_name       = infer_session_name(session_path);
    rec.control_set_idx    = ctl_idx_info; % S5:[1,2], S7(T1):[1]
    rec.control_set_name   = ctl_name;     % S5:'control_mean', S7(T1):'paired_U1'
    rec.trained_cat        = 'trained';    
    rec.trained_idx        = c_idx;        
    rec.n_days_total       = n_files;
    rec.baseline_idx       = base_idx;
    rec.after_idx          = aft_idx;
    rec.alphaM_base_mean   = a_base_mean;
    rec.alphaM_after_mean  = a_after_mean;
    rec.alphaM_DiD         = did_alpha;
    rec.RM_base_mean       = R_base_mean;
    rec.RM_after_mean      = R_after_mean;
    rec.RM_DiD             = did_R;
    rec.DM_base_mean       = D_base_mean;
    rec.DM_after_mean      = D_after_mean;
    rec.DM_DiD             = did_D;
    
    DiD(end+1) = rec; %#ok<AGROW>
    
    rows(end+1,:) = { ...
        string(rec.session_name), string(rec.control_set_name), ...
        string(rec.trained_cat), numel(base_idx), numel(aft_idx), ...
        a_base_mean, a_after_mean, did_alpha, ...
        R_base_mean, R_after_mean, did_R, ...
        D_base_mean, D_after_mean, did_D }; %#ok<AGROW>
end

function name = infer_session_name(pth)
    [~, base, ~] = fileparts(pth);
    name = base;
end
