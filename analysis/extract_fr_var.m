function extract_fr_var(sessions, CDT_ROOT, outRoot, neuron_selector)
% Objective:
% 1. Extract FR, format [N x C x R] (neuron x condition x trial),
%    which is the same core logic as FR_s7.
% 2. Group by G/H/R/S/T categories, "slicing" the [N x C x R] matrix into 5 groups.
% 3. Output a Cell Array (5x1),
%    where cell{1} is an [N x C_g x R] matrix (G category)
%    where cell{2} is an [N x C_h x R] matrix (H category)
%    ...etc.
if nargin < 1, sessions = 7; end
% neuron_selector(session, filename) must return the original 2-by-N
% channel/unit list. No SNR or area filtering is applied here.
assert(nargin == 4, 'Supply sessions, CDT directory, output directory, and neuron selector.');
% Paths

% *** New output path ***

if ~exist(outRoot,'dir'), mkdir(outRoot); end
% Stimulus end time (ms) for each session
stim_end_ms_map = containers.Map({1,2,4,5,6,7},{870,870,870,1000,1000,1000});
% --- Time windows (ms) ---
STIM_ON_MS      = 630;
BASELINE_ON_MS  = 350; BASELINE_OFF_MS = 500; 
BASELINE_LEN_MS = BASELINE_OFF_MS - BASELINE_ON_MS;
% --- ADDED: Define mapping from category to Condition index ---
category_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
% s1 // 4 * 45 //S T R H
category_map(1) = struct('T1', 1:45, 'U1', 46:90, 'U2', 91:135, 'T2', 136:180);
% s2 // 4 * 45 //S T R H
category_map(2) = category_map(1);
% s4 // 4 * 27 // 4S
category_map(4) = struct('T1', 1:27, 'U1', 28:54, 'T2', 55:81, 'U2', 82:108);
% s5 // 4 * 27 // ?
category_map(5) = struct('T1', 1:27, 'T2', 28:54, 'U1', 55:81, 'U2', 82:108);
% s6 // 4 * 27 // 4S
category_map(6) = struct('T1', 1:27, 'U1', 28:54, 'U2', 55:81, 'T2', 82:108);
% s7 // 8 * 47 // H R S T
category_map(7) = struct('T1', 1:45, 'U1', 46:90, 'T2', 91:135, 'U2', 136:180, 'T3',181:225,'U3',226:270,'T4',271:315,'U4',316:360);
% ------------------------------------------------

% --- MODIFIED LINE: Updated cat_names to match struct fields ---
cat_names = {'T1', 'T2', 'T3', 'T4', 'U1', 'U2', 'U3', 'U4'};
num_categories = numel(cat_names);
for si = 1:numel(sessions)
    sess = sessions(si);
    % Check if session is in our stimulus and category definitions
    if ~isKey(stim_end_ms_map, sess) || ~isKey(category_map, sess)
        fprintf('Skipping session %d (no stim_end or category_map entry).\n', sess);
        continue; 
    end
    
    STIM_OFF_MS = stim_end_ms_map(sess);
    WIN_LEN_MS  = STIM_OFF_MS - STIM_ON_MS; 
    if WIN_LEN_MS <= 0, continue; end
    
    % Get category definitions for the current session
    categories = category_map(sess);
    
    dataDir = fullfile(CDT_ROOT, sprintf('s%d',sess), 'var');
    if ~exist(dataDir,'dir'), continue; end
    
    M = dir(fullfile(dataDir,'*.mat'));
    for fi = 1:numel(M)
        fpath = fullfile(dataDir, M(fi).name);
        S = load(fpath,'CDTTables');
        if ~isfield(S,'CDTTables') || isempty(S.CDTTables), continue; end
        
        cdt = S.CDTTables{1};
        cond_vec = cdt.condition(:);
        if isempty(cond_vec) || ~isnumeric(cond_vec), continue; end
        
        maxc = max(cond_vec); if isempty(maxc) || maxc<=0, continue; end
        
        % s4 only keeps the first 108 conditions (from FR_s7)
        if sess==4
            kept_old_ids = 1:min(108, maxc);
        else
            kept_old_ids = 1:maxc;
        end
        cnd_num = numel(kept_old_ids); % This is C (total conditions)
        
        % Trial count per condition, use the minimum (from FR_s7)
        counts = zeros(cnd_num,1);
        for j=1:cnd_num, counts(j) = sum(cond_vec==kept_old_ids(j)); end
        have = counts(counts>0); if isempty(have), continue; end
        num_trials = min(have); % This is R (trials)
        
        % Select trial indices to be used for each condition (from FR_s7)
        trials_idx_by_cond = cell(cnd_num,1);
        for j=1:cnd_num
            idx = find(cond_vec==kept_old_ids(j));
            trials_idx_by_cond{j} = idx(1:min(num_trials,numel(idx)));
        end
        
        % Get neuron list (V2 only) (from FR_s7)
        this_date = M(fi).name;
        Neurons = neuron_selector(sess, this_date);
        if isempty(Neurons), continue; end
        if size(Neurons,1)==1, Neurons=[Neurons; zeros(1,size(Neurons,2))]; end
        
        pairs = Neurons.';
        num_neurons = size(pairs,1); % This is N (neurons)
        
        % Pre-allocate
        % *** This is the original FR_mat from FR_s7 ***
        FR_mat = zeros(num_neurons, cnd_num, num_trials, 'single');
        
        baseline_per_neuron = zeros(num_neurons,1);
        crop_on  = STIM_ON_MS/1000; crop_off = STIM_OFF_MS/1000;
        
        % -- Pre-build spikeTimes access handles by trial (from FR_s7)
        ele_all = cdt.spikeElectrode; 
        uni_all = cdt.spikeUnit;      
        tim_all = cdt.spikeTimes;     
        
        % baseline (350–500 ms) (from FR_s7)
        for ni=1:num_neurons
            ele = pairs(ni,1); ut = pairs(ni,2);
            acc = [];
            for j=1:cnd_num
                tr_use = trials_idx_by_cond{j};
                for r = tr_use(:)'
                    st_cell = tim_all{r}(ele_all{r}==ele & uni_all{r}==ut);
                    if iscell(st_cell), sp = horzcat(st_cell{:}); else, sp = st_cell; end
                    cnt_b = nnz(sp > BASELINE_ON_MS/1000 & sp < BASELINE_OFF_MS/1000);
                    acc(end+1,1) = cnt_b / BASELINE_LEN_MS * 1000;
                end
            end
            if ~isempty(acc), baseline_per_neuron(ni) = mean(acc); end
        end
        
        % FR Main Loop (from FR_s7)
        for ni=1:num_neurons
            ele = pairs(ni,1); ut = pairs(ni,2);
            base = baseline_per_neuron(ni);
            for j=1:cnd_num % Iterate over all conditions
                tr_use = trials_idx_by_cond{j};
                for rr=1:num_trials % Iterate over all trials
                    r = tr_use(rr);
                    st_cell = tim_all{r}(ele_all{r}==ele & uni_all{r}==ut);
                    if iscell(st_cell), sp = horzcat(st_cell{:}); else, sp = st_cell; end
                    cnt = nnz(sp > crop_on & sp < crop_off);
                    FR_mat(ni,j,rr) = single(cnt / WIN_LEN_MS * 1000 - base);
                end
            end
        end
        
        % --- ADDED: Slice FR_mat by category ---
        
        % FR_by_category will be a 8x1 cell array
        FR_by_category = cell(num_categories, 1);
        % id_to_matrix_idx maps each "kept_old_id" (e.g., 1 to 108) to
        % the index of the second dimension of FR_mat (also 1 to 108)
        % (This is a simple 1-to-1 mapping, because kept_old_ids is just 1:cnd_num)
        id_to_matrix_idx = containers.Map(kept_old_ids, 1:cnd_num);
        
        for ci = 1:num_categories % Iterate through '1'...'8'
            cat_name = cat_names{ci};
            
            % **********************************
            % Additional fix
            if ~isfield(categories, cat_name)
                % This category (e.g., '5') does not exist for this session (e.g., s4).
                % FR_by_category{ci} will remain empty, and we skip to the next.
                continue;
            end
            % **********************************
            
            % 1. Find the conditions included in this category (G/H/..) that *actually exist*
            conds_in_this_cat_raw = getfield(categories, cat_name);
            % *Must* intersect with the original kept_old_ids
            conds_to_use = intersect(conds_in_this_cat_raw, kept_old_ids);
            
            if isempty(conds_to_use)
                fprintf('  Sess %d, File %s: Category %s has no valid conditions. Skipping.\n', sess, M(fi).name, cat_name);
                continue; % Data for this category will be an empty cell
            end
            
            % 2. Find the indices of these conditions in the FR_mat matrix
            matrix_indices = cell2mat(values(id_to_matrix_idx, num2cell(conds_to_use)));
            
            % 3. "Slice"
            % From FR_mat (N x C_total x R),
            % extract all N, all R, but only the corresponding indices in the C dimension
            FR_by_category{ci} = FR_mat(:, matrix_indices, :);
        end
        
        % --- Save (Modified) ---
        tag6 = M(fi).name(1:min(6,numel(M(fi).name)));
        tag  = sprintf('s%d_var_%s', sess, tag6);
        % Note: save path has been changed (outRoot)
        save_path = fullfile(outRoot, ['FR_' tag '.mat']);
        
        % Update meta information
        meta = struct();
        meta.session = sess; meta.file = M(fi).name;
        meta.num_neurons = num_neurons;
        meta.num_trials_balanced = num_trials;
        meta.category_names = cat_names; % Category names {'1', '2', ... '8'}
        meta.category_cond_idx_original = categories; % Original condition indices
        meta.stim_on_ms = STIM_ON_MS; meta.stim_off_ms = STIM_OFF_MS;
        meta.baseline_on_ms = BASELINE_ON_MS; meta.baseline_off_ms = BASELINE_OFF_MS;
        meta.neurons_used = Neurons; 
        
        % Save FR_by_category, this is the N groups of 3D matrices you wanted
        save(save_path, 'FR_by_category', 'meta', 'baseline_per_neuron', '-v7.3');
        fprintf('Saved FR by Category: %s\n', save_path);
        
    end % File loop (fi)
end % Session loop (si)
fprintf('Done.\n');
end
