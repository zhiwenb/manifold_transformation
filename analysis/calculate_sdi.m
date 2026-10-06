function calculate_sdi(input_directory, output_directory)
% Compute saved metrics from firing rates without neuron reselection.
%% ---------------------------------------------------------------------
% BATCH Pairwise SDI Calculation Script (Between Conditions)
%
% Final output fields:
%   1. SDI (signal/noise) for all pairs.
%   2. T/U trial variance (noise).
%   3. T/U interstimulus distance (signal).
% ---------------------------------------------------------------------

clc;
%% 1. Session settings must match the analysis script.
session_config = struct();
session_config.s1.trained_ids = [7;8;9;12;13;14;17;18;19];
session_config.s1.control_ids = [27;28;29;32;33;34;37;38;39];
session_config.s2.trained_ids = [7;8;9;12;13;14;17;18;19];
session_config.s2.control_ids = [27;28;29;32;33;34;37;38;39];
session_config.s7.trained_ids = [7;8;9;12;13;14;17;18;19];
session_config.s7.control_ids = [27;28;29;32;33;34;37;38;39];
session_config.s4.trained_ids = [3;11;12;13;21];
session_config.s4.control_ids = [7;15;16;17;25];
session_config.s6.trained_ids = [3;11;12;13;21];
session_config.s6.control_ids = [7;15;16;17;25];
session_config.s5.trained_ids = [10:18];
session_config.s5.control_ids = [1:9]';


%% --- 2. Parameters & Setup ---



% Ensure that the output directory exists.
if ~exist(output_directory, 'dir')
    fprintf('Creating output directory: %s\n', output_directory);
    mkdir(output_directory);
end
%% 3. Find all MAT files.
fprintf('Searching for .mat files in: %s\n', input_directory);
mat_files = dir(fullfile(input_directory, '*.mat'));
if isempty(mat_files)
    error('No .mat files found in the specified directory: %s', input_directory);
end
fprintf('Found %d .mat files. Results will be saved to: %s\n', ...
    numel(mat_files), output_directory);
%% Start the file loop.
for f_idx = 1:numel(mat_files)
    
    fileName = mat_files(f_idx).name;
    filePath = fullfile(input_directory, fileName);
    
    fprintf('\n\n======================================================\n');
    fprintf('Processing file (%d/%d): %s\n', f_idx, numel(mat_files), fileName);
    fprintf('======================================================\n');
    
    session_match = regexp(fileName, 's(\d+)', 'tokens');
    if isempty(session_match)
        fprintf('  -> ERROR: Cannot parse session name (e.g., "s4") from file name. Skipping.\n');
        continue;
    end
    session_name = sprintf('s%d', str2double(session_match{1}{1}));
    
    if ~isfield(session_config, session_name)
        fprintf('  -> ERROR: No session_config found for %s. Skipping.\n', session_name);
        continue;
    end
    trained_ids = session_config.(session_name).trained_ids;
    control_ids = session_config.(session_name).control_ids;
    
    % 2. Load data.
    try
        S = load(filePath, 'FR_by_category', 'meta');
    catch e
        fprintf('  -> ERROR loading file: %s. Skipping.\n', e.message);
        continue;
    end
    
    if ~isfield(S, 'FR_by_category') || ~isfield(S, 'meta')
        fprintf('  -> ERROR: File is missing ''FR_by_category'' or ''meta''. Skipping.\n');
        continue;
    end
    
    FR_by_category = S.FR_by_category;
    original_meta = S.meta; 
    category_names = original_meta.category_names;
    num_categories = numel(category_names);
    
    % 3. Initialize the SDI result structure.
    SDI_results = struct();
    
    % 4. Compute results for each category.
    for cat_idx = 1:num_categories
        
        cat_name = category_names{cat_idx};
        fprintf('\n--- Processing Category: %s ---\n', cat_name);
        
        FR_cat = FR_by_category{cat_idx};
        if isempty(FR_cat)
            fprintf('  -> WARNING: Empty data for category %s. Skipping.\n', cat_name);
            continue;
        end
        
        [N_neurons, N_conditions, N_repetitions] = size(FR_cat);
        
        fprintf('  N_neurons = %d, N_conditions = %d, N_repetitions = %d\n', ...
            N_neurons, N_conditions, N_repetitions);
        
        if N_conditions < 2
            fprintf('  -> WARNING: Less than 2 conditions. Skipping calculation.\n');
            continue;
        end
        
        % A. Compute trial variance for trained and control conditions (noise).
        valid_trained_ids = trained_ids(trained_ids <= N_conditions);
        valid_control_ids = control_ids(control_ids <= N_conditions);
        
        FR_T_subset = FR_cat(:, valid_trained_ids, :); 
        var_T_per_cond = var(FR_T_subset, 0, 3); 
        mean_T_var = mean(var_T_per_cond, 'all', 'omitnan');
        
        FR_U_subset = FR_cat(:, valid_control_ids, :); 
        var_U_per_cond = var(FR_U_subset, 0, 3); 
        mean_U_var = mean(var_U_per_cond, 'all', 'omitnan');
        
        fprintf('  Mean "Trained" Trial Variance = %.4f\n', mean_T_var);
        fprintf('  Mean "Control" Trial Variance = %.4f\n', mean_U_var);
        
        
        % B. Compute SDI and distances for condition pairs.
        num_pairs = nchoosek(N_conditions, 2);
        fprintf('  Calculating Pairwise SDI & Dist for %d pairs...\n', num_pairs);
        
        SDI_matrix = nan(N_conditions, N_conditions);
        pair_list = zeros(num_pairs, 2);
        SDI_values = zeros(num_pairs, 1);
        
        Euclidean_Dist_matrix = nan(N_conditions, N_conditions);
        Euclidean_Dist_values = zeros(num_pairs, 1);
        
        pair_count = 0;
        
        for cond_i = 1:N_conditions
            for cond_j = (cond_i + 1):N_conditions
                pair_count = pair_count + 1;
                response_i = squeeze(FR_cat(:, cond_i, :));
                response_j = squeeze(FR_cat(:, cond_j, :));
                
                mu_i = mean(response_i, 2);
                mu_j = mean(response_j, 2);
                
                var_i = var(response_i, 0, 2);
                var_j = var(response_j, 0, 2);
                pooled_var = (var_i + var_j) / 2;
                
                euc_dist_value = norm(mu_i - mu_j);
                
                pooled_var(pooled_var < eps) = eps;
                diff_normalized = (mu_i - mu_j).^2 ./ pooled_var;
                SDI_value = sqrt(sum(diff_normalized));
                
                SDI_matrix(cond_i, cond_j) = SDI_value;
                SDI_matrix(cond_j, cond_i) = SDI_value; 
                pair_list(pair_count, :) = [cond_i, cond_j];
                SDI_values(pair_count) = SDI_value;
                
                Euclidean_Dist_matrix(cond_i, cond_j) = euc_dist_value;
                Euclidean_Dist_matrix(cond_j, cond_i) = euc_dist_value;
                Euclidean_Dist_values(pair_count) = euc_dist_value;
            end
        end
        
        fprintf('    Processed %d/%d pairs\n', pair_count, num_pairs);

        % C. Compute mean interstimulus distance for trained and control conditions (signal).
        fprintf('  (New) Calculating mean distances for T and U subsets...\n');
        
        % Use the helper to extract T-T and U-U pairs from the full matrix.
        trained_dist_values = extract_pair_values(Euclidean_Dist_matrix, valid_trained_ids);
        control_dist_values = extract_pair_values(Euclidean_Dist_matrix, valid_control_ids);
        
        mean_T_dist = mean(trained_dist_values, 'omitnan');
        mean_U_dist = mean(control_dist_values, 'omitnan');
        
        fprintf('  Mean "Trained" Pairwise Dist = %.4f\n', mean_T_dist);
        fprintf('  Mean "Control" Pairwise Dist = %.4f\n', mean_U_dist);
        
        
        % D. Store all results for this category in the structure.
        SDI_results.(cat_name).SDI_matrix = SDI_matrix;
        SDI_results.(cat_name).pair_list = pair_list;
        SDI_results.(cat_name).SDI_values = SDI_values;
        SDI_results.(cat_name).Mean_SDI_AllPairs = mean(SDI_values, 'omitnan');
        
        % Save noise metrics.
        SDI_results.(cat_name).Mean_Trained_Variance = mean_T_var;
        SDI_results.(cat_name).Mean_Control_Variance = mean_U_var;
        
        % Save signal metrics.
        SDI_results.(cat_name).Euclidean_Dist_matrix = Euclidean_Dist_matrix;
        SDI_results.(cat_name).Euclidean_Dist_values = Euclidean_Dist_values;
        SDI_results.(cat_name).Mean_Euclidean_Dist_AllPairs = mean(Euclidean_Dist_values, 'omitnan');
        SDI_results.(cat_name).Mean_Trained_Distance = mean_T_dist;
        SDI_results.(cat_name).Mean_Control_Distance = mean_U_dist;

        SDI_results.(cat_name).N_conditions = N_conditions;
        SDI_results.(cat_name).N_neurons = N_neurons;
        SDI_results.(cat_name).num_pairs = num_pairs;
        
        fprintf('  -> Category %s processing complete.\n', cat_name);
        
    end % end category loop
    
    % 5. Save SDI results for all categories.
    save_name = ['SDI_Results_' fileName];
    save_path = fullfile(output_directory, save_name);
    
    fprintf('\nSaving SDI results to: %s\n', save_path);
    save(save_path, 'SDI_results', 'original_meta', '-v7.3');
    
    fprintf('File processing complete.\n');
    
end % end file loop
fprintf('\n\n======================================================\n');
fprintf('All files processed. SDI calculation finished.\n');
fprintf('======================================================\n');

%% Helper function
function values = extract_pair_values(VALUE_matrix, ids)
% Extract values for the specified IDs from the cond_i-by-cond_j matrix.
    values = [];
    ids = ids(:); 
    if isempty(ids), return; end
    
    % Create all unique ID pairs.
    pairs = nchoosek(ids, 2);
    if isempty(pairs), return; end
    
    for k = 1:size(pairs, 1)
        id_i = pairs(k, 1);
        id_j = pairs(k, 2);
        
        % Check that indices are within the matrix bounds.
        if id_i <= size(VALUE_matrix, 1) && id_j <= size(VALUE_matrix, 2)
            val = VALUE_matrix(id_i, id_j);
            if ~isnan(val)
                values = [values; val];
            end
        elseif id_j <= size(VALUE_matrix, 1) && id_i <= size(VALUE_matrix, 2)
             % Try reversed indices; this is unnecessary for symmetric matrices.
             val = VALUE_matrix(id_j, id_i);
             if ~isnan(val)
                values = [values; val];
             end
        end
    end
end
end
