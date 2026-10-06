function calculate_global_sdi(input_directory, output_directory)
% Compute saved metrics from firing rates without neuron reselection.
%% ---------------------------------------------------------------------
% BATCH Pairwise SDI Calculation Script (FOR "GLOBAL" DATA)
%
% Updated logic:
% 1. Load FR_by_category as a 1-by-5 or 5-by-1 cell array.
% 2. The five cells correspond to G/H/R/S/T classes.
% 3. Compute the 5-by-5 matrix of all category pairs, including G-H, G-R and H-R.
%
% Save:
%   1. SDI_matrix (5x5)
%   2. Euclidean_Dist_matrix: 5-by-5 category-center distances.
%   3. Mean_Class_Variance: one noise value per class, as a 1-by-5 vector.
% ---------------------------------------------------------------------

clc;
%% 1. Session-specific configuration removed.
% It is unnecessary for comparisons between G/H/R and other categories.

%% --- 2. Parameters & Setup ---
% Use the global data path.



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
    
    if ~iscell(FR_by_category)
        error('FATAL: FR_by_category is not a cell array. This script is for "Global" data.');
    end
    
    % 3. Precompute the mean Mu and variance Var for each class.
    num_categories = numel(FR_by_category); % Expected value: 5.
    if num_categories < 2
        warning('  -> Less than 2 categories found in cell array. Skipping.');
        continue;
    end
    
    fprintf('  Found %d categories (classes).\n', num_categories);
    
    % Get N_neurons from meta.
    if isfield(original_meta, 'N_neurons')
        N_neurons = original_meta.N_neurons;
    else
        % Detect the dimensions using the first nonempty cell.
        first_valid_cell = find(~cellfun(@isempty, FR_by_category), 1);
        if isempty(first_valid_cell)
            warning('  -> All categories are empty. Skipping.');
            continue;
        end
        N_neurons = size(FR_by_category{first_valid_cell}, 1);
    end
    
    mu_all_classes = nan(N_neurons, num_categories);
    var_all_classes = nan(N_neurons, num_categories);
    mean_var_per_class = nan(1, num_categories); % Store the mean noise for each class.
    
    fprintf('  Pre-calculating Mean (Signal) and Variance (Noise) for %d classes...\n', num_categories);
    
    for c_idx = 1:num_categories
        FR_cat = FR_by_category{c_idx};
        if isempty(FR_cat)
             fprintf('    -> Class %d is empty. Skipping.\n', c_idx);
             continue;
        end
        
        % FR_cat has dimensions N_neurons x N_Stimuli x N_repetitions.
        % The example shows N_rep ranging from 1 to 12 and N_stim from 4 to 36.
        
        % 1. Compute the mean signal.
        % Average across all stimuli and repetitions.
        mu_all_classes(:, c_idx) = mean(FR_cat, [2, 3], 'omitnan');
        
        % 2. Compute the noise variance.
        % Compute trial variance within each stimulus, then average across stimuli.
        var_per_stim = var(FR_cat, 0, 3, 'omitnan'); % [N_neur x N_stim]
        var_all_classes(:, c_idx) = mean(var_per_stim, 2, 'omitnan'); % [N_neur x 1]
        
        % Store a single class noise value averaged across neurons.
        mean_var_per_class(c_idx) = mean(var_all_classes(:, c_idx), 'all', 'omitnan');
        
        fprintf('    - Class %d: Mean Var = %.4f\n', c_idx, mean_var_per_class(c_idx));
    end
    
    % 4. Compute pairwise SDI and distances between classes.
    num_pairs = nchoosek(num_categories, 2);
    fprintf('  Calculating Pairwise SDI & Dist for %d class pairs...\n', num_pairs);
    
    SDI_matrix = nan(num_categories, num_categories);
    Euclidean_Dist_matrix = nan(num_categories, num_categories);
    
    for c_i = 1:num_categories
        for c_j = (c_i + 1):num_categories
            
            mu_i = mu_all_classes(:, c_i);
            mu_j = mu_all_classes(:, c_j);
            
            var_i_per_neuron = var_all_classes(:, c_i);
            var_j_per_neuron = var_all_classes(:, c_j);
            
            % Check data validity.
            if any(isnan(mu_i)) || any(isnan(mu_j)) || any(isnan(var_i_per_neuron)) || any(isnan(var_j_per_neuron))
                fprintf('    -> Skipping pair (%d, %d) due to NaN values in mu or var.\n', c_i, c_j);
                continue;
            end
            
            % 1. Compute Euclidean distance (signal).
            euc_dist_value = norm(mu_i - mu_j);
            
            % 2. Compute SDI (signal/noise).
            pooled_var = (var_i_per_neuron + var_j_per_neuron) / 2;
            pooled_var(pooled_var < eps) = eps;
            
            diff_normalized = (mu_i - mu_j).^2 ./ pooled_var;
            SDI_value = sqrt(sum(diff_normalized));
            
            % Store SDI.
            SDI_matrix(c_i, c_j) = SDI_value;
            SDI_matrix(c_j, c_i) = SDI_value; 
            
            % Store distances.
            Euclidean_Dist_matrix(c_i, c_j) = euc_dist_value;
            Euclidean_Dist_matrix(c_j, c_i) = euc_dist_value;
        end
    end
    
    fprintf('    Processed all %d pairs\n', num_pairs);

    % 5. Save results.
    SDI_results = struct();
    
    % Rename SDI_results.Global.
    SDI_results.Global.SDI_matrix = SDI_matrix;
    SDI_results.Global.Euclidean_Dist_matrix = Euclidean_Dist_matrix; % Category-center distance.
    SDI_results.Global.Mean_Class_Variance = mean_var_per_class; % One noise value per class, as a 1-by-5 vector.
    SDI_results.Global.N_neurons = N_neurons;
    SDI_results.Global.N_categories = num_categories;

    % Print a summary.
    fprintf('  Mean Pairwise SDI (all classes) = %.4f\n', mean(SDI_matrix(triu(true(num_categories),1)), 'omitnan'));
    fprintf('  Mean Pairwise Dist (all classes) = %.4f\n', mean(Euclidean_Dist_matrix(triu(true(num_categories),1)), 'omitnan'));
        
    save_name = ['SDI_Results_' fileName];
    save_path = fullfile(output_directory, save_name);
    
    fprintf('\nSaving SDI results to: %s\n', save_path);
    % Save the updated SDI_results structure.
    save(save_path, 'SDI_results', 'original_meta', '-v7.3');
    
    fprintf('File processing complete.\n');
    
end % end file loop
fprintf('\n\n======================================================\n');
fprintf('All files processed. SDI calculation finished.\n');
fprintf('======================================================\n');
end
