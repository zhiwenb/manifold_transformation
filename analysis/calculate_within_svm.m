function calculate_within_svm(input_directory, output_directory)
% Compute saved metrics from firing rates without neuron reselection.
%% =====================================================================
%  BATCH Pairwise SVM-LOO Calculation Script
%  Replacement for SDI.m.
%
%  Objective:
%  1. Loop over all MAT files.
%  2. Compute SVM decoding accuracy for every condition pair within each category.
%  3. Use leave-one-out cross-validation.
%  4. Save SVM_Results_*.mat with a structure matching SDI_Results_*.mat.
% =====================================================================

clc;

%% --- 1. Parameters & Setup ---


% Output directory for SVM results.


% Ensure that the output directory exists.
if ~exist(output_directory, 'dir')
    fprintf('Creating output directory: %s\n', output_directory);
    mkdir(output_directory);
end
%% 0. Find all MAT files.
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
    category_names = original_meta.category_names;
    num_categories = numel(category_names);
    
    % 3. Initialize the SVM result structure.
    SVM_results = struct();
    
    % 4. Compute pairwise SVM accuracy between conditions in each category.
    for cat_idx = 1:num_categories
        
        cat_name = category_names{cat_idx};
        fprintf('\n--- Processing Category: %s ---\n', cat_name);
        
        FR_cat = FR_by_category{cat_idx};
        if isempty(FR_cat)
            fprintf('  -> WARNING: Empty data for category %s. Skipping.\n', cat_name);
            continue;
        end
        
        % FR_cat: [N_neurons × N_conditions × N_repetitions]
        [N_neurons, N_conditions, N_repetitions] = size(FR_cat);
        
        fprintf('  N_neurons = %d, N_conditions = %d, N_repetitions = %d\n', ...
            N_neurons, N_conditions, N_repetitions);
        
        if N_conditions < 2
            fprintf('  -> WARNING: Less than 2 conditions. Skipping SVM calculation.\n');
            continue;
        end
        
        num_pairs = nchoosek(N_conditions, 2);
        fprintf('  Total condition pairs: %d\n', num_pairs);
        
        % Initialize result storage.
        SVM_matrix = nan(N_conditions, N_conditions);
        pair_list = zeros(num_pairs, 2); 
        SVM_values = zeros(num_pairs, 1);
        
        pair_count = 0;
        
        % Compute condition-pair SVM accuracy using leave-one-out cross-validation.
        tic_cat = tic;
        for cond_i = 1:N_conditions
            for cond_j = (cond_i + 1):N_conditions
                
                pair_count = pair_count + 1;
                
                % 1. Extract single-trial responses for the two conditions.
                % response_i: [N_neurons × N_repetitions]
                response_i = squeeze(FR_cat(:, cond_i, :));
                response_j = squeeze(FR_cat(:, cond_j, :));
                
                % 2. Prepare SVM features X and labels Y.
                % Transpose to N_trials x N_neurons.
                data_i = response_i'; 
                data_j = response_j';
                
                % Remove NaNs used to pad unequal trial counts.
                data_i(any(isnan(data_i), 2), :) = [];
                data_j(any(isnan(data_j), 2), :) = [];
                
                % Create class labels 1 and -1.
                labels_i =  ones(size(data_i, 1), 1);
                labels_j = -ones(size(data_j, 1), 1);
                
                X = [data_i; data_j]; % [Total_trials x N_neurons]
                Y = [labels_i; labels_j]; % [Total_trials x 1]
                
                N_total_trials = size(X, 1);
                
                if N_total_trials < 2
                    fprintf('    -> Skipping pair (%d, %d): Not enough total trials.\n', cond_i, cond_j);
                    continue;
                end
                
                % 3. Perform leave-one-out cross-validation.
                predictions = zeros(N_total_trials, 1);
                
                for k = 1:N_total_trials
                    % Split training and test sets.
                    X_train = X;
                    X_train(k, :) = [];
                    Y_train = Y;
                    Y_train(k) = [];
                    
                    X_test = X(k, :);
                    
                    % Train a linear SVM.
                    % Use Standardize=true to standardize features.
                    svm_model = fitcsvm(X_train, Y_train, ...
                        'KernelFunction', 'linear', 'Standardize', true);
                    
                    % Predict.
                    predictions(k) = predict(svm_model, X_test);
                end
                
                % 4. Compute accuracy.
                SVM_acc_value = mean(predictions == Y);
                
                % 5. Store results.
                SVM_matrix(cond_i, cond_j) = SVM_acc_value;
                SVM_matrix(cond_j, cond_i) = SVM_acc_value; % Symmetric matrix
                
                pair_list(pair_count, :) = [cond_i, cond_j];
                SVM_values(pair_count) = SVM_acc_value;
                
                if mod(pair_count, 10) == 0 || pair_count == num_pairs
                    fprintf('    Processed %d/%d pairs... (Last acc=%.3f)\n', ...
                        pair_count, num_pairs, SVM_acc_value);
                end
                
            end
        end % end cond_j
        
        toc_cat = toc(tic_cat);
        fprintf('  -> Category %s finished in %.2f sec.\n', cat_name, toc_cat);
        
        % Store category results in the structure.
        SVM_results.(cat_name).SVM_matrix = SVM_matrix;
        SVM_results.(cat_name).pair_list = pair_list;
        SVM_results.(cat_name).SVM_values = SVM_values;
        SVM_results.(cat_name).N_conditions = N_conditions;
        SVM_results.(cat_name).N_neurons = N_neurons;
        SVM_results.(cat_name).num_pairs = num_pairs;
        
        fprintf('  -> Category %s: Mean SVM Acc = %.4f, Std SVM Acc = %.4f\n', ...
            cat_name, mean(SVM_values, 'omitnan'), std(SVM_values, 'omitnan'));
        
    end % end category loop
    
    % 5. Save SVM results for all categories.
    save_name = ['SVM_Results_' fileName]; % New addition
    save_path = fullfile(output_directory, save_name);
    
    fprintf('\nSaving SVM results to: %s\n', save_path);
    save(save_path, 'SVM_results', 'original_meta', '-v7.3'); % New addition
    
    fprintf('File processing complete.\n');
    
end % end file loop

fprintf('\n\n======================================================\n');
fprintf('All files processed. SVM calculation finished.\n');
fprintf('======================================================\n');
end
