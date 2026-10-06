function calculate_internal_svm(input_directory, output_directory)
% Compute saved metrics from firing rates without neuron reselection.
%% ---------------------------------------------------------------------
% Batch internal pairwise SVM decoding, version 5 with fixes.
%
% Objective:
% 1. Loop through all 'FR_...' files in the input folder.
% 2. For each file, load 'FR_by_category'.
% 3. Loop through each category (k = 1 to 5, e.g., 'G', 'H', ...).
% 4.   Inside 'G' (size [N x C_g x R]), loop through all internal
%      condition pairs (i vs j).
% 5.   [FIXED] Use a manual, NON-STRATIFIED cvpartition to handle
%      very low trial counts (R) and avoid warnings.
% 6.   Calculate the average accuracy of all C_g*(C_g-1)/2 pairs.
% 7. Save the 5 average accuracies (G_avg, H_avg, ...) for this file.
% ---------------------------------------------------------------------

clc;

%% --- 1. Parameters & Setup ---
% Directory containing the original FR files.


% Directory for within-category decoding results.

if ~exist(output_directory, 'dir'), mkdir(output_directory); end

% SVM & CV Parameters
CV_FOLDS = 10; % 10-fold CV. 
SVM_KERNEL = 'linear';
RANDOM_SEED = 1; % For reproducible SVM fits and CV partitions.

%% 0. Find all MAT files.
fprintf('Searching for .mat files in: %s\n', input_directory);
mat_files = dir(fullfile(input_directory, '*.mat'));

if isempty(mat_files)
    error('No .mat files found: %s', input_directory);
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
    
    FR_by_category = S.FR_by_category;
    original_meta = S.meta;
    category_names = original_meta.category_names;
    num_categories = numel(category_names);
    
    % 3. Loop over G/H/R/S/T categories.
    
    % internal_avg_accuracies stores the means for G/H/R/S/T.
    internal_avg_accuracies = nan(num_categories, 1);
    
    for k = 1:num_categories
        cat_name = category_names{k};
        FR_k = FR_by_category{k};
        
        if isempty(FR_k)
            fprintf('Category ''%s'' has no data. Skipping.\n', cat_name);
            continue;
        end
        
        [N_neurons, C_conditions, R_trials] = size(FR_k);
        
        fprintf('\n--- Processing INTERNAL OVO for: %s ---\n', cat_name);
        fprintf('  (N=%d, C_conditions=%d, R_trials=%d)\n', ...
            N_neurons, C_conditions, R_trials);
        
        % Check that there are enough data for one-versus-one decoding.
        if C_conditions < 2
            fprintf('  -> Skipping: Not enough conditions (C=%d) for OVO.\n', C_conditions);
            continue;
        end
        
        % Check the trial count R.
        if R_trials < 2
            fprintf('  -> Skipping: Not enough trials (R=%d) for CV.\n', R_trials);
            continue;
        end
        
        % Store all pairwise accuracies within this category.
        all_pair_accs = [];
        
        % 4. Internal one-versus-one loop: Cond_i versus Cond_j.
        for i = 1:(C_conditions - 1)
            for j = (i + 1):C_conditions
                
                % 4.1. Prepare X and Y.
                % [N x 1 x R] -> [N x R]' -> [R x N]
                X_i = squeeze(FR_k(:, i, :))';
                X_j = squeeze(FR_k(:, j, :))';
                
                if N_neurons == 1
                    X_i = X_i';
                    X_j = X_j';
                end
                
                X_total = [X_i; X_j]; % [2*R x N]
                Y_total = [ones(R_trials, 1); ones(R_trials, 1) * 2]; % [2*R x 1]
                
                total_samples_M = size(X_total, 1); % M = 2 * R_trials
                current_K = min(CV_FOLDS, total_samples_M);
                
                % Skip the pair if K < 2, because cross-validation is not possible.
                if current_K < 2
                    fprintf('  -> Skipping pair (%d vs %d): Not enough samples (M=%d) for CV.\n', i, j, total_samples_M);
                    continue; 
                end

                % Version 5 fix: create a nonstratified CV partition manually.
                % This handles very small R, such as R=2, M=4 and K=4,
                % where stratified cross-validation fails or issues warnings.
                
                rng(RANDOM_SEED); % Make the partition reproducible.
                
                % Check for leave-one-out cross-validation, K == M.
                if current_K == total_samples_M
                    c = cvpartition(total_samples_M, 'LeaveOut');
                else
                    % Pass the total sample count M, total_samples_M,
                    % rather than the label vector Y, Y_total,
                    % to create a nonstratified partition.
                    c = cvpartition(total_samples_M, 'KFold', current_K);
                end
                % End of fix.

                
                % 4.2. Train the SVM.
                % Pass CVPartition object c instead of KFold.
                SVM_Model = fitcsvm(X_total, Y_total, ...
                    'KernelFunction', SVM_KERNEL, ...
                    'Standardize', true, ...
                    'CVPartition', c); % Use the manually created partition c.
                
                % 4.3. Compute and store accuracy.
                class_error = kfoldLoss(SVM_Model, 'LossFun', 'ClassifError');
                accuracy = (1 - class_error) * 100;
                
                all_pair_accs(end+1) = accuracy;
                
            end % end OVO loop (j)
        end % end OVO loop (i)
        
        % 5. Compute mean accuracy within this category.
        if isempty(all_pair_accs)
            fprintf('  -> No pairs were decoded.\n');
        else
            avg_acc_k = mean(all_pair_accs);
            internal_avg_accuracies(k) = avg_acc_k;
            fprintf('  -> Avg Internal Accuracy for %s: %.2f%% (over %d pairs)\n', ...
                cat_name, avg_acc_k, numel(all_pair_accs));
        end
        
    end % end category loop (k)

    % 6. Save this file's results.
    save_name = ['Internal_Decoder_Acc_' fileName];
    save_path = fullfile(output_directory, save_name);
    
    fprintf('\nSaving internal results to: %s\n', save_path);
    
    % Store the five means in a readable table.
    internal_accuracy_table = array2table(internal_avg_accuracies, ...
        'VariableNames', {'Average_Internal_Accuracy'}, ...
        'RowNames', category_names);
    
    disp(internal_accuracy_table);
    
    save(save_path, 'internal_accuracy_table', ...
         'internal_avg_accuracies', 'original_meta', '-v7.3');
    
end % End of the file loop.

fprintf('\n\n======================================================\n');
fprintf('All files processed. Script finished.\n');
end
