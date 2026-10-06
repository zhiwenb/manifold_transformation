function calculate_cross_svm(input_directory, output_directory)
% Compute saved metrics from firing rates without neuron reselection.
%% ---------------------------------------------------------------------
% BATCH Pairwise SVM Decoding Script (V3.1: With SAVING, safer version)
%
% Objective:
% 1. Find all .mat files in an input folder.
% 2. Loop through each file.
% 3. Run pairwise (OVO) decoding with undersampling.
% 4. Print the accuracy matrix.
% 5. Save the resulting accuracy_table and matrix to a new .mat file.
% ---------------------------------------------------------------------

clc;

%% --- 1. Parameters & Setup ---



if ~exist(output_directory, 'dir')
    fprintf('Creating output directory: %s\n', output_directory);
    mkdir(output_directory);
end

% SVM & CV Parameters
CV_FOLDS = 5;
SVM_KERNEL = 'linear';
RANDOM_SEED = 42;

% Set random seed ONCE
rng(RANDOM_SEED);

%% --- 0. Find all .mat files ---
fprintf('Searching for .mat files in: %s\n', input_directory);
mat_files = dir(fullfile(input_directory, '*.mat'));

if isempty(mat_files)
    error('No .mat files found in the specified directory: %s', input_directory);
end

fprintf('Found %d .mat files. Results will be saved to: %s\n', ...
    numel(mat_files), output_directory);

%% --- START FILE LOOP ---
for f_idx = 1:numel(mat_files)

    fileName = mat_files(f_idx).name;
    filePath = fullfile(input_directory, fileName);

    fprintf('\n\n======================================================\n');
    fprintf('Processing file (%d/%d): %s\n', f_idx, numel(mat_files), fileName);
    fprintf('======================================================\n');

    %% --- 2. Load data ---
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

    if ~isfield(original_meta, 'category_names')
        fprintf('  -> ERROR: meta is missing ''category_names''. Skipping.\n');
        continue;
    end

    category_names = original_meta.category_names;
    num_categories = numel(category_names);

    %% --- 3. Pre-process Data: Reshape ---
    X_all = cell(num_categories, 1);
    Y_all = cell(num_categories, 1);
    sample_counts = zeros(num_categories, 1);
    N_neurons = 0;
    is_data_valid = true;

    for i = 1:num_categories
        FR_k = FR_by_category{i};

        if isempty(FR_k)
            fprintf('  -> WARNING: Category %s is empty.\n', category_names{i});
            continue;
        end

        [N, C, R] = size(FR_k);

        if N_neurons == 0
            N_neurons = N;
        end

        if N ~= N_neurons || N == 0
            is_data_valid = false;
            break;
        end

        num_samples_k = C * R;
        sample_counts(i) = num_samples_k;

        X_all{i} = reshape(FR_k, N, num_samples_k)';
        Y_all{i} = ones(num_samples_k, 1) * i;
    end

    if ~is_data_valid
        fprintf('  -> ERROR: Invalid data (neuron mismatch or 0 neurons). Skipping file.\n');
        continue;
    end

    fprintf('Data reshaped. Found %d neurons.\n', N_neurons);
    fprintf('Samples per category:\n');
    for i = 1:num_categories
        fprintf('  %s: %d samples\n', category_names{i}, sample_counts(i));
    end

    %% --- 4. Pairwise Decoding Loop ---
    fprintf('\n--- Starting Pairwise (One-vs-One) Decoding ---\n');

    accuracy_matrix = nan(num_categories, num_categories);

    for i = 1:num_categories
        for j = (i + 1):num_categories

            X_i = X_all{i}; Y_i = Y_all{i};
            X_j = X_all{j}; Y_j = Y_all{j};

            if isempty(X_i) || isempty(X_j)
                fprintf('  -> Skipping %s vs %s (empty category).\n', ...
                    category_names{i}, category_names{j});
                continue;
            end

            num_samples_i = size(X_i, 1);
            num_samples_j = size(X_j, 1);
            min_samples = min(num_samples_i, num_samples_j);

            if min_samples < CV_FOLDS
                fprintf('  -> Skipping %s vs %s (min_samples=%d < CV_FOLDS=%d).\n', ...
                    category_names{i}, category_names{j}, min_samples, CV_FOLDS);
                continue;
            end

            idx_i = randperm(num_samples_i, min_samples);
            idx_j = randperm(num_samples_j, min_samples);

            X_balanced = [X_i(idx_i, :); X_j(idx_j, :)];
            Y_balanced = [Y_i(idx_i, :); Y_j(idx_j, :)];

            try
                SVM_Model = fitcsvm(X_balanced, Y_balanced, ...
                    'KernelFunction', SVM_KERNEL, ...
                    'Standardize', true, ...
                    'CrossVal', 'on', ...
                    'KFold', CV_FOLDS);

                class_error = kfoldLoss(SVM_Model, 'LossFun', 'ClassifError');
                accuracy = (1 - class_error) * 100;

                accuracy_matrix(i, j) = accuracy;
                accuracy_matrix(j, i) = accuracy;

                fprintf('  -> %s vs %s: %.2f%%\n', ...
                    category_names{i}, category_names{j}, accuracy);

            catch e
                fprintf('  -> ERROR in %s vs %s: %s\n', ...
                    category_names{i}, category_names{j}, e.message);
            end
        end
    end

    fprintf('\n--- Results for %s ---\n', fileName);

    %% --- 5. Display Results ---
    valid_var_names = matlab.lang.makeValidName(category_names, 'ReplacementStyle', 'delete');
    valid_var_names = matlab.lang.makeUniqueStrings(valid_var_names);

    accuracy_table = array2table(accuracy_matrix, ...
        'VariableNames', valid_var_names, ...
        'RowNames', category_names);

    disp('Pairwise Decoding Accuracy Matrix (%):');
    disp(accuracy_table);

    %% --- 6. Save Results ---
    save_name = ['Decoder_Acc_' fileName];
    save_path = fullfile(output_directory, save_name);

    decoding_params.CV_FOLDS = CV_FOLDS;
    decoding_params.SVM_KERNEL = SVM_KERNEL;
    decoding_params.RANDOM_SEED = RANDOM_SEED;
    decoding_params.input_directory = input_directory;
    decoding_params.input_file = filePath;
    decoding_params.analysis_type = 'between_category_pairwise_svm_trial_level';

    fprintf('Saving results to: %s\n', save_path);

    save(save_path, 'accuracy_table', 'accuracy_matrix', ...
        'original_meta', 'decoding_params', '-v7.3');
end

fprintf('\n\n======================================================\n');
fprintf('All files processed. Script finished.\n');
end
