function calculate_geometry(input_directory, output_directory)
% Compute saved metrics from firing rates without neuron reselection.
%% =====================================================================
%  BATCH Geometry Calculation: Radius (Local) & Separation (Global)
%
%  Calculates:
%  1. Within-Category Radius (Expansion metric): 
%     - Mean Euclidean distance from variations to category centroid.
%     - Mean Pairwise distance between all variations.
%  2. Between-Category Distance (Compression metric):
%     - Euclidean distance between global centroids of different categories.
% =====================================================================

clc;

%% --- 1. Parameters & Setup ---
% Input data directory containing four prototypes with 45 variations each.



if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

fprintf('Searching for .mat files in: %s\n', input_directory);
mat_files = dir(fullfile(input_directory, '*.mat'));

%% --- File Loop ---
for f_idx = 1:numel(mat_files)
    
    fileName = mat_files(f_idx).name;
    filePath = fullfile(input_directory, fileName);
    
    fprintf('\nProcessing (%d/%d): %s\n', f_idx, numel(mat_files), fileName);
    
    % Load data.
    try
        S = load(filePath, 'FR_by_category', 'meta');
    catch
        continue;
    end
    
    FR_by_category = S.FR_by_category;
    category_names = S.meta.category_names;
    num_cats = numel(category_names);
    
    % Initialize the result structure.
    Geo_Res = struct();
    
    % Store the global centroid of each category for between-category distances.
    All_Centroids = cell(num_cats, 1); 
    
    %% === Part 1: Within-Category Metrics (Local Expansion) ===
    fprintf('  [Level 2] Calculating Local Radius/Expansion...\n');
    
    for c = 1:num_cats
        cat_name = category_names{c};
        FR_cat = FR_by_category{c}; % [Neurons x Conditions x Reps]
        
        if isempty(FR_cat), continue; end
        
        % 1. Reduce to condition means by averaging over trials.
        % Result: [Neurons x Conditions]
        FR_conds_mean = mean(FR_cat, 3, 'omitnan');
        
        % Handle NaNs when a condition has no trials.
        valid_cols = ~any(isnan(FR_conds_mean), 1);
        FR_conds_mean = FR_conds_mean(:, valid_cols);
        
        [N_neurons, N_conds] = size(FR_conds_mean);
        
        % 2. Compute the category centroid across all variations.
        % Result: [Neurons x 1]
        Cat_Centroid = mean(FR_conds_mean, 2);
        All_Centroids{c} = Cat_Centroid; % Store for use in Part 2.
        
        % 3. Metric A: Radial Radius (Distance from Centroid)
        % Compute each variation's distance to the centroid.
        dists_to_center = zeros(N_conds, 1);
        for i = 1:N_conds
            dists_to_center(i) = norm(FR_conds_mean(:, i) - Cat_Centroid);
        end
        mean_radius = mean(dists_to_center);
        
        % 4. Metric B: Pairwise Spread (Average distance between all variations)
        % A more robust dispersion measure.
        p_dists = pdist(FR_conds_mean', 'euclidean'); % Transpose to observations x variables.
        mean_pairwise_dist = mean(p_dists);
        
        % 5. Metric C: Total Variance (Trace of Covariance)
        % Analogous to radius in MFT.
        cov_mat = cov(FR_conds_mean');
        total_var = trace(cov_mat);
        
        % --- Store Results ---
        Geo_Res.Within.(cat_name).Mean_Radius_ToCentroid = mean_radius;
        Geo_Res.Within.(cat_name).Mean_Pairwise_Dist = mean_pairwise_dist;
        Geo_Res.Within.(cat_name).Total_Variance = total_var;
        Geo_Res.Within.(cat_name).N_variations = N_conds;
        
        fprintf('    -> %s: Radius=%.2f, PairDist=%.2f\n', ...
            cat_name, mean_radius, mean_pairwise_dist);
    end
    
    %% === Part 2: Between-Category Metrics (Global Compression) ===
    fprintf('  [Level 3] Calculating Global Separation/Compression...\n');
    
    pairs = nchoosek(1:num_cats, 2);
    
    for p = 1:size(pairs, 1)
        idx1 = pairs(p, 1);
        idx2 = pairs(p, 2);
        
        name1 = category_names{idx1};
        name2 = category_names{idx2};
        
        C1 = All_Centroids{idx1};
        C2 = All_Centroids{idx2};
        
        if isempty(C1) || isempty(C2)
            dist_val = NaN;
        else
            % Compute the distance between the centers of two category clouds.
            dist_val = norm(C1 - C2);
        end
        
        % Construct a field name such as Hyperbolic_vs_Spiral.
        pair_name = sprintf('%s_vs_%s', name1, name2);
        
        Geo_Res.Between.(pair_name) = dist_val;
        
        fprintf('    -> %s: Dist = %.2f\n', pair_name, dist_val);
    end
    
    %% --- Save ---
    save_name = ['Geometry_Results_' fileName];
    save_path = fullfile(output_directory, save_name);
    save(save_path, 'Geo_Res', 'category_names');
    
end

fprintf('\n=== Done. Check results in %s ===\n', output_directory);
end
