function angles = orthogonality_daily(FR_by_category)
% Preserve the original Fig5B estimator; data and stage are selected upstream.
[P,~,~]=getPara(FR_by_category);
angles=nan(1,3); ids={[3 4 5],1,2}; px={'ang_freq','ori_deg','ori_deg'};py={'radial_freq','freq_cpd','radial_freq'};
for g=1:3
 [X,x,y]=stitch_data(FR_by_category,P,ids{g},px{g},py{g});
 if isempty(X),continue;end
 if g==1,[xn,yn]=normalize_params(x,y,'linear');angles(g)=calculate_cartesian_orthogonality(X,xn,yn);
 else,angles(g)=calculate_cylinder_orthogonality_original(X,x,y);end
end
end
function [X_combined, p_x_raw, p_y_raw] = stitch_data(FR_by_category, P, indices, p_x_name, p_y_name)

    X_combined = [];
    p_x_raw = [];
    p_y_raw = [];

    table_names = {'G', 'H', 'R', 'S', 'T'};

    for c = indices
        if c > length(FR_by_category)
            continue;
        end

        X_c = mean(FR_by_category{c}, 3)';  % [conditions x neurons]
        if isempty(X_c) || size(X_c, 1) < 2
            X_combined = [];
            p_x_raw = [];
            p_y_raw = [];
            return;
        end

        T_name = table_names{c};
        T_params = P.(T_name);

        if ~ismember(p_x_name, T_params.Properties.VariableNames) || ...
           ~ismember(p_y_name, T_params.Properties.VariableNames)
            X_combined = [];
            p_x_raw = [];
            p_y_raw = [];
            return;
        end

        X_combined = [X_combined; X_c];
        p_x_raw = [p_x_raw; T_params.(p_x_name)];
        p_y_raw = [p_y_raw; T_params.(p_y_name)];
    end
end

function [p_x_norm, p_y_norm] = normalize_params(p_x_raw, p_y_raw, p_type) %#ok<INUSD>
    p_x_norm = (p_x_raw - min(p_x_raw)) / (range(p_x_raw) + eps);
    p_y_norm = (p_y_raw - min(p_y_raw)) / (range(p_y_raw) + eps);
end

function ortho_val = calculate_cartesian_orthogonality(X_combined, p_x_norm, p_y_norm)

    X_centered = X_combined - mean(X_combined, 1);
    [~, score, ~, explained] = pca(X_centered);

    cum_var = cumsum(explained);
    k_90 = find(cum_var >= 90, 1);
    k_final = max(k_90, 3);
    k_final = min(k_final, size(score, 2));

    X_pca_z = zscore(score(:, 1:k_final));
    B_x = pinv(X_pca_z) * p_x_norm;
    B_y = pinv(X_pca_z) * p_y_norm;

    if norm(B_x) < eps || norm(B_y) < eps
        ortho_val = NaN;
        return;
    end

    cos_theta = dot(B_x, B_y) / (norm(B_x) * norm(B_y));
    cos_theta = max(min(cos_theta, 1), -1);
    ortho_val = acosd(abs(cos_theta));
end

function ortho_val = calculate_cylinder_orthogonality_original(X_combined, p_x_raw, p_y_raw)

    ori_rad = deg2rad(p_x_raw);
    cos_comp = cos(2 * ori_rad);
    sin_comp = sin(2 * ori_rad);
    p_linear = p_y_raw;

    X_centered = X_combined - mean(X_combined, 1);
    [~, score, ~, explained] = pca(X_centered);

    cum_var = cumsum(explained);
    k_90 = find(cum_var >= 90, 1);
    k_final = max(k_90, 3);
    k_final = min(k_final, size(score, 2));

    X_pca_z = zscore(score(:, 1:k_final));
    B_F = pinv(X_pca_z) * zscore(p_linear);
    B_C = pinv(X_pca_z) * cos_comp;
    B_S = pinv(X_pca_z) * sin_comp;

    if norm(B_F) < eps || norm(B_C) < eps || norm(B_S) < eps
        ortho_val = NaN;
        return;
    end

    cos_FC = dot(B_F, B_C) / (norm(B_F) * norm(B_C));
    cos_FS = dot(B_F, B_S) / (norm(B_F) * norm(B_S));

    max_cos = max(abs(cos_FC), abs(cos_FS));
    max_cos = max(min(max_cos, 1), -1);

    ortho_val = acosd(max_cos);
end