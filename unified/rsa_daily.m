function values = rsa_daily(FR_by_category)
% Use the shared getPara mapping for all representation metrics.
[P,~,~]=getPara(FR_by_category);
values=nan(1,3);ids={[3 4 5],1,2};px={'ang_freq','ori_deg','ori_deg'};py={'radial_freq','freq_cpd','radial_freq'};
for g=1:3
 [X,x,y]=stitch_data(FR_by_category,P,ids{g},px{g},py{g});
 if ~isempty(X)&&size(X,1)>1,values(g)=calculate_rsa_mixed_model(X,x,y,px{g},py{g});end
end
end
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

