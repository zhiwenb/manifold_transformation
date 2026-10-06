function calculate_cvpca_distances(base_dir, opts)
% generate_distance_matrices_batch
%
% Step 1 (cvPCA version):
%   For each session s*/global/:
%     1. Load FR_*.mat (FR_by_category: {G,H,R,S,T}, each [N x S x R])
%     2. For each class (G,H,R,S,T):
%          - Denoise with cvPCA in stimulus dimension
%          - Compute neural distance matrix (Neural)
%          - Compute parameter distance matrix (Parameter)
%     3. Build a merged Polar manifold:
%          - Concatenate R, S, T in stimulus dimension
%          - Compute DistMatrices.Polar.Neural and DistMatrices.Polar.Parameter
%
% Inputs:
%   base_dir : root folder containing s*/global/FR_*.mat
%              default: .../Data_FR
%   opts     : struct with fields
%       .metric   : 'euclidean' (default) or 'mahal'
%       .cv       : struct, cvPCA options
%           .split_mode : 'odd-even' (default) or 'random'
%           .seed       : rng seed for random split
%           .keep_rule  : 'positive' (default) | 'topK' | 'var90'
%           .topK       : K for 'topK'
%
% Output:
%   Saves DistMatrices & CVInfo into
%     ../Results_Distance_Matrices_cvPCA_YYYYMMDD/sX/Distances_FR_*.mat
%
% NOTE on parameter distance:
%   All frequency-like parameters are transformed with
%       f(v) = log2(1 + max(v,0))
%   before computing absolute differences in parameter space.

    % ---------- defaults ----------
    if nargin < 1 || isempty(base_dir)
        error('Provide the firing-rate input directory.');
    end
    if nargin < 2, opts = struct(); end
    opts = fill_default_opts(opts);

    assert(isfolder(base_dir), 'Base dir not found: %s', base_dir);

    out_root = opts.output_root;
    if ~exist(out_root,'dir'), mkdir(out_root); end

    sess_dirs = list_subdirs(base_dir, '^s\d+$');

    % ---------- loop over sessions ----------
    for si = 1:numel(sess_dirs)
        sess = sess_dirs{si};
        gdir = fullfile(base_dir, sess, 'global');
        files = list_files(gdir, '^FR_.*\.mat$');
        fprintf('== Session %s: %d files ==\n', sess, numel(files));

        if isempty(files), continue; end

        sess_out = fullfile(out_root, sess);
        if ~exist(sess_out,'dir'), mkdir(sess_out); end

        % ------ loop over FR files within session ------
        for fi = 1:numel(files)
            fpath = fullfile(gdir, files{fi});
            [~, fname] = fileparts(files{fi});
            out_fpath = fullfile(sess_out, ['Distances_' fname '.mat']);
            fprintf('Processing %s ...\n', files{fi});

            % --- load FR_by_category {G,H,R,S,T} ---
            S = load(fpath, 'FR_by_category');
            assert(isfield(S,'FR_by_category'), ...
                'Missing FR_by_category in %s', fpath);
            FRc = S.FR_by_category;   % {G,H,R,S,T}, each [N x S x R]

            % build parameter tables (uses your getPara.m)
            [P, ~, ~] = distance_parameters(FRc);

            cats = {'G','H','R','S','T'};
            DistMatrices = struct();  % output container
            CVInfo       = struct();  % cvPCA metadata

            % ---------- per-class manifolds (G/H/R/S/T) ----------
            for k = 1:5
                Cname   = cats{k};
                X_all   = FRc{k};      % [N x S x R]
                Tab_all = P.(Cname);   % table(S x p)

                if isempty(X_all) || isempty(Tab_all) || size(X_all,2) < 2
                    continue;
                end

                % optional hand grouping (for 159 case: S has hand L/R)
                use_hand = (ismember('hand', Tab_all.Properties.VariableNames) && ...
                            numel(unique(Tab_all.hand)) > 1);

                if use_hand
                    fprintf('    Found ''hand'' groups in %s. Processing L/R separately.\n', Cname);
                    groups = unique(string(Tab_all.hand));

                    for gi = 1:numel(groups)
                        g_name = strtrim(groups(gi));
                        if g_name == "", continue; end

                        idx_g  = (string(Tab_all.hand) == g_name);
                        Tab_g  = Tab_all(idx_g, :);
                        X_g    = X_all(:, idx_g, :);  % [N x Sg x R]

                        if size(X_g,2) < 2, continue; end

                        % cvPCA denoising: [N x Sg]
                        [X_denoised, meta] = cvpca_denoise_block(X_g, opts.cv);

                        % neural distances on denoised trial-avg
                        v_neural = neural_pdist(X_denoised, opts.metric);
                        DistMatrices.(Cname).(g_name).Neural = squareform(v_neural);

                        % parameter distances (with log2(1+freq) transform)
                        DistMatrices.(Cname).(g_name).Parameter = ...
                            compute_parameter_distance_matrix(Tab_g);

                        CVInfo.(Cname).(g_name) = meta;
                    end

                else
                    % no hand -> whole class as one block
                    [X_denoised, meta] = cvpca_denoise_block(X_all, opts.cv);

                    v_neural = neural_pdist(X_denoised, opts.metric);
                    DistMatrices.(Cname).Neural = squareform(v_neural);

                    DistMatrices.(Cname).Parameter = ...
                        compute_parameter_distance_matrix(Tab_all);

                    CVInfo.(Cname) = meta;
                end
            end  % end per-class loop

            % ---------- Polar manifold (R+S+T) ----------
            % FR_by_category: {G,H,R,S,T} -> 3=R, 4=S, 5=T
            X_R = FRc{3};
            X_S = FRc{4};
            X_T = FRc{5};

            if ~isempty(X_R) && ~isempty(X_S) && ~isempty(X_T) ...
                    && size(X_R,1)==size(X_S,1) && size(X_R,1)==size(X_T,1) ...
                    && size(X_R,3)==size(X_S,3) && size(X_R,3)==size(X_T,3)

                % concatenate stimuli [N x (S_R+S_S+S_T) x R]
                X_Polar = cat(2, X_R, X_S, X_T);

                % build parameter table with only radial_freq & ang_freq
                % NOTE: R/T should have dummy freq = 0 set in getPara.m
                Tab_R = unify_param_table_for_polar(P.R);
                Tab_S = unify_param_table_for_polar(P.S);
                Tab_T = unify_param_table_for_polar(P.T);

                if height(Tab_R) ~= size(X_R,2) || ...
                   height(Tab_S) ~= size(X_S,2) || ...
                   height(Tab_T) ~= size(X_T,2)

                    warning('  Polar: parameter table size mismatch with FR. Skipping Polar.');
                else
                    Tab_Polar = [Tab_R; Tab_S; Tab_T];

                    % cvPCA on combined R+S+T
                    [X_Polar_denoised, meta_Polar] = cvpca_denoise_block(X_Polar, opts.cv);

                    % neural distances
                    v_neural_Polar = neural_pdist(X_Polar_denoised, opts.metric);
                    DistMatrices.Polar.Neural = squareform(v_neural_Polar);

                    % parameter distances (two frequency dims only)
                    DistMatrices.Polar.Parameter = ...
                        compute_parameter_distance_matrix(Tab_Polar);

                    CVInfo.Polar = meta_Polar;

                    fprintf('  Built Polar manifold (R+S+T, %d stimuli).\n', size(X_Polar,2));
                end
            else
                fprintf('  Skip Polar: R/S/T empty or size mismatch.\n');
            end

            % ---------- save ----------
            if isempty(fieldnames(DistMatrices))
                warning('  No valid manifold data found in %s. Skipping save.', fname);
                continue;
            end

            save(out_fpath, 'DistMatrices', 'CVInfo', '-v7.3');
            fprintf('  Saved: %s\n', out_fpath);
        end
    end

    fprintf('\nALL DONE. Distance matrices saved to:\n  %s\n', out_root);
end

%% ===================== cvPCA denoising core =====================
function [M_denoised, meta] = cvpca_denoise_block(X, cvopt)
% X: [N x S x R] (neurons x stimuli x repeats)
% Returns:
%   M_denoised : [N x S] trial-averaged responses projected into
%                the "reliable signal" subspace
%   meta       : struct with K, cv_eigs, frac, split indices, etc.

    [N,S,R] = size(X);

    % no repeats: fall back to plain trial average
    if R < 1
        warning('R<1: no repeats; using plain trial-avg.');
        M_denoised = mean(X, 3, 'omitnan'); % [N x S]
        meta = struct('K', 0, 'cv_eigs', [], 'frac', [], 'split', []);
        return;
    elseif R == 1
        % cannot cross-validate, do simple PCA denoising (no component cut)
        M  = mean(X,3,'omitnan');   % [N x S]
        M0 = center_neuronwise(M);
        [U,~,~] = svd(M0,'econ');
        K = min(S,N);
        U_keep = U(:,1:K);
        M_denoised = U_keep*(U_keep'*M0);
        meta = struct('K',K,'cv_eigs',[],'frac',[],'split','single');
        return;
    end

    % ---- split repeats into A / B ----
    idx = 1:R;
    if strcmpi(cvopt.split_mode,'random')
        rng(cvopt.seed);
        idx = idx(randperm(R));
    end
    A = idx(1:floor(R/2));
    B = setdiff(idx, A);
    if isempty(B)
        B = idx(ceil(R/2)+1:end);
    end

    % stimulus-averaged responses for each half + neuron-wise centering
    RA  = mean(X(:,:,A), 3, 'omitnan');   % [N x S]
    RB  = mean(X(:,:,B), 3, 'omitnan');
    RA0 = center_neuronwise(RA);
    RB0 = center_neuronwise(RB);

    % SVD on half A
    [U,~,~] = svd(RA0, 'econ');          % U: [N x Kmax]

    % project A/B onto U
    TA = U' * RA0;                       % [Kmax x S]
    TB = U' * RB0;

    % cross-validated eigenvalues (reliable variance)
    lam_cv = sum(TA .* TB, 2) ./ max(S-1,1);   % [Kmax x 1]

    % choose components
    [keep_idx, K, frac] = choose_components(lam_cv, cvopt);
    if K == 0
        % no positive reliable variance -> plain average (no denoising)
        M_denoised = mean(X,3,'omitnan');
        meta = struct('K',0,'cv_eigs',lam_cv,'frac',frac, ...
                      'split',struct('A',A,'B',B));
        return;
    end
    U_keep = U(:, keep_idx);

    % project all repeats into signal subspace and reconstruct
    Mall  = mean(X, 3, 'omitnan');   % [N x S]
    Mall0 = center_neuronwise(Mall);
    M_denoised = U_keep * (U_keep' * Mall0);

    meta = struct('K',K,'cv_eigs',lam_cv,'frac',frac, ...
                  'split',struct('A',A,'B',B));
end

function X0 = center_neuronwise(X)
% Subtract mean over stimuli for each neuron (no variance scaling).
    mu = mean(X,2,'omitnan');   % [N x 1]
    X0 = X - mu;
end

function [keep_idx, K, frac] = choose_components(lam_cv, cvopt)
% Choose which components to keep based on cv eigenvalues.
% lam_cv: [Kmax x 1] reliable variance spectrum.

    lam = lam_cv(:);
    switch lower(cvopt.keep_rule)
        case 'positive'
            keep_idx = find(lam > 0);
        case 'topk'
            k = min(cvopt.topK, numel(lam));
            [~,ord] = sort(lam,'descend');
            keep_idx = ord(1:k);
        case 'var90'
            [v,ord] = sort(max(lam,0),'descend');
            csum = cumsum(v);
            thr  = 0.90 * sum(v);
            k = find(csum>=thr,1,'first');
            if isempty(k), k = 0; end
            keep_idx = ord(1:k);
        otherwise
            keep_idx = find(lam > 0);
    end

    keep_idx = keep_idx(:).';
    K = numel(keep_idx);

    pos = max(lam,0);
    if sum(pos)>0
        frac = cumsum(sort(pos,'descend')) / sum(pos);
    else
        frac = [];
    end
end

%% ===================== Polar parameter table helper =====================
function T2 = unify_param_table_for_polar(T)
% For the Polar manifold we only keep:
%   class, radial_freq, ang_freq
% (Two "frequency" dimensions; handness is ignored here.)
    allvars = {'class','radial_freq','ang_freq'};

    % ensure columns exist
    for i = 1:numel(allvars)
        nm = allvars{i};
        if ~ismember(nm, T.Properties.VariableNames)
            if strcmp(nm, 'class')
                T.(nm) = strings(height(T),1);
            else
                T.(nm) = NaN(height(T),1);
            end
        end
    end

    % fixed column order
    T2 = T(:, allvars);
end

%% ===================== parameter distance matrix =====================
function D_out = compute_parameter_distance_matrix(Tclass)
% Compute parameter-space distance matrices based on any subset of:
%   freq_cpd, radial_freq, ang_freq, ori_deg
%
% NOTE:
%   For all frequency-like parameters we use
%       f(v) = log2(1 + max(v,0))
%   so that dummy = 0 and large frequencies are compressed (more perceptual).

    D_out = struct();
    params = intersect({'freq_cpd','radial_freq','ang_freq','ori_deg'}, ...
                       Tclass.Properties.VariableNames, 'stable');
    S = height(Tclass);
    if isempty(params) || S < 2
        D_out.Composite = nan(S, S);
        return;
    end

    use_ori_180 = ismember('ori_deg', params);
    if ismember('class', Tclass.Properties.VariableNames)
        cls = unique(string(Tclass.class));
        use_ori_180 = any(ismember(upper(cls), ["G","H"])) && ...
                      ismember('ori_deg', params);
    end

    comps_vec = {};
    Tsub = Tclass;

    % grating spatial frequency (cpd)
    if ismember('freq_cpd', params)
        v = double(Tsub.freq_cpd(:));
        if ~isempty(v)
            v = max(v,0);
            L = log2(1 + v);      % log2(1 + freq)
            D = abs(L - L');
            D_out.freq_cpd = D;
            comps_vec{end+1} = zvec(squareform(D)); %#ok<AGROW>
        end
    end

    % radial frequency
    if ismember('radial_freq', params)
        v = double(Tsub.radial_freq(:));
        if ~isempty(v)
            v = max(v,0);
            L = log2(1 + v);      % log2(1 + radial)
            D = abs(L - L');
            D_out.radial_freq = D;
            comps_vec{end+1} = zvec(squareform(D)); %#ok<AGROW>
        end
    end

    % angular frequency
    if ismember('ang_freq', params)
        v = double(Tsub.ang_freq(:));
        if ~isempty(v)
            v = max(v,0);
            L = log2(1 + v);      % log2(1 + angular)
            D = abs(L - L');
            D_out.ang_freq = D;
            comps_vec{end+1} = zvec(squareform(D)); %#ok<AGROW>
        end
    end

    % orientation (circular distance)
    if ismember('ori_deg', params)
        th = double(Tsub.ori_deg(:));
        period = 360;
        if use_ori_180, period = 180; end
        Dth = abs(th - th');
        Dth = min(Dth, period - Dth);
        D_out.ori_deg = Dth;
        comps_vec{end+1} = zvec(squareform(Dth)); %#ok<AGROW>
    end

    if isempty(comps_vec)
        D_out.Composite = nan(S, S);
        return;
    end

    v_param = zeros(size(comps_vec{1}));
    for c = 1:numel(comps_vec)
        if all(isnan(comps_vec{c})), continue; end
        v_param = v_param + comps_vec{c};
    end
    D_out.Composite = squareform(v_param);
end

%% ===================== distances & small helpers =====================
function d = neural_pdist(X, metric)
% X: [N x S] (columns = stimuli)
% returns vectorized upper-triangular distances

    switch lower(metric)
        case 'euclidean'
            d = pdist(X', 'euclidean');
        case 'mahal'
            C   = cov(X');
            lam = 0.1;
            Csh = (1-lam)*C + lam*mean(diag(C))*eye(size(C));
            [U,Sv] = eig((Csh+Csh')/2);
            Sv = max(Sv, 1e-6);
            W  = U*diag(1./sqrt(diag(Sv)))*U';
            Xw = W*X;
            d  = pdist(Xw', 'euclidean');
        otherwise
            error('Unknown metric %s', metric);
    end
end

function v = zvec(x)
% z-score a vector, robust to zero std
    mu = mean(x,'omitnan');
    sd = std(x,0,'omitnan');
    if sd==0 || isnan(sd)
        v = x - mu;
    else
        v = (x - mu) ./ sd;
    end
end

function names = list_subdirs(root, pattern)
% list subdirectories whose names match a regexp
    d = dir(root);
    d = d([d.isdir]);
    names = {};
    for i = 1:numel(d)
        n = d(i).name;
        if any(strcmp(n,{'.','..'})), continue; end
        if ~isempty(regexp(n, pattern, 'once'))
            names{end+1} = n; %#ok<AGROW>
        end
    end
end

function files = list_files(root, pattern)
% list files whose names match a regexp
    d = dir(root);
    files = {};
    for i = 1:numel(d)
        if d(i).isdir, continue; end
        if ~isempty(regexp(d(i).name, pattern, 'once'))
            files{end+1} = d(i).name; %#ok<AGROW>
        end
    end
    files = sort(files);
end

function opts = fill_default_opts(opts)
% fill missing fields of opts with defaults
    if ~isfield(opts,'metric')
        opts.metric = 'euclidean';
    end
    if ~isfield(opts,'cv')
        opts.cv = struct();
    end
    if ~isfield(opts.cv,'split_mode')
        opts.cv.split_mode = 'odd-even';   % or 'random'
    end
    if ~isfield(opts.cv,'seed')
        opts.cv.seed = 7;
    end
    if ~isfield(opts.cv,'keep_rule')
        opts.cv.keep_rule = 'positive';    % 'positive' | 'topK' | 'var90'
    end
    if ~isfield(opts.cv,'topK')
        opts.cv.topK = 32;
    end
end
