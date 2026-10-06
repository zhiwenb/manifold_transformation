function [P, scheme, cfg] = distance_parameters(FR_by_category, scheme_in)
% DISTANCE_PARAMETERS  Build parameter mapping tables for {G,H,R,S,T}.
% Usage:
%   load('.../FR_s1_global_041216.mat','FR_by_category');
%   P = distance_parameters(FR_by_category);                  % auto-detect 72/159
%   [P, scheme, cfg] = distance_parameters(FR_by_category);   % also return scheme/cfg
%
% Inputs:
%   FR_by_category : 5x1 cell {G,H,R,S,T}, each [Neuron x Stimuli x Trial]
%   scheme_in      : 'auto' (default) | 'p72' | 'p159'
%
% Outputs:
%   P: struct with tables P.G, P.H, P.R, P.S, P.T, and P.all (merged)
%      Each table has columns:
%        class, stim_idx (1..nStim_in_class), class_local_idx (same as stim_idx),
%        and the parameters you care about:
%        - G:  freq_cpd, ori_deg
%        - H:  radial_freq, ori_deg
%        - R:  ang_freq, radial_freq (dummy = 0; only ang_freq is meaningful)
%        - S:  ang_freq, radial_freq, (hand if applicable)
%        - T:  radial_freq, ang_freq (dummy = 0; only radial_freq is meaningful)
%   scheme: 'p72' or 'p159' (auto-detected if not provided)
%   cfg:    the numeric/string grids actually used to instantiate tables

    if nargin < 2, scheme_in = 'auto'; end
    assert(iscell(FR_by_category) && numel(FR_by_category)==5, ...
        'FR_by_category must be 5x1 cell in order {G,H,R,S,T}.');

    % --------- auto detect 72 vs 159 by per-class counts ----------
    nStim = cellfun(@(X) size(X,2), FR_by_category(:)');
    % G=36,H=12,R=8,S=12,T=4  (p72)
    % G=48,H=24,R=6,S=72,T=6  (p159)
    if strcmpi(scheme_in,'auto')
        if isequal(nStim, [36 12 8 12 4])
            scheme = 'p72';
        elseif isequal(nStim, [48 24 6 72 6])
            scheme = 'p159';
        else
            % fallback heuristic
            if nStim(4) >= 40 || nStim(1) >= 40
                scheme = 'p159';
            else
                scheme = 'p72';
            end
            warning('Unusual per-class counts %s. Heuristically using %s.', ...
                mat2str(nStim), scheme);
        end
    else
        scheme = scheme_in;
    end

    % --------- pick parameter grids according to scheme ----------
    switch lower(scheme)
        case 'p72'
            cfg = cfg_72();
        case 'p159'
            cfg = cfg_159();
        otherwise
            error('Unknown scheme: %s', scheme);
    end

    % check counts match
    expected = [numel(cfg.G.freq)*numel(cfg.G.ori_deg), ...
                numel(cfg.H.radial_freq)*numel(cfg.H.ori_deg), ...
                numel(cfg.R.ang_freq), ...
                numel(cfg.S.ang_freq)*numel(cfg.S.radial_freq)* ...
                    max(1,numel(getfielddefault(cfg.S,'hand',[]))), ...
                numel(cfg.T.radial_freq)];
    assert(all(nStim == expected), ...
        'Per-class stimulus counts %s do not match expected %s for %s.', ...
        mat2str(nStim), mat2str(expected), scheme);

    % --------- build per-class tables ----------
    % G: (order depends on scheme)
    P.G = build_table_block(nStim(1), "G", ...
        struct('freq_cpd', cfg.G.freq, 'ori_deg', cfg.G.ori_deg), cfg.G.order);

    % H: radial first
    P.H = build_table_block(nStim(2), "H", ...
        struct('radial_freq', cfg.H.radial_freq, 'ori_deg', cfg.H.ori_deg), cfg.H.order);

    % R: angular + dummy radial_freq = 0 (only ang_freq is meaningful)
    P.R = table((1:nStim(3))', repmat("R",nStim(3),1), cfg.R.ang_freq(:), ...
        'VariableNames', {'stim_idx','class','ang_freq'});
    P.R.class_local_idx = P.R.stim_idx;
    P.R.radial_freq     = zeros(nStim(3),1);  % dummy radial frequency = 0

    % S: ang, radial, optional hand
    fieldsS = struct('ang_freq', cfg.S.ang_freq, 'radial_freq', cfg.S.radial_freq);
    if isfield(cfg.S,'hand') && ~isempty(cfg.S.hand)
        fieldsS.hand = cfg.S.hand;
    end
    P.S = build_table_block(nStim(4), "S", fieldsS, cfg.S.order);

    % T: radial + dummy ang_freq = 0 (only radial_freq is meaningful)
    P.T = table((1:nStim(5))', repmat("T",nStim(5),1), cfg.T.radial_freq(:), ...
        'VariableNames', {'stim_idx','class','radial_freq'});
    P.T.class_local_idx = P.T.stim_idx;
    P.T.ang_freq        = zeros(nStim(5),1);  % dummy angular frequency = 0

    % --------- merged table in a canonical column order ----------
    tabs = standardizeVarOrder({P.G,P.H,P.R,P.S,P.T});
    P.all = vertcat(tabs{:});
end

% ================= helpers & configs =================

function cfg = cfg_159()
    six_list = [1 2 4 8 12 16];

    % S (72): radial x angular x hand(2)
    % order: hand -> ang_freq -> radial_freq  (hand changes fastest)
    cfg.S.radial_freq = six_list;
    cfg.S.ang_freq    = six_list;
    cfg.S.hand        = ["L","R"];
    cfg.S.order       = 'hand->ang_freq->radial_freq';

    % R (6): angular-only (file order)
    cfg.R.ang_freq    = six_list;

    % T (6): radial-only
    cfg.T.radial_freq = six_list;

    % H (24): radial x ori, radial first
    cfg.H.radial_freq = six_list;
    cfg.H.ori_deg     = [0 45 90 135];
    cfg.H.order       = 'radial_freq->ori_deg';

    % G (48): freq x ori, freq first
    cfg.G.freq        = six_list;
    cfg.G.ori_deg     = (0:7)*(180/8);    % 0:22.5:157.5
    cfg.G.order       = 'freq_cpd->ori_deg';
end

function cfg = cfg_72()
    % S1 (72 prototypes)
    % G: 12 orientations x 3 spatial frequencies (ori first)
    cfg.G.ori_deg     = 0:15:165;           % 12
    cfg.G.freq        = [4 8 16];           % 3
    cfg.G.order       = 'ori_deg->freq_cpd';

    % H: 3 radial x 4 orientations (radial first)
    cfg.H.radial_freq = [4 6 8];            % 3
    cfg.H.ori_deg     = [0 45 90 135];      % 4
    cfg.H.order       = 'radial_freq->ori_deg';

    % R: 8 variants; only angular frequency matters
    cfg.R.ang_freq    = [1 1 1 1 2 4 5 8];

    % S: angular x radial (4 x 3), angular first; no handness
    cfg.S.ang_freq    = [2 4 6 8];
    cfg.S.radial_freq = [2 4 6];
    cfg.S.order       = 'ang_freq->radial_freq';

    % T: 4 radial frequencies
    cfg.T.radial_freq = [4 6 8 10];
end

function T = build_table_block(nStim, classTag, fieldsStruct, orderStr)
% Build a table for one class given the parameter grids and an order string.
% orderStr defines unfolding order, e.g.:
%   'ori_deg->freq_cpd' or 'hand->ang_freq->radial_freq'
% The first dimension in order changes fastest.

    names = fieldnames(fieldsStruct);

    % parse order string
    order = split(strrep(orderStr,' ',''),'->');
    order = string(order(:))';
    assert(all(ismember(order, string(names))), ...
        'order contains undefined field(s).');

    % expected stimulus count
    nexp = 1;
    for k = 1:numel(order)
        nexp = nexp * numel(fieldsStruct.(order{k}));
    end
    assert(nStim==nexp, ...
        'Count mismatch for class %s: nStim=%d vs expected=%d (order=%s).', ...
        classTag, nStim, nexp, orderStr);

    % linear unravel: first term in order changes fastest
    lin = 0:(nStim-1);
    idxMap = struct(); stride = 1;
    for k = 1:numel(order)
        nm = order{k};
        Lk = numel(fieldsStruct.(nm));
        idxMap.(nm) = mod( floor(lin/stride), Lk ) + 1;
        stride = stride * Lk;
    end

    % assemble table: stim_idx, class
    T = table((1:nStim)', repmat(string(classTag), nStim,1), ...
        'VariableNames', {'stim_idx','class'});

    % add parameter columns following original fieldnames order
    for i = 1:numel(names)
        nm  = names{i};
        vec = fieldsStruct.(nm);
        ix  = idxMap.(string(nm));
        if isstring(vec) || ischar(vec)
            T.(nm) = vec(ix)';
        else
            T.(nm) = vec(ix)';
        end
    end
    T.class_local_idx = T.stim_idx;
end

function tabs = standardizeVarOrder(tablist)
% Standardize column order so that vertcat across classes is consistent.

    pref = {'class','stim_idx','class_local_idx', ...
            'freq_cpd','ori_deg','radial_freq','ang_freq','hand'};

    % If any table has 'hand', we will add missing 'hand' as string in others.
    hasHand = any(cellfun(@(T) any(strcmp('hand', T.Properties.VariableNames)), tablist));

    tabs = cell(size(tablist));
    for i = 1:numel(tablist)
        Ti = tablist{i};

        % ensure preferred columns exist
        for nm = pref
            nm = nm{1};
            if ~ismember(nm, Ti.Properties.VariableNames)
                if strcmp(nm,'hand') && hasHand
                    Ti.(nm) = strings(height(Ti),1);
                else
                    Ti.(nm) = NaN(height(Ti),1);
                end
            end
        end

        % preferred columns first, then the remaining columns in original order
        keep = pref(ismember(pref, Ti.Properties.VariableNames));
        rest = setdiff(Ti.Properties.VariableNames, keep, 'stable');
        Ti   = Ti(:, [keep, rest]);
        tabs{i} = Ti;
    end
end

function val = getfielddefault(S, field, defaultVal)
% Simple helper: get S.field if exists, otherwise defaultVal.
    if isfield(S, field)
        val = S.(field);
    else
        val = defaultVal;
    end
end
