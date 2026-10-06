function [output] = analyze_manifold_file(data_file_path, options)
% Load and analyze a single file

data_in = load(data_file_path);
FR_by_category = data_in.FR_by_category;

% Reshape
P_raw = numel(FR_by_category);
XtotT = cell(P_raw, 1);

for ii = 1:P_raw
    data_3D = FR_by_category{ii};
    if isempty(data_3D)
        continue;
    end
    [N, C_i, R] = size(data_3D);
    M_i = C_i * R;
    XtotT{ii} = reshape(data_3D, [N, M_i]);
end

% Clean
XtotT_clean = XtotT(~cellfun('isempty', XtotT));
P_final = numel(XtotT_clean);
fprintf('  Loaded %d non-empty categories.\n', P_final);

if P_final == 0
    error('No analyzable data in this file.');
end

% Run MFT analysis
output = manifold_stable_analysis_corr(XtotT_clean, options);
end
