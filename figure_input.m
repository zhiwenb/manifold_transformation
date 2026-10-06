function path = figure_input(relative)
% Resolve the plotting view or output path inside this checkout.
p = figure_paths();
if startsWith(relative, 'figure_exports/')
    path = figure_output(extractAfter(relative, 'figure_exports/'));
elseif contains(relative, '/results/old/')
    path = figure_output(extractAfter(relative, '/results/old/'));
else
    path = fullfile(p.output_root, 'inputs', relative);
end
end
