function path = figure_input(relative)
% Read packaged historical inputs; no screening or external data is used.
p = figure_paths();
if startsWith(relative, 'figure_exports/')
    path = figure_output(extractAfter(relative, 'figure_exports/'));
elseif contains(relative, '/results/old/')
    path = figure_output(extractAfter(relative, '/results/old/'));
else
    path = fullfile(p.data_root, relative);
end
end
