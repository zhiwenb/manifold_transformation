function path = figure_output(relative)
p = figure_paths();
path = fullfile(p.output_root, relative);
[parent,~,ext] = fileparts(path);
if isempty(ext), parent = path; end
if ~exist(parent, 'dir'), mkdir(parent); end
end
