function p = figure_paths()
% All inputs and outputs are resolved relative to this checkout.
p.project_root = fileparts(mfilename('fullpath'));
p.data_root = fullfile(p.project_root, 'data');
p.output_root = fullfile(p.project_root, 'outputs_simple');
end
