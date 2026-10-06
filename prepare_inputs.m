function root = prepare_inputs()
% Materialize lightweight legacy plotting views from the simple session layout.
% Data is stored once in data/sN; these disposable copies stay in outputs.
p = figure_paths();
root = fullfile(p.output_root, 'inputs');
index = jsondecode(fileread(fullfile(p.project_root,'data','input_index.json')));
for k = 1:numel(index)
    source = fullfile(p.project_root, index(k).file);
    target = fullfile(root, index(k).view);
    assert(isfile(source), 'Missing packaged input: %s', source);
    if ~isfile(target)
        parent = fileparts(target);
        if ~exist(parent,'dir'), mkdir(parent); end
        copyfile(source, target);
    end
end
end
