function output_dir = run_analysis(task, session, max_files)
% Recompute metrics from packaged firing rates; preserve all stored neurons.
% Examples: run_analysis('sdi','s1'); run_analysis('cross_svm','s2',1)
% Outputs are separate from the cached inputs used by run_figures.
if nargin < 1, task = 'sdi'; end
if nargin < 2, session = 's1'; end
if nargin < 3, max_files = Inf; end
session = char(session); task = char(task);
valid_sessions = {'s1','s2','s3','s4','s5','extra_var'};
assert(ismember(session,valid_sessions),'Unknown session.');
assert(isscalar(max_files) && max_files >= 1 && (isinf(max_files) || fix(max_files)==max_files), 'max_files must be a positive integer or Inf.');
valid_tasks = {'sdi','within_svm','cross_svm','internal_svm','geometry', ...
    'category_sdi','global_sdi','global_svm','manifold_var','manifold_global','cvpca'};
assert(ismember(task,valid_tasks),'Unknown task. See help run_analysis.');
root = fileparts(mfilename('fullpath')); addpath(root,fullfile(root,'analysis'));
mode = 'var';
if ismember(task,{'global_sdi','global_svm','manifold_global','cvpca'}), mode = 'global'; end
assert(~(strcmp(session,'extra_var') && ismember(mode,{'global','prototype'})), ...
    'extra_var has variation data only.');
if ismember(task,{'category_sdi','cross_svm'}), mode = 'variation'; end
if ismember(task,{'global_sdi','global_svm'}), mode = 'prototype'; end
input_dir = fullfile(root,'data',session,'fr',mode);
files = dir(fullfile(input_dir,'*.mat'));
assert(~isempty(files),'No packaged firing rates for %s/%s.',session,mode);
[~,order] = sort({files.name}); files = files(order);
files = files(1:min(numel(files),max_files));
% Internal legacy IDs select the original stimulus definitions only.
legacy_ids = [1 4 5 6 7 2]; legacy_id = legacy_ids(strcmp(session,valid_sessions));
output_dir = fullfile(root,'outputs_analysis',session,task);
if ~exist(output_dir,'dir'), mkdir(output_dir); end
stage_root = tempname; mkdir(stage_root);
cleanup = onCleanup(@() rmdir(stage_root,'s'));
stage = fullfile(stage_root,session,mode); mkdir(stage);
for k = 1:numel(files)
    old_name = regexprep(files(k).name,'^FR_(s\d+|extra_var)_',sprintf('FR_s%d_',legacy_id));
    copyfile(fullfile(files(k).folder,files(k).name),fullfile(stage,old_name));
end
switch task
    case 'sdi', calculate_sdi(stage,output_dir);
    case 'within_svm', calculate_within_svm(stage,output_dir);
    case {'cross_svm','global_svm'}, calculate_cross_svm(stage,output_dir);
    case 'internal_svm', calculate_internal_svm(stage,output_dir);
    case 'geometry', calculate_geometry(stage,output_dir);
    case {'category_sdi','global_sdi'}, calculate_global_sdi(stage,output_dir);
    case {'manifold_var','manifold_global'}
        options = struct('kappa',0,'n_t',1000,'flag_NbyM',1,'center_scale',1);
        session_results = struct('session',session,'stimulus_definition',legacy_id,'num_files',numel(files),'files',[]);
        for k = 1:numel(files)
            name = regexprep(files(k).name,'^FR_(s\d+|extra_var)_',sprintf('FR_s%d_',legacy_id));
            output = analyze_manifold_file(fullfile(stage,name),options);
            save(fullfile(output_dir,['Manifold_' files(k).name]),'output','options');
            session_results.files(k).filename = files(k).name;
            session_results.files(k).manifold_details = output.manifold_details;
        end
        save(fullfile(output_dir,'session_results.mat'),'session_results','options');
    case 'cvpca'
        calculate_cvpca_distances(stage_root,struct('output_root',output_dir));
end
% Restore release IDs in generated filenames without altering result values.
renumber_outputs(output_dir,legacy_id,session);
run_info = struct('task',task,'session',session,'input_files',{{files.name}}, ...
    'neuron_filter','none','generated_at',datestr(now,30));
save(fullfile(output_dir,'run_info.mat'),'run_info');
fprintf('Analysis output: %s\n',output_dir);
end

function renumber_outputs(folder,legacy_id,session)
files = dir(folder);
for k = 1:numel(files)
    name = files(k).name;
    if ismember(name,{'.','..'}), continue; end
    p = fullfile(folder,name);
    if files(k).isdir, renumber_outputs(p,legacy_id,session); continue; end
    new_name = regexprep(name,sprintf('(?<=_)s%d(?=_)',legacy_id),session);
    if ~strcmp(new_name,name), movefile(p,fullfile(folder,new_name)); end
end
end
