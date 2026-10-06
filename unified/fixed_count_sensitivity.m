function fixed_count_sensitivity()
% Random count matching is a sensitivity analysis, never a primary neuron filter.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'unified'));
fs=dir(fullfile(root,'unified','inputs','union_common240','s*','**','*.mat'));n=Inf;
for k=1:numel(fs),S=load(fullfile(fs(k).folder,fs(k).name),'meta');n=min(n,S.meta.num_neurons);end
seeds=[11 29 47];dest=fullfile(root,'unified','daily','fixed_count_common240');mkdir(dest);drawroot=fullfile(root,'unified','fixed_count_draws');mkdir(drawroot);
for k=1:numel(fs)
 parts=strsplit(fs(k).folder,filesep);session=parts{end-1};mode=parts{end};name=[session '_' mode '_' fs(k).name];
 if isfile(fullfile(dest,name)),continue;end
 S=load(fullfile(fs(k).folder,fs(k).name));draws=cell(1,3);indices=cell(1,3);tic;
 for di=1:3
  rng(seeds(di),'twister');indices{di}=randperm(S.meta.num_neurons,n);meta=S.meta;meta.num_neurons=n;meta.neurons_used=meta.neurons_used(:,indices{di});FR_by_category=S.FR_by_category;
  for g=1:numel(FR_by_category),if ~isempty(FR_by_category{g}),FR_by_category{g}=FR_by_category{g}(indices{di},:,:);end;end
  tmp=[tempname '.mat'];save(tmp,'FR_by_category','meta');cleanup=onCleanup(@()delete(tmp));draws{di}=compute_daily_metrics(tmp,session,mode);draws{di}.file=fullfile(fs(k).folder,fs(k).name);draws{di}.file=draws{di}.file(numel(root)+2:end);draws{di}.date=regexp(fs(k).name,'\d{6}(?=\.mat)','match','once');clear cleanup;
 end
 result=draws{1};fields=fieldnames(result);
 for f=1:numel(fields),field=fields{f};v=result.(field);if isnumeric(v)&&~isempty(v),result.(field)=mean(cat(ndims(v)+1,draws{1}.(field),draws{2}.(field),draws{3}.(field)),ndims(v)+1);end;end
 result.candidate='fixed_count_common240';result.count_match_seeds=seeds;save(fullfile(dest,name),'result');save(fullfile(drawroot,name),'draws','indices','seeds','n');fprintf('Count match %s %s %d units %.1f sec\n',session,result.date,n,toc);
end
fid=fopen(fullfile(root,'unified','count_matching.json'),'w');fprintf(fid,'%s',jsonencode(struct('units_per_recording',n,'seeds',seeds,'repetitions',3,'primary_neuron_filter','none','purpose','sensitivity only; each original recording retains all units in primary analysis')));fclose(fid);
end
