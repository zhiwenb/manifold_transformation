function rebuild_candidates(raw_root)
% Rebuild identical category layouts from raw CDT using original stored units.
% Union adds only units already present in one historical input version.
root=fileparts(fileparts(mfilename('fullpath')));
legacy=[1 4 5 6 7]; sessions={'s1','s2','s3','s4','s5'}; checks={}; kk=0;
for si=1:5
 for mode={'var','global'}
  kind=mode{1}; alt='variation'; if strcmp(kind,'global'),alt='prototype';end
  fs=dir(fullfile(root,'data',sessions{si},'fr',kind,'*.mat'));
  for fi=1:numel(fs)
   T=load(fullfile(fs(fi).folder,fs(fi).name)); m=T.meta; date=regexp(fs(fi).name,'\d{6}','match','once');
   af=dir(fullfile(root,'data',sessions{si},'fr',alt,['*' date '.mat'])); units=m.neurons_used';
   if ~isempty(af),A=load(fullfile(af(1).folder,af(1).name),'meta');units=unique([units;A.meta.neurons_used'],'rows','stable');end
   raw=fullfile(raw_root,['s' num2str(legacy(si))],kind,m.file); assert(isfile(raw),'Missing CDT: %s',raw);
   R=load(raw,'CDTTables');c=R.CDTTables{1};
   % Use exactly the release input condition grouping and repetition count.
   ids=cell(1,numel(T.FR_by_category)); trial_ids=[];
   for g=1:numel(ids)
    name=m.category_names{g};
    if isempty(T.FR_by_category{g}),ids{g}=[];continue;end
    ids{g}=m.category_cond_idx_original.(name);
    assert(numel(ids{g})==size(T.FR_by_category{g},2),'Condition layout mismatch');
    for cc=ids{g},ix=find(c.condition==cc);assert(numel(ix)>=m.num_trials_balanced,'Trial deficit');trial_ids=[trial_ids;ix(1:m.num_trials_balanced)];end
   end
   trial_ids=[]; max_condition=max(c.condition);
   if legacy(si)==4 && strcmp(kind,'var'),max_condition=min(108,max_condition);end
   for cc=1:max_condition,ix=find(c.condition==cc);trial_ids=[trial_ids;ix(1:min(numel(ix),m.num_trials_balanced))];end
   trial_ids=unique(trial_ids,'stable');
   rates=zeros(size(units,1),numel(c.condition),3); base=zeros(size(units,1),numel(c.condition));
   for ti=trial_ids'
    ele=c.spikeElectrode{ti};un=c.spikeUnit{ti};tm=c.spikeTimes{ti};
    for ni=1:size(units,1)
     times=tm(ele==units(ni,1)&un==units(ni,2));if iscell(times),times=horzcat(times{:});end
     base(ni,ti)=nnz(times>.350 & times<.500)/.150;
     rates(ni,ti,1)=nnz(times>m.stim_on_ms/1000 & times<m.stim_off_ms/1000)/((m.stim_off_ms-m.stim_on_ms)/1000);
     rates(ni,ti,2)=nnz(times>.630 & times<m.stim_off_ms/1000)/((m.stim_off_ms-630)/1000);
     rates(ni,ti,3)=nnz(times>.630 & times<.870)/.240;
    end
   end
   % Match the historical pooled baseline subtraction within each recording.
   rates=rates-mean(base(:,trial_ids),2);
   variants={'release_rebuild','union_native','union_common240'};
   for vi=1:3
    FR_by_category=cell(size(T.FR_by_category));
    for g=1:numel(ids)
     if isempty(ids{g}),continue;end
     X=zeros(size(units,1),numel(ids{g}),m.num_trials_balanced,'single');
     for ci=1:numel(ids{g}),ix=find(c.condition==ids{g}(ci));X(:,ci,:)=reshape(rates(:,ix(1:m.num_trials_balanced),vi),size(units,1),1,[]);end
     if vi==1,[~,ii]=ismember(m.neurons_used',units,'rows');X=X(ii,:,:);end
     FR_by_category{g}=X;
    end
    meta=m;
    if vi~=1,meta.neurons_used=units';meta.num_neurons=size(units,1);meta.stim_on_ms=630;if vi==3,meta.stim_off_ms=870;end;end
    if vi>1
     folder=fullfile(root,'unified','inputs',variants{vi},sessions{si},kind);if ~exist(folder,'dir'),mkdir(folder);end
     save(fullfile(folder,fs(fi).name),'FR_by_category','meta','-v7.3');
    end
    if vi==1
     errs=cellfun(@(x,y) max(abs(double(x(:))-double(y(:))),[],'omitnan'),FR_by_category,T.FR_by_category,'UniformOutput',false);
     errs=[errs{:}];kk=kk+1;checks{kk}=struct('session',sessions{si},'mode',kind,'date',date,'max_abs_fr_difference',max(errs),'release_units',size(m.neurons_used,2),'union_units',size(units,1),'raw_file',m.file);
    end
   end
   fprintf('%s %s %s: original %d, union %d, rebuild error %.6g\n',sessions{si},kind,date,size(m.neurons_used,2),size(units,1),checks{kk}.max_abs_fr_difference);
  end
 end
end
fid=fopen(fullfile(root,'unified','raw_rebuild_checks.json'),'w');fprintf(fid,'%s',jsonencode(checks));fclose(fid);
end
