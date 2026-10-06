function fig2_baseline_audit()
root=fileparts(fileparts(mfilename('fullpath')));out={};checks={};postidx={ [3 4 5 6],[2 3],[2 3],[],[3 4 5 6]};sidlist=[1 2 3 5];
ids={[7 8 9 12 13 14 17 18 19],[3 11 12 13 21],10:18,[3 11 12 13 21],[7 8 9 12 13 14 17 18 19]};
for task={'sdi','svm'}
 metric=task{1};for si=sidlist
  session=['s' num2str(si)];fs=dir(fullfile(root,'data',session,'metrics',['fig2_' metric],'*.mat'));
  dates=arrayfun(@(x)regexp(x.name,'\d{6}','match','once'),fs,'UniformOutput',false);nums=cellfun(@date_num,dates);[~,order]=sort(nums);fs=fs(order);dates=dates(order);daily=nan(numel(fs),8);
  for k=1:numel(fs)
   S=load(fullfile(fs(k).folder,fs(k).name));if strcmp(metric,'sdi'),R=S.SDI_results;field='SDI_matrix';else,R=S.SVM_results;field='SVM_matrix';end
   labels={'T1','T2','T3','T4','U1','U2','U3','U4'};
   for g=1:8
    if ~isfield(R,labels{g}),continue;end
    A=R.(labels{g}).(field);pp=nchoosek(ids{si},2);val=A(sub2ind(size(A),pp(:,1),pp(:,2)));daily(k,g)=mean(val,'omitnan');
   end
  end
  if strcmp(metric,'svm'),daily=100*daily;end
  df=dir(fullfile(root,'unified','daily','canonical',[session '_var_*.mat']));dd=cell(1,numel(df));for k=1:numel(df),S=load(fullfile(df(k).folder,df(k).name));dd{k}=S.result;end
  nums=cellfun(@(x)date_num(x.date),dd);[~,ii]=sort(nums);dd=dd(ii);
  field='within_sdi';if strcmp(metric,'svm'),field='within_svm_loo';end
  assert(all(cellfun(@(x)isfield(x,field),dd)),'LOO not yet complete.');
  first=dd{1}.(field);baseline=daily(1,:);t=1:2;u=5:6;if si==5,t=1:4;u=5:8;end
  for k=1:numel(fs)
   ind=find(cellfun(@(x)strcmp(x.date,dates{k}),dd));if ~isempty(ind),checks{end+1}=struct('metric',metric,'session',session,'date',dates{k},'max_abs_difference',max(abs(daily(k,:)-dd{ind}.(field)),[],'omitnan'));end
  end
  for b={'cache_baseline','fresh_cache_baseline','fresh_earliest_fr_baseline'}
   base=baseline;basedate=dates{1};daily_used=daily;
   if ~strcmp(b{1},'cache_baseline')
    for dk=1:numel(dates),ix=find(cellfun(@(x)strcmp(x.date,dates{dk}),dd));daily_used(dk,:)=dd{ix}.(field);end
    base=daily_used(1,:);
   end
   if strcmp(b{1},'fresh_earliest_fr_baseline'),base=first;basedate=dd{1}.date;end
   for k=postidx{si}
    if k>size(daily,1),continue;end
    eT=mean(base(t),'omitnan');eU=mean(base(u),'omitnan');lT=mean(daily_used(k,t),'omitnan');lU=mean(daily_used(k,u),'omitnan');effect=100*((lT-eT)/eT-(lU-eU)/eU);
    out{end+1}=struct('metric',metric,'session',session,'animal',sprintf('M%d',1+(si>2)),'baseline_rule',b{1},'baseline_date',basedate,'post_date',dates{k},'effect',effect);
   end
  end
 end
end
summ={};for metric={'sdi','svm'},for rule={'cache_baseline','fresh_cache_baseline','fresh_earliest_fr_baseline'},for group={'M1','M2','pooled'}
 a=out(cellfun(@(x)strcmp(x.metric,metric{1})&&strcmp(x.baseline_rule,rule{1})&&(strcmp(group{1},'pooled')||strcmp(x.animal,group{1})),out));v=cellfun(@(x)x.effect,a);[~,p]=ttest(v);ss=unique(cellfun(@(x)x.session,a,'UniformOutput',false));sv=nan(1,numel(ss));for k=1:numel(ss),sv(k)=mean(v(cellfun(@(x)strcmp(x.session,ss{k}),a)));end
 sp=NaN;if numel(sv)>1,[~,sp]=ttest(sv);end
 summ{end+1}=struct('metric',metric{1},'baseline_rule',rule{1},'group',group{1},'n_days',numel(v),'day_mean',mean(v),'day_ttest_p',p,'n_sessions',numel(sv),'session_mean',mean(sv),'session_ttest_p',sp,'session_effects',sv);
end;end;end
write(fullfile(root,'unified','fig2_baseline_effects.json'),out);write(fullfile(root,'unified','fig2_baseline_summary.json'),summ);write(fullfile(root,'unified','fig2_cache_checks.json'),checks);
end
function n=date_num(s),n=datenum(2000+str2double(s(5:6)),str2double(s(1:2)),str2double(s(3:4)));end
function write(p,v),fid=fopen(p,'w');fprintf(fid,'%s',jsonencode(v,PrettyPrint=true));fclose(fid);end
