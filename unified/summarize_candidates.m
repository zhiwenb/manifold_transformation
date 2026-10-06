function summarize_candidates()
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'unified'); rows={};stages={};effects={};
candidates={'stored_cross','stored_cross_common_labels','canonical','union_native','union_common240','fixed_count_common240'}; datasets=struct();
for ci=1:numel(candidates)
 fs=dir(fullfile(out,'daily',candidates{ci},'*.mat'));arr={};
 for fi=1:numel(fs),S=load(fullfile(fs(fi).folder,fs(fi).name));arr{end+1}=S.result;end
 datasets.(candidates{ci})=arr;
end
% The matched cohort intersects recording dates across source families.
for ci=1:numel(candidates)
 for cohort={'matched','all'}
  if strcmp(cohort{1},'all') && ~ismember(candidates{ci},{'canonical','union_common240','fixed_count_common240'}),continue;end
  for rule={'thirds','halves','first_last','calendar_thirds'}
   tag=strjoin({candidates{ci},cohort{1},rule{1}},'__'); per={};
   for si=1:5
    session=['s' num2str(si)];a=datasets.(candidates{ci});
    calendar=[];
    for k=1:numel(a),if strcmp(a{k}.session,session),calendar(end+1)=date_number(a{k}.date);end;end
    calendar=unique(sort(calendar));cut=max(1,floor(numel(calendar)/3));
    for mode={'var','global'}
     rr=a(cellfun(@(x) strcmp(x.session,session)&&strcmp(x.mode,mode{1}),a));
     if strcmp(cohort{1},'matched')
      b=datasets.stored_cross;bd=cellfun(@(x) x.date,b(cellfun(@(x) strcmp(x.session,session)&&strcmp(x.mode,mode{1}),b)),'UniformOutput',false);
      rr=rr(cellfun(@(x) ismember(x.date,bd),rr));
     end
     if isempty(rr),continue;end
     nums=cellfun(@(x) date_number(x.date),rr);[nums,idx]=sort(nums);rr=rr(idx);n=numel(rr);w=max(1,floor(n/3));E=1:w;L=n-w+1:n;
     if strcmp(rule{1},'halves'),w=floor(n/2);E=1:w;L=n-w+1:n;end
     if strcmp(rule{1},'first_last'),E=1;L=n;end
     if strcmp(rule{1},'calendar_thirds'),E=find(nums<=calendar(cut));L=find(nums>=calendar(end-cut+1));end
     stages{end+1}=struct('scenario',tag,'session',session,'mode',mode{1},'available_dates',{{rr{:}}},'early_dates',{cellfun(@(x)x.date,rr(E),'UniformOutput',false)},'late_dates',{cellfun(@(x)x.date,rr(L),'UniformOutput',false)});
     stages{end}.available_dates=cellfun(@(x)x.date,rr,'UniformOutput',false);
     if isempty(E)||isempty(L)||any(ismember(E,L)),continue;end
     metrics=stage_metrics(rr,E,L,si,mode{1});
     for k=1:numel(metrics)
      rec=metrics{k};rec.scenario=tag;rec.candidate=candidates{ci};rec.cohort=cohort{1};rec.stage_rule=rule{1};rec.session=session;rec.animal='M2';if si<=2,rec.animal='M1';end;rec.mode=mode{1};per{end+1}=rec;effects{end+1}=rec;
     end
    end
   end
   if isempty(per),continue;end
   metrics=unique(cellfun(@(x)x.metric,per,'UniformOutput',false));
   for mi=1:numel(metrics)
    for group={'pooled','M1','M2'}
     pp=per(cellfun(@(x)strcmp(x.metric,metrics{mi})&&(strcmp(group{1},'pooled')||strcmp(x.animal,group{1})),per));
     if isempty(pp),continue;end
     v=cellfun(@(x)x.effect,pp);valid=isfinite(v);v=v(valid);pp=pp(valid);if isempty(v),continue;end
     r=statistics(v);r.scenario=tag;r.candidate=candidates{ci};r.cohort=cohort{1};r.stage_rule=rule{1};r.metric=metrics{mi};r.group=group{1};r.sessions=cellfun(@(x)x.session,pp,'UniformOutput',false);r.early=cellfun(@(x)x.early,pp);r.late=cellfun(@(x)x.late,pp);rows{end+1}=r;
    end
   end
  end
 end
end
% Reanalyse original pair observations, preserving their original stage selection.
ref={'fig3E','fig3F','fig4E','fig4F'};
for k=1:4
 T=readtable(fullfile(out,'reference_values',[ref{k} '.csv']));v=T.delta;
 pooled=statistics(v);pooled.metric=ref{k};pooled.unit='category_pair';refstats{k}.pooled=pooled;
 for group={'pooled','M1','M2'}
  TT=T;if ~strcmp(group{1},'pooled'),TT=TT(strcmp(TT.monkey,group{1}),:);end
  ss=unique(TT.session);sv=nan(numel(ss),1);for si=1:numel(ss),sv(si)=mean(TT.delta(strcmp(TT.session,ss{si})),'omitnan');end
  r=statistics(sv);r.scenario='historical_cache__historical_stage__session_statistics';r.candidate='historical_cache';r.cohort='historical';r.stage_rule='historical';r.metric=ref{k};r.group=group{1};r.sessions=ss;r.early=[];r.late=[];rows{end+1}=r;
 end
end
write_json(fullfile(out,'summary.json'),rows);write_json(fullfile(out,'session_effects.json'),effects);write_json(fullfile(out,'stage_manifest.json'),stages);write_json(fullfile(out,'reference_statistics.json'),refstats);
% CSV with scalar results; JSON retains full per-session vectors.
C=cell(numel(rows),12);
for k=1:numel(rows),r=rows{k};C(k,:)={r.scenario,r.metric,r.group,r.n,r.mean,r.sem,r.ttest_p,r.signflip_p,r.ci_low,r.ci_high,r.n_positive,r.n_negative};end
T=cell2table(C,'VariableNames',{'scenario','metric','group','n_sessions','mean','sem','ttest_p','exact_signflip_p','ci_low','ci_high','positive_sessions','negative_sessions'});writetable(T,fullfile(out,'summary.csv'));
fprintf('Saved %d summary comparisons and %d session effects\n',numel(rows),numel(effects));
end

function rr=stage_metrics(d,E,L,sid,mode)
rr={};names=d{1}.category_names;P=numel(names);
if strcmp(mode,'var'),ti=find(startsWith(names,'T'));ui=find(startsWith(names,'U'));else,ti=[3 4 5];ui=1;if ismember(sid,[2 4]),ti=4;ui=[1 2 3 5];end;end
for f={'within_sdi','within_svm','within_svm_loo','signal_distance','noise_variance','rms_radius','signal_radius','signal_dimension','rms_radius_raw'}
 field=f{1};if ~all(cellfun(@(x)isfield(x,field),d)),continue;end;if strcmp(mode,'global')&&ismember(field,{'within_sdi','within_svm','within_svm_loo','signal_distance','noise_variance'}),continue;end
 A=cell2mat(cellfun(@(x)x.(field),d,'UniformOutput',false)');
 good=all(isfinite(A([E L],:)),1);tt=ti(good(ti));uu=ui(good(ui));if isempty(tt)||isempty(uu),continue;end
 eT=mean(A(E,tt),'all');eU=mean(A(E,uu),'all');lT=mean(A(L,tt),'all');lU=mean(A(L,uu),'all');
 effect=(lT-eT)-(lU-eU);
 if ismember(field,{'within_sdi','within_svm','within_svm_loo','signal_distance','noise_variance'}),if abs(eT)<eps||abs(eU)<eps,continue;end;effect=100*((lT-eT)/eT-(lU-eU)/eU);end
 rr{end+1}=struct('metric',[mode '_' field '_DiD'],'effect',effect,'early',eT-eU,'late',lT-lU);
end
for f={'between_sdi','between_svm'}
 field=f{1};V=nan(numel(d),P*(P-1)/2);mask=triu(true(P),1);
 for k=1:numel(d),V(k,:)=d{k}.(field)(mask);end
 good=all(isfinite(V([E L],:)),1);if ~any(good),continue;end
 early=mean(V(E,good),'all');late=mean(V(L,good),'all');
 rr{end+1}=struct('metric',[mode '_' field],'effect',late-early,'early',early,'late',late);
end
if strcmp(mode,'global')
 A=cell2mat(cellfun(@(x)x.orthogonality,d,'UniformOutput',false)');e=mean(A(E,:),1);l=mean(A(L,:),1);
 rr{end+1}=struct('metric','global_orthogonality_Polar_minus_Grating','effect',(l(1)-e(1))-(l(2)-e(2)),'early',e(1)-e(2),'late',l(1)-l(2));
 if ismember(sid,[1 5]),rr{end+1}=struct('metric','global_orthogonality_Hyper_minus_Grating','effect',(l(3)-e(3))-(l(2)-e(2)),'early',e(3)-e(2),'late',l(3)-l(2));end
 A=cell2mat(cellfun(@(x)x.rsa,d,'UniformOutput',false)');e=mean(A(E,:),1);l=mean(A(L,:),1);rel=(l-e)./(l+e);rel(abs(l+e)<eps)=NaN;
 rr{end+1}=struct('metric','global_RSA_Polar_net','effect',rel(1)-rel(2),'early',e(1),'late',l(1));
 if ismember(sid,[1 5]),rr{end+1}=struct('metric','global_RSA_Hyper_net','effect',rel(3)-rel(2),'early',e(3),'late',l(3));end
end
end

function r=statistics(v)
v=v(isfinite(v));n=numel(v);r=struct('n',n,'mean',mean(v),'sem',NaN,'ttest_p',NaN,'signflip_p',NaN,'ci_low',NaN,'ci_high',NaN,'effects',v,'n_positive',sum(v>0),'n_negative',sum(v<0));
if n>1,[~,r.ttest_p,ci]=ttest(v,0);r.ci_low=ci(1);r.ci_high=ci(2);r.sem=std(v)/sqrt(n);end
if n<=15&&n>0,signs=2*double(dec2bin(0:2^n-1,n)-'0')-1;null=mean(signs.*v(:)',2);r.signflip_p=mean(abs(null)>=abs(mean(v))-1e-12);end
end
function n=date_number(s)
n=datenum(2000+str2double(s(5:6)),str2double(s(1:2)),str2double(s(3:4)));
end
function write_json(path,value)
fid=fopen(path,'w');fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));fclose(fid);
end
