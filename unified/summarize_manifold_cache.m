function summarize_manifold_cache()
% Stage sensitivity of historical daily estimates; not raw-input validated MFT.
root=fileparts(fileparts(mfilename('fullpath')));rows={};skips={};
for mode={'var','global'}
 for rule={'thirds','halves','first_last'}
  per={};
  for si=1:5
   session=['s' num2str(si)];S=load(fullfile(root,'unified','legacy_manifold',[session '_' mode{1} '.mat']));D=S.session_results.files;D=D(arrayfun(@(x)~isempty(x.filename)&&(ischar(x.filename)||isstring(x.filename)),D));
   nums=arrayfun(@(x)date_num(regexp(x.filename,'\d{6}','match','once')),D);[~,ii]=sort(nums);D=D(ii);n=numel(D);w=max(1,floor(n/3));E=1:w;L=n-w+1:n;
   if strcmp(rule{1},'halves'),w=floor(n/2);E=1:w;L=n-w+1:n;end
   if strcmp(rule{1},'first_last'),E=1;L=n;end
   names={};values=cell(n,1);
   for d=1:n
    F=D(d).manifold_details;
    if strcmp(mode{1},'var')
     % Original estimator drops empty classes; infer retained labels from FR.
     date=regexp(D(d).filename,'\d{6}','match','once');fs=dir(fullfile(root,'data',session,'fr','var',['*' date '.mat']));
     if isempty(fs),continue;end
     X=load(fullfile(fs(1).folder,fs(1).name));labels=X.meta.category_names(~cellfun(@isempty,X.FR_by_category));
    else,labels={'G','H','R','S','T'};end
    if numel(labels)~=numel(F),skips{end+1}=struct('session',session,'mode',mode{1},'file',D(d).filename,'reason','category layout differs');continue;end
    values{d}=[[F.R_M];[F.D_M];[F.a_M]]';names=labels;
   end
   if any(cellfun(@isempty,values([E L]))),continue;end
   if numel(unique(cellfun(@(x)size(x,1),values([E L]))))~=1,continue;end
   if strcmp(mode{1},'var'),ti=find(startsWith(names,'T'));ui=find(startsWith(names,'U'));else,ti=[3 4 5];ui=1;if ismember(si,[2 4]),ti=4;ui=[1 2 3 5];end;end
   A=cat(3,values{E});B=cat(3,values{L});e=mean(A,3);l=mean(B,3);v=mean(l(ti,:)-e(ti,:),1)-mean(l(ui,:)-e(ui,:),1);
   for k=1:3,per{end+1}=struct('session',session,'metric',k,'effect',v(k));end
  end
  for k=1:3
   v=cellfun(@(x)x.effect,per(cellfun(@(x)x.metric==k,per)));v=v(isfinite(v));if isempty(v),continue;end
   p=NaN;if numel(v)>1,[~,p]=ttest(v);end
   signs=2*double(dec2bin(0:2^numel(v)-1,numel(v))-'0')-1;sp=mean(abs(mean(signs.*v(:)',2))>=abs(mean(v))-1e-12);
   labels={'radius','dimension','capacity'};rows{end+1}=struct('mode',mode{1},'rule',rule{1},'metric',labels{k},'n',numel(v),'mean',mean(v),'effects',v,'ttest_p',p,'signflip_p',sp,'lineage','historical daily MFT caches; raw input identity not verified');
  end
 end
end
fid=fopen(fullfile(root,'unified','manifold_cache_sensitivity.json'),'w');fprintf(fid,'%s',jsonencode(rows,PrettyPrint=true));fclose(fid);
fid=fopen(fullfile(root,'unified','manifold_cache_skips.json'),'w');fprintf(fid,'%s',jsonencode(skips,PrettyPrint=true));fclose(fid);
end
function n=date_num(s),n=datenum(2000+str2double(s(5:6)),str2double(s(1:2)),str2double(s(3:4)));end
