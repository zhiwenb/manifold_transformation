function recompute_loo(candidate)
% Retain original Fig2 leave-one-out estimator on precisely the same daily FR.
root=fileparts(fileparts(mfilename('fullpath')));ids={[7 8 9 12 13 14 17 18 19],[3 11 12 13 21],10:18,[3 11 12 13 21],[7 8 9 12 13 14 17 18 19]};
fs=dir(fullfile(root,'unified','daily',candidate,'*var*.mat'));
for f=1:numel(fs)
 path=fullfile(fs(f).folder,fs(f).name);S=load(path);result=S.result;if isfield(result,'loo_seed_version')&&result.loo_seed_version==1,continue;end
 source_file=result.file;if ~isfile(source_file),source_file=fullfile(root,result.file);end;D=load(source_file);F=D.FR_by_category;sid=str2double(result.session(2:end));vv=nan(1,numel(F));tic;
 for g=1:numel(F)
  X=double(F{g});if isempty(X),continue;end;[N,C,R]=size(X);ii=ids{sid};ii=ii(ii<=C);pairs=nchoosek(ii,2);av=nan(size(pairs,1),1);
  for p=1:size(pairs,1)
   rng(42+1000*g+p,'twister');A=reshape(X(:,pairs(p,1),:),N,R)';B=reshape(X(:,pairs(p,2),:),N,R)';A(any(isnan(A),2),:)=[];B(any(isnan(B),2),:)=[];features=[A;B];labels=[ones(size(A,1),1);-ones(size(B,1),1)];pred=zeros(size(labels));
   for k=1:numel(labels)
    train=true(size(labels));train(k)=false;model=fitcsvm(features(train,:),labels(train),'KernelFunction','linear','Standardize',true);pred(k)=predict(model,features(k,:));
   end
   av(p)=100*mean(pred==labels);
  end
  vv(g)=mean(av,'omitnan');
 end
 result.within_svm_loo=vv;result.loo_seed_version=1;save(path,'result');fprintf('LOO %s %s %s %.1f sec\n',candidate,result.session,result.date,toc);
end
end
