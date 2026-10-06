function result = compute_daily_metrics(input_file, session, mode)
% Identical estimators across candidates; CV=5 folds, fixed seed per pair.
S=load(input_file);root=fileparts(fileparts(mfilename('fullpath')));portable_file=input_file;if startsWith(input_file,[root filesep]),portable_file=input_file(numel(root)+2:end);end;F=S.FR_by_category;m=S.meta;names=m.category_names;
result=struct('file',portable_file,'session',session,'mode',mode,'date',regexp(input_file,'\d{6}(?=\.mat)','match','once'),'neurons',m.num_neurons,'trials',m.num_trials_balanced,'category_names',{names});
P=numel(F);result.between_sdi=nan(P);result.between_svm=nan(P);result.within_sdi=nan(1,P);result.within_svm=nan(1,P);result.signal_distance=nan(1,P);result.noise_variance=nan(1,P);result.rms_radius=nan(1,P);result.signal_radius=nan(1,P);result.signal_dimension=nan(1,P);result.rms_radius_raw=nan(1,P);
mu=cell(1,P);vv=cell(1,P);
trained_ids={ [7 8 9 12 13 14 17 18 19], [3 11 12 13 21], 10:18, [3 11 12 13 21], [7 8 9 12 13 14 17 18 19]};sid=str2double(session(2:end));
for g=1:P
 X=double(F{g});if isempty(X),continue;end
 [N,C,R]=size(X);cm=mean(X,3);cv=var(X,0,3);mu{g}=mean(cm,2);vv{g}=mean(cv,2);
 cloud=reshape(X,N,[]);cloud=cloud-mean(cloud,2);
 result.rms_radius_raw(g)=sqrt(mean(sum(cloud.^2,1)));result.rms_radius(g)=result.rms_radius_raw(g)/sqrt(N);
 centers=cm-mean(cm,2);mean_radius2=mean(sum(centers.^2,1));noise_bias=(1-1/C)*sum(mean(cv,2))/R;
 result.signal_radius(g)=sqrt(max(mean_radius2-noise_bias,0)/N);
 if C>1
  cov_signal=centers*centers'/(C-1)-diag(mean(cv,2)/R);eigv=max(eig((cov_signal+cov_signal')/2),0);result.signal_dimension(g)=sum(eigv)^2/(sum(eigv.^2)+eps);
 end
 result.noise_variance(g)=mean(cv,'all');
 if strcmp(mode,'var')
  ids=trained_ids{sid};ids=ids(ids<=C);pairs=nchoosek(ids,2);sdi=nan(size(pairs,1),1);acc=sdi;dist=sdi;
  for p=1:size(pairs,1)
   i=pairs(p,1);j=pairs(p,2);pv=(cv(:,i)+cv(:,j))/2;pv=max(pv,eps);sdi(p)=sqrt(sum((cm(:,i)-cm(:,j)).^2./pv));dist(p)=norm(cm(:,i)-cm(:,j));
   acc(p)=decode_pair(reshape(X(:,i,:),N,R)',reshape(X(:,j,:),N,R)',42+p);
  end
  result.within_sdi(g)=mean(sdi,'omitnan');result.within_svm(g)=mean(acc,'omitnan');result.signal_distance(g)=mean(dist,'omitnan');
 end
end
for i=1:P
 if isempty(F{i}),continue;end
 for j=i+1:P
  if isempty(F{j}),continue;end
  pv=max((vv{i}+vv{j})/2,eps);result.between_sdi(i,j)=sqrt(sum((mu{i}-mu{j}).^2./pv));result.between_sdi(j,i)=result.between_sdi(i,j);
  Ni=size(F{i},1);Ai=reshape(double(F{i}),Ni,[])';Aj=reshape(double(F{j}),Ni,[])';
  result.between_svm(i,j)=decode_pair(Ai,Aj,4200+100*i+j);result.between_svm(j,i)=result.between_svm(i,j);
 end
end
if strcmp(mode,'global')
 result.orthogonality=orthogonality_daily(F(1:5));result.rsa=rsa_daily(F(1:5));
else,result.orthogonality=nan(1,3);result.rsa=nan(1,3);end
end

function acc=decode_pair(A,B,seed)
A=A(all(isfinite(A),2),:);B=B(all(isfinite(B),2),:);n=min(size(A,1),size(B,1));acc=NaN;if n<5,return;end
rng(seed,'twister');A=A(randperm(size(A,1),n),:);B=B(randperm(size(B,1),n),:);X=[A;B];Y=[ones(n,1);-ones(n,1)];
model=fitcsvm(X,Y,'KernelFunction','linear','Standardize',true,'CrossVal','on','KFold',5);acc=100*(1-kfoldLoss(model));
end
