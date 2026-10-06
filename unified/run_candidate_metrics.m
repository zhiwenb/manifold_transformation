function run_candidate_metrics(candidate)
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'unified'));
if ismember(candidate,{'canonical','stored_cross'}),input_root=fullfile(root,'data');else,input_root=fullfile(root,'unified','inputs',candidate);end
out=fullfile(root,'unified','daily',candidate);if ~exist(out,'dir'),mkdir(out);end
for si=1:5
 session=['s' num2str(si)];for mode={'var','global'}
  input_mode=mode{1};if strcmp(candidate,'stored_cross'),if strcmp(input_mode,'var'),input_mode='variation';else,input_mode='prototype';end;end
  fs=dir(fullfile(input_root,session,'fr',input_mode,'*.mat')); 
  if isempty(fs),fs=dir(fullfile(input_root,session,mode{1},'*.mat'));end
  for fi=1:numel(fs)
   destination=fullfile(out,[session '_' mode{1} '_' fs(fi).name]);
   if isfile(destination),fprintf('Cached %s\n',fs(fi).name);continue;end
   tic;result=compute_daily_metrics(fullfile(fs(fi).folder,fs(fi).name),session,mode{1});
   result.candidate=candidate;save(destination,'result');fprintf('%s %s %s %s %.1f sec\n',candidate,session,mode{1},result.date,toc);
  end
 end
end
end
