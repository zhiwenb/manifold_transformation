function finalize_metadata()
% Convert historical local paths to portable paths; restore sample dates.
root=fileparts(fileparts(mfilename('fullpath')));
fs=dir(fullfile(root,'unified','daily','**','*.mat'));
for k=1:numel(fs)
 p=fullfile(fs(k).folder,fs(k).name);S=load(p);result=S.result;
 result.file=portable(result.file,root);result.date=regexp(fs(k).name,'\d{6}(?=\.mat)','match','once');assert(~isempty(result.date));save(p,'result');
end
fs=dir(fullfile(root,'unified','fixed_count_draws','*.mat'));
for k=1:numel(fs)
 p=fullfile(fs(k).folder,fs(k).name);S=load(p);for d=1:numel(S.draws),S.draws{d}.file=portable(S.draws{d}.file,root);S.draws{d}.date=regexp(fs(k).name,'\d{6}(?=\.mat)','match','once');end;save(p,'-struct','S');
end
end
function p=portable(p,root)
if startsWith(p,[root filesep]),p=p(numel(root)+2:end);end
end
