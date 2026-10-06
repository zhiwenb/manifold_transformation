function inventory_inputs()
root=fileparts(fileparts(mfilename('fullpath'))); rows=struct([]); k=0;
sessions={'s1','s2','s3','s4','s5','extra_var'}; modes={'var','global','variation','prototype'};
for si=1:numel(sessions)
 for mi=1:numel(modes)
  fs=dir(fullfile(root,'data',sessions{si},'fr',modes{mi},'*.mat'));
  for fi=1:numel(fs)
   S=load(fullfile(fs(fi).folder,fs(fi).name)); k=k+1;
   rows(k).session=sessions{si}; rows(k).mode=modes{mi}; rows(k).file=fs(fi).name;
   rows(k).date=regexp(fs(fi).name,'\d{6}','match','once');
   rows(k).fields=fieldnames(S); rows(k).sizes=cellfun(@size,S.FR_by_category,'UniformOutput',false);
   if isfield(S,'meta'), rows(k).meta=S.meta; else,rows(k).meta=struct(); end
  end
 end
end
fid=fopen(fullfile(root,'unified','input_inventory.json'),'w'); fprintf(fid,'%s',jsonencode(rows)); fclose(fid);
S=load(fullfile(root,'data','s1','fr','var','FR_s1_var_041116.mat')); disp(S.meta);
end
