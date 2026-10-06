function restrict_cross_labels()
root=fileparts(fileparts(mfilename('fullpath')));source=fullfile(root,'unified','daily','stored_cross');dest=fullfile(root,'unified','daily','stored_cross_common_labels');if ~exist(dest,'dir'),mkdir(dest);end
fs=dir(fullfile(source,'*.mat'));
for k=1:numel(fs)
 S=load(fullfile(fs(k).folder,fs(k).name));result=S.result;
 keep=1:5;if strcmp(result.mode,'var'),keep=1:8;end
 vector_fields={'within_sdi','within_svm','signal_distance','noise_variance','rms_radius','signal_radius','signal_dimension','rms_radius_raw'};
 for f=vector_fields,result.(f{1})=result.(f{1})(keep);end
 result.between_sdi=result.between_sdi(keep,keep);result.between_svm=result.between_svm(keep,keep);result.category_names=result.category_names(keep);result.candidate='stored_cross_common_labels';
 save(fullfile(dest,fs(k).name),'result');
end
end
