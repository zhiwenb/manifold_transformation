function plot_comparisons()
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'unified','figures');if ~exist(out,'dir'),mkdir(out);end
R=jsondecode(fileread(fullfile(root,'unified','summary.json')));
metrics={'var_within_sdi_DiD','var_within_svm_DiD','var_between_svm','global_between_svm','global_orthogonality_Polar_minus_Grating','global_RSA_Polar_net'};
titles={'Within-category SDI DiD (%)','Within-category SVM DiD (%; 5-fold)','Variation category decoding (pp)','Global category decoding (pp)','Polar - Grating angle change (deg)','Polar net RSA change'};
candidates={'stored_cross','canonical','union_native','union_common240'};labels={'Cross caches FR','Canonical FR','Union, native','Union, 240 ms'};
f=figure('Color','w','Position',[50 50 1000 950]);tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
for k=1:6
 nexttile;mu=nan(1,4);se=mu;
 for c=1:4
  mask=arrayfun(@(x) strcmp(x.scenario,[candidates{c} '__matched__thirds'])&&strcmp(x.metric,metrics{k})&&strcmp(x.group,'pooled'),R);a=R(mask);mu(c)=a.mean;se(c)=a.sem;
 end
 bar(1:4,mu,'FaceColor',[.3 .5 .7]);hold on;errorbar(1:4,mu,se,'k.');yline(0,'k:');xticks(1:4);xticklabels(labels);xtickangle(35);title(titles{k});box off;set(gca,'FontSize',12);
end
sgtitle('Matched recording dates; early/late thirds; session mean +/- SEM');exportgraphics(f,fullfile(out,'source_comparison.pdf'),'ContentType','vector');close(f);
metrics={'var_rms_radius_DiD','var_signal_radius_DiD','global_rms_radius_DiD','global_signal_radius_DiD','var_signal_distance_DiD','var_noise_variance_DiD'};
titles={'Variation RMS-radius DiD (Hz/unit RMS)','Variation signal-radius DiD','Global RMS-radius DiD','Global signal-radius DiD','Trained-pair distance DiD (%)','Trial-noise variance DiD (%)'};
f=figure('Color','w','Position',[50 50 1000 950]);tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
for k=1:6
 nexttile;mu=nan(1,4);se=mu;for c=1:4,mask=arrayfun(@(x) strcmp(x.scenario,[candidates{c} '__matched__thirds'])&&strcmp(x.metric,metrics{k})&&strcmp(x.group,'pooled'),R);a=R(mask);mu(c)=a.mean;se(c)=a.sem;end
 bar(1:4,mu,'FaceColor',[.35 .6 .45]);hold on;errorbar(1:4,mu,se,'k.');yline(0,'k:');xticks(1:4);xticklabels(labels);xtickangle(35);title(titles{k});box off;set(gca,'FontSize',12);
end
sgtitle('Direct geometry is a separate estimator from MFT radius / capacity');exportgraphics(f,fullfile(out,'geometry_comparison.pdf'),'ContentType','vector');close(f);
metrics={'var_within_sdi_DiD','var_within_svm_loo_DiD','var_between_svm','global_between_sdi','global_orthogonality_Polar_minus_Grating','global_RSA_Polar_net'};
titles={'Within-category SDI DiD (%)','Within-category SVM DiD (%; LOO)','Variation decoding change (pp)','Global SDI change','Polar - Grating angle change (deg)','Polar net RSA change'};
f=figure('Color','w','Position',[50 50 1000 950]);tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
for k=1:6
 nexttile;a=R(arrayfun(@(x)strcmp(x.scenario,'union_common240__all__thirds')&&strcmp(x.metric,metrics{k})&&strcmp(x.group,'pooled'),R));v=a.effects;
 colors=[repmat([.2 .4 .8],2,1);repmat([.8 .3 .2],3,1)];scatter(1:numel(v),v,55,colors,'filled');hold on;yline(0,'k:');yline(mean(v),'k--');xticks(1:5);xticklabels({'s1','s2','s3','s4','s5'});title(titles{k});box off;set(gca,'FontSize',12);
end
sgtitle('Shared 240-ms window, original-unit union, all recordings, stage thirds');exportgraphics(f,fullfile(out,'primary_session_effects.pdf'),'ContentType','vector');close(f);
B=jsondecode(fileread(fullfile(root,'unified','fig2_baseline_summary.json')));f=figure('Color','w','Position',[50 50 850 400]);tiledlayout(1,2);
for k=1:2
 nexttile;tasks={'sdi','svm'};rule={'cache_baseline','fresh_cache_baseline','fresh_earliest_fr_baseline'};values=nan(1,3);
 for j=1:3,a=B(arrayfun(@(x)strcmp(x.metric,tasks{k})&&strcmp(x.baseline_rule,rule{j})&&strcmp(x.group,'M2'),B));values(j)=a.day_mean;end
 bar(values,'FaceColor',[.65 .35 .3]);yline(0,'k:');xticks(1:3);xticklabels({'Original cache','Fresh, old baseline','Fresh, earliest FR'});xtickangle(20);title(['M2 ' upper(tasks{k}) ' DiD (%)']);box off;
end
sgtitle('Same post dates; changing the baseline reverses M2 improvement');exportgraphics(f,fullfile(out,'baseline_effect.pdf'),'ContentType','vector');close(f);
end
