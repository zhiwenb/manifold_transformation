function benchmark_manifold()
% Full estimator benchmark only; no whole-cohort capacity inference.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'analysis'));
options=struct('kappa',0,'n_t',1000,'flag_NbyM',1,'center_scale',1);rng(42);benchmark_timer=tic;
output=analyze_manifold_file(fullfile(root,'unified','inputs','union_common240','s1','var','FR_s1_var_041116.mat'),options);elapsed=toc(benchmark_timer);
save(fullfile(root,'unified','manifold_benchmark.mat'),'output','options','elapsed');
write_manifold_benchmark();
end
