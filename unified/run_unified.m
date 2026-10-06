function run_unified(rebuild, raw_root)
% Recompute the complete comparison from packaged FR, without SNR screening.
% Set rebuild=true and supply the original CDT root to rebuild candidates.
if nargin<1,rebuild=false;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'unified'),fullfile(root,'analysis'));
if rebuild,assert(nargin==2,'Provide raw CDT root.');rebuild_candidates(raw_root);end
for candidate={'canonical','stored_cross','union_native','union_common240'},run_candidate_metrics(candidate{1});end
for candidate={'canonical','union_common240'},recompute_loo(candidate{1});end
restrict_cross_labels();fixed_count_sensitivity();finalize_metadata();summarize_candidates();fig2_baseline_audit();summarize_manifold_cache();plot_comparisons();
end
