function write_manifold_benchmark()
root=fileparts(fileparts(mfilename('fullpath')));S=load(fullfile(root,'unified','manifold_benchmark.mat'));
r=struct('session','s1','mode','var','date','041116','candidate','union_common240','n_t',S.options.n_t,'elapsed_seconds',S.elapsed,'output',S.output,'scope','one-recording full benchmark; not cohort-level validation');
fid=fopen(fullfile(root,'unified','manifold_benchmark.json'),'w');fprintf(fid,'%s',jsonencode(r,PrettyPrint=true));fclose(fid);
end
