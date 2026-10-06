function run_figures(panels)
% Draw available figure panels from original FR and saved historical results.
% Examples: run_figures({'fig4BCD','fig4E','fig4F','fig5B','fig5C'})
%           run_figures('fig4BCD_mixed')
root = fileparts(mfilename('fullpath'));
addpath(root);
prepare_inputs();
if nargin == 0
    panels = {'fig2D','fig2E','fig2F','fig3A','fig3BCD','fig3E','fig3F', ...
              'fig4BCD','fig4E','fig4F','fig5B','fig5C'};
end
if ischar(panels) || isstring(panels), panels = cellstr(panels); end
for k = 1:numel(panels)
    panel = char(panels{k});
    assert(~isempty(regexp(panel, '^fig[2-5][A-F]+(_mixed|_global|_recomputed|_delta|_checked)?$', 'once')), ...
        'Invalid panel name.');
    script = fullfile(root, 'panels', [panel '.m']);
    assert(isfile(script), 'Unavailable panel: %s', panel);
    destination = figure_output(panel);
    fprintf('Drawing %s from historical inputs.\n', panel);
    execute_panel(script, destination);
end
end

function execute_panel(script, destination)
close all; % Export only figures created by this panel.
% Isolate legacy scripts that clear their workspace and use relative exports.
previous = pwd;
cleanup = onCleanup(@() cd(previous));
cd(destination);
temporary = [tempname(destination) '.m'];
fid = fopen(temporary,'w'); assert(fid >= 0); fprintf(fid,'%s',fileread(script)); fclose(fid);
script_cleanup = onCleanup(@() delete(temporary));
run_isolated(temporary);
figures = findall(groot, 'Type', 'figure');
for k = 1:numel(figures)
    exportgraphics(figures(k), fullfile(destination, sprintf('panel_%d.png',k)), 'Resolution',300);
    exportgraphics(figures(k), fullfile(destination, sprintf('panel_%d.pdf',k)), 'ContentType','vector');
end
end

function run_isolated(script)
% The temporary script is in the output directory; clear is isolated here.
run(script);
end
