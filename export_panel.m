function export_panel(fh, name, outdir)
% EXPORT_PANEL  Write one panel as vector PDF + 300 dpi PNG.
if nargin < 3 || isempty(outdir), outdir = fullfile('..','output'); end
if ~exist(outdir,'dir'), mkdir(outdir); end
exportgraphics(fh, fullfile(outdir,[name '.pdf']), 'ContentType','vector');
exportgraphics(fh, fullfile(outdir,[name '.png']), 'Resolution',300);
fprintf('  wrote %s.pdf / .png\n', name);
end

