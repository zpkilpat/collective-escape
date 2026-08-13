function out = make_figS6(csvfile, xlsfile, outdir)
%MAKE_FIGS6  Area scaling and shoal density: what Data S2 does and does not
%            license.
%
%  A  first response time against shoal area, log-log, with the file|bout
%     clustered interval on the slope. The point estimate is negative, as a
%     first-passage minimum requires, but the clustered interval crosses
%     zero, so this is a SIGN CHECK and not a bound on M.
%  B  Data S2, headcount against area, pooled and within site. Pooled the
%     slope looks like constant density, but that is entirely between-site;
%     within a site area does not track headcount at all. Shoal area
%     measures spread, not pool size.
%  C  density from the nearest-neighbor distance rather than from the
%     tracked count, and the resulting physical headcount against Mhat.
%     Ntags_total/area is NOT a density: it contradicts the NND in the same
%     file by two orders of magnitude, so Ntags is a tracked subset.
%
%  Usage:  out = make_figS6('sciadv_adt8600_data_s1.csv', ...
%                           'sciadv_adt8600_data_s2.xlsx', 'output');
%
%  ZPK 2026

if nargin < 1 || isempty(csvfile), csvfile = 'sciadv_adt8600_data_s1.csv';  end
if nargin < 2 || isempty(xlsfile), xlsfile = 'sciadv_adt8600_data_s2.xlsx'; end
if nargin < 3 || isempty(outdir),  outdir  = 'output'; end
P.outdir = outdir;
P.LW = 3.5;  P.FSZ = 26;  P.FIG = [22 19];
P.green = [0.20 0.55 0.32];  P.blue = [0.08 0.24 0.48];
P.grey  = [0.45 0.45 0.45];  P.est  = [0.35 0.10 0.55];
P.B = 4000;
% Mhat is derived from the loaded latencies rather than from a hardcoded 3.93,
% so this panel cannot drift from the data the rest of the paper fits.
P.Mhat = [];   % filled after the load, below

rng(31);
D = read_pacher(csvfile, xlsfile);
assert_key(D);
td = D.tdet;                                    % attacks only, 125 latencies
P.stat  = (1/(mean(1./td) - 1/mean(td)))/mean(td);
P.Mhat  = calib_pool('invert', P.stat);
fprintf('  shape %.4f -> Mhat %.2f (expect 3.928 -> 13.52)\n', P.stat, P.Mhat);
out.A = panelA(D, P);
out.B = panelB(D, P);
out.C = panelC(D, P);
fprintf('\n');
end

%% ============== A: latency against area, clustered ======================
function A = panelA(D, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
sel = D.attack & ~isnan(D.tfirst) & ~isnan(D.area);   % detections only
a = D.area(sel); t = D.tfirst(sel);
plot(ax, a, t, 'o', 'MarkerEdgeColor', P.grey, 'MarkerSize', 7, 'LineWidth', 1.5);
b = [ones(numel(a),1) log(a)] \ log(t);
xg = linspace(min(a), max(a), 200)';
plot(ax, xg, exp(b(1))*xg.^b(2), '-', 'Color', P.green, 'LineWidth', P.LW);

% clustered bootstrap on the log-log slope
[~, ~, ci] = unique(D.clust(sel));
cu = unique(ci); sl = nan(P.B,1);
for k = 1:P.B
    pick = cu(randi(numel(cu), numel(cu), 1));
    idx = [];
    for j = 1:numel(pick), idx = [idx; find(ci == pick(j))]; end %#ok<AGROW>
    if numel(unique(a(idx))) < 5, continue; end
    bb = [ones(numel(idx),1) log(a(idx))] \ log(t(idx));
    sl(k) = bb(2);
end
sl = sl(~isnan(sl)); ciS = prctile(sl, [2.5 97.5]);
set(ax,'XScale','log','YScale','log','TickDir','out', ...
       'TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'shoal area $A$ (m$^2$)','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'first response time $t_{(1)}$ (s)','Interpreter','latex','FontSize',P.FSZ+8);
fprintf(['  panel A: n = %d timed attacks over %d clusters\n' ...
         '           log-log slope %.3f, clustered 95%% CI [%.3f, %.3f]\n'], ...
        numel(a), numel(cu), b(2), ciS(1), ciS(2));
fprintf('           interval crosses zero: %d\n', ciS(1) < 0 && ciS(2) > 0);
% Expect n = 125 over 42 clusters, slope -0.429, CI [-0.841, 0.071]. If n
% comes back 148, the selection has picked up the 23 answered flybys, which
% are false-alarm latencies rather than detections, and the slope drifts to
% -0.382. Shoal area is near-constant within a recording-bout (within-cluster
% SD of log A is 0.018 against 0.234 overall, under 1% of the variance), so
% this slope is identified across recordings, not across events, and the
% width of the interval reflects that rather than thin timing data.
%  POST: state the slope and the CLUSTERED interval, and say in the caption
%        that the interval crosses zero, so the panel checks the sign of the
%        area dependence and does not bound M.
export(f,'figs6A',P);
A = struct('slope',b(2),'ci',ciS,'boot',sl);
end

%% ============== B: area does not track headcount within site ===========
function B = panelB(D, P)
if isempty(D.s2)
    fprintf('  panel B skipped (Data S2 not supplied)\n'); B = []; return
end
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
s2 = D.s2; sites = unique(s2.site);
cols = [P.green; P.blue; P.grey];
if numel(sites) > size(cols,1)       % three sites expected; do not index past it
    cols = repmat(cols, ceil(numel(sites)/size(cols,1)), 1);
    warning('collesc:sites','%d sites found, expected 3', numel(sites));
end
sl = nan(numel(sites),1);
for i = 1:numel(sites)
    m = s2.site == sites(i);
    plot(ax, s2.area(m), s2.N(m), 'o', 'MarkerFaceColor', cols(i,:), ...
         'MarkerEdgeColor', cols(i,:), 'MarkerSize', 8);
    bb = [ones(sum(m),1) log(s2.area(m))] \ log(s2.N(m));
    sl(i) = bb(2);
    xg = linspace(min(s2.area(m)), max(s2.area(m)), 50)';
    plot(ax, xg, exp(bb(1))*xg.^bb(2), '-', 'Color', cols(i,:), 'LineWidth', P.LW);
end
bp = [ones(numel(s2.area),1) log(s2.area)] \ log(s2.N);
xg = linspace(min(s2.area), max(s2.area), 100)';
plot(ax, xg, exp(bp(1))*xg.^bp(2), '--', 'Color',[0.1 0.1 0.1], 'LineWidth', P.LW);
set(ax,'XScale','log','YScale','log','TickDir','out', ...
       'TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'shoal area $A$ (m$^2$)','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'tracked count','Interpreter','latex','FontSize',P.FSZ+8);
fprintf('  panel B: pooled slope %.3f; within-site slopes', bp(2));
fprintf(' %+.2f', sl); fprintf('\n');
%  POST: dashed black is the pooled fit, which looks like constant density.
%        The solid within-site fits are what matters, and they are flat or
%        negative, so the pooled slope is a between-site artifact and shoal
%        area is not a proxy for pool size.
export(f,'figs6B',P);
B = struct('pooled_slope',bp(2),'site_slopes',sl);
end

%% ============== C: density from spacing, not from tracking =============
function C = panelC(D, P)
if isempty(D.s2)
    fprintf('  panel C skipped (Data S2 not supplied)\n'); C = []; return
end
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
s2 = D.s2; nnd = s2.nnd_mm/1000;                    % m
rho_pois = 1./(4*nnd.^2);                           % mean NND = 1/(2 sqrt(rho))
rho_hex  = 2./(sqrt(3)*nnd.^2);                     % close packing
rho_trk  = s2.N ./ s2.area;                         % NOT a density

Ag = linspace(30, 80, 100)';
patch(ax, [Ag; flipud(Ag)], [median(rho_pois)*Ag; flipud(median(rho_hex)*Ag)], ...
      P.green, 'FaceAlpha', 0.22, 'EdgeColor','none');
plot(ax, Ag, median(rho_pois)*Ag, '-', 'Color', P.green, 'LineWidth', P.LW);
plot(ax, Ag, median(rho_hex)*Ag,  '-', 'Color', P.green, 'LineWidth', P.LW);
plot(ax, Ag, median(rho_trk)*Ag,  '--', 'Color', P.grey, 'LineWidth', P.LW);
plot(ax, Ag, P.Mhat*ones(size(Ag)), '-', 'Color', P.est, 'LineWidth', P.LW);
set(ax,'YScale','log','TickDir','out','TickLabelInterpreter','latex', ...
       'FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'shoal area $A$ (m$^2$)','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'number of fish','Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[30 80]); ylim(ax,[5 1e6]);
fprintf('  panel C: NND %.1f mm -> rho = %.0f/m2 (Poisson) to %.0f/m2 (hex)\n', ...
        1000*median(nnd), median(rho_pois), median(rho_hex));
fprintf('           at A = 56.2 m2: %.2e to %.2e fish, against Mhat = %.2f\n', ...
        56.2*median(rho_pois), 56.2*median(rho_hex), P.Mhat);
fprintf('           pooling count is %.1e of the spacing-implied headcount\n', ...
        P.Mhat/(56.2*median(rho_pois)));
fprintf('           tracked-count ratio understates density by %.0fx\n', ...
        median(rho_pois)/median(rho_trk));
fprintf('           median BL %.1f mm -> spacing %.2f body lengths\n', ...
        median(s2.bl_mm), median(nnd)*1000/median(s2.bl_mm));
%  POST: green band is the headcount implied by the measured nearest-neighbor
%        distance, between the Poisson and close-packed readings. Purple is
%        the fitted pooling count. Grey dashed is Ntags_total/area, which is
%        two orders of magnitude below the spacing-implied density and is
%        therefore a tracked subset rather than a census. The pooling count
%        is about 1e-4 of the fish present, which is the quantitative form of
%        the correlated-and-occluded-visual-fields argument.
export(f,'figs6C',P);
C = struct('rho_poisson',median(rho_pois),'rho_hex',median(rho_hex), ...
           'rho_tracked',median(rho_trk));
end

%% ============================ helpers ==================================
function export(f, name, P)
if exist(P.outdir,'dir')
    exportgraphics(f, fullfile(P.outdir,[name '.pdf']), 'ContentType','vector');
    exportgraphics(f, fullfile(P.outdir,[name '.png']), 'Resolution',400);
    fprintf('           wrote %s.pdf / .png\n', name);
else
    fprintf('           outdir "%s" not found -- %s not exported\n', P.outdir, name);
end
end

function assert_key(D)
%ASSERT_KEY  Guard against the column-index error that keyed clusters on
%   location|bout instead of file|bout. That gave 83/47/17 and silently
%   narrowed the flyby stratum, which every q-derived interval below rests
%   on. Correct is 73 clusters, 47 carrying attacks, 26 carrying flybys.
[~, ~, ci] = unique(D.clust);
n  = numel(unique(ci));
nA = numel(unique(ci(D.attack)));
nF = numel(unique(ci(~D.attack)));
fprintf('  cluster key: %d clusters, %d attack, %d flyby\n', n, nA, nF);
if ~isequal([n nA nF], [73 47 26])
    warning('collesc:clusterkey', ...
        ['cluster key gives %d/%d/%d, expected 73/47/26 -- read_pacher is ' ...
         'keying on the wrong column. Every clustered interval below is ' ...
         'suspect; fix the loader before using this output.'], n, nA, nF);
end
end
