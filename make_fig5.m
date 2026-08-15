function out = make_fig5(csvfile, opts)
%MAKE_FIG5  Figure 5 for "Social Discounting Enables Fast and Reliable
%           Collective Escape" -- three standalone panels off Pacher Data S1.
%
%   A  response rates vs shoal area (attacks vs flybys)   -- where the rates come from
%   B  first-response-time distribution + M calibration   -- where the pooling count comes from
%   C  alpha-hat(M) against the admissible optimum        -- the verdict
%
%   Usage:  out = make_fig5('DataS1.csv');
%           out = make_fig5('DataS1.csv', 'output');        % directory
%           out = make_fig5('DataS1.csv', struct('nboot',1000));
%
%   Reproduction targets (Aug 2026, self-consistent calibration):
%       177 attacks / 127 responses     TP    = 0.7175
%        81 flybys  /  26 responses     q_obs = 0.3210
%       125 timings                     mean  = 4.92 s
%       lam/mean = 3.93  ->  Mhat = 13.53     alpha-hat = 0.949
%       q_1 = 0.0282     theta_1 = 3.568      K_max = 35.93
%       L_inf(Mhat) = 0.573
%       73 clusters, 47 attack / 26 flyby     4000 replicates, pOver = 1
%       Mci = [8.825, 29.403]                 alphaci = [0.906, 0.982]
%
%   THREE DIFFERENT GAPS, all correct, do not interchange them:
%       0.376  alpha-hat - L_inf(Mhat) at the POINT ESTIMATE (panel C)
%       0.204  min over the 4000 bootstrap replicates of alpha_b - L_inf(M_b)
%       0.26   min over the admissible (M,K) wedge, attained on K = M
%              (that one lives in make_figS4, not here)
%   The gap DECREASES in M along the identification curve (0.417 at M = 8.8
%   to 0.292 at M = 29.4), so the binding replicates are the large-M ones.
%
%   ESTIMATOR is the dimensionless moment ratio lambda/Tbar via calib_pool.
%   The model carries no timescale of its own, so only scale-free features
%   of the latencies identify M; a likelihood on raw seconds is not
%   scale-invariant and returns 4.3 here by absorbing the units mismatch.
%   With a free scale profiled out the likelihood is scale-invariant but has
%   no advantage over the ratio (recovery [10.1, 22.6] against [10.3, 22.5]
%   at n = 125). fit_pool_mle.m runs that version as a cross-check.
%
%   The calibration lives in calib_pool.m and is run at theta_1(M) at each
%   grid point rather than at a fixed threshold, because lambda/Tbar is NOT
%   threshold-free at this M (11% spread over theta in [1,6] at M = 14).

if nargin < 1 || isempty(csvfile), csvfile = find_datas1(); end
if nargin < 2, opts = struct; end
opts = defaults(opts);

D = load_pacher(csvfile);
S = summarise(D);
B = bootstrap_joint(D, opts);

out = struct('data',D,'stats',S,'boot',B);

panelA(D, S, opts);
panelB(D, S, B, opts);
panelC(S, B, opts);

end

% =========================================================================
% data
% =========================================================================
function f = find_datas1()
%FIND_DATAS1  Locate the Pacher Data S1 CSV when none is passed in.
pat = {'DataS1.csv','Data_S1.csv','*data_s1*.csv','*Data*S1*.csv','*pacher*.csv','*adt8600*.csv'};
for i = 1:numel(pat)
    d = dir(pat{i});
    d = d(~[d.isdir]);
    if ~isempty(d), f = fullfile(d(1).folder, d(1).name); return; end
end
error('make_fig5:nofile', ...
    ['No Data S1 CSV found in %s.\n' ...
     'Call it with an explicit path, e.g.  make_fig5(''~/Dropbox/escape/data/DataS1.csv'')'], pwd);
end

function D = load_pacher(csvfile)
%LOAD_PACHER  Parse Data S1. Semicolon-delimited, 19 lines of legend, header
%   on line 20, 17 columns. Line 2 claims the sheet is attacks only -- it is
%   not; the 81 flybys sit in the same block with risk = 'flyby'. Below the
%   real data is a second stacked block (26 blanks, a repeated header, and 47
%   rows shifted one column left) which must be dropped: filtering rows whose
%   risk field is exactly 'predator_attack' or 'flyby' removes it.
%
%   Columns used: 1 file, 10 risk, 11 area (sqm), 12 true_positive,
%                 13 t_first (FRAMES at 25 fps), 14 bout_id.
%   Missing values are the literal string 'NA'.

HDRLINE = 20; FPS = 25;
fid = fopen(csvfile,'r');
if fid < 0, error('make_fig5:csv','cannot open %s', csvfile); end
raw = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
raw = raw{1};
rows = raw(HDRLINE+1:end);

nm = {}; risk = {}; area = []; tp = []; tf = []; bout = {};
for i = 1:numel(rows)
    f = strsplit(rows{i}, ';', 'CollapseDelimiters', false);
    if numel(f) < 14, continue; end
    r = strtrim(f{10});
    if ~(strcmp(r,'predator_attack') || strcmp(r,'flyby')), continue; end
    nm{end+1,1}   = strtrim(f{1});                              %#ok<AGROW>
    risk{end+1,1} = r;                                          %#ok<AGROW>
    area(end+1,1) = str2num_na(f{11});                          %#ok<AGROW>
    tp(end+1,1)   = str2num_na(f{12});                          %#ok<AGROW>
    tf(end+1,1)   = str2num_na(f{13}) / FPS;                    %#ok<AGROW>
    bout{end+1,1} = strtrim(f{14});                             %#ok<AGROW>
end

D.file   = nm;
D.attack = strcmp(risk,'predator_attack');
D.area   = area;
D.resp   = tp;                     % 1 responded, 0 not, NaN unscored
D.tfirst = tf;                     % seconds, NaN where NA -- ALL events
D.tdet   = tf(D.attack & ~isnan(tf));  % detection latencies, attacks only
D.clust  = clusterkey(nm, bout);   % see below

% D.tfirst carries timings for every scored event, including 23 flybys that
% were answered. Those are false-alarm response latencies, not detections,
% and must not enter the pooling-count fit -- including them gives 148
% latencies rather than 125. Use D.tdet for anything modeling
% E[T_(1) | H = 1].

% bout_id restarts within each recording (only 10 distinct values across 18
% recordings), so the key must be file|bout rather than bout alone. That
% gives 73 clusters, 47 carrying attacks and 26 carrying flybys. The NA
% branch in clusterkey never fires on this dataset -- flybys DO carry bout
% labels ('kk' 49 times, 'other' 32, otherwise '1' to '9').
end

function key = clusterkey(file, bout)
n = numel(file); key = cell(n,1); k = 0;
for i = 1:n
    if isempty(bout{i}) || strcmpi(bout{i},'NA')
        k = k + 1; key{i} = sprintf('%s|solo%04d', file{i}, k);
    else
        key{i} = sprintf('%s|%s', file{i}, bout{i});
    end
end
end

function v = str2num_na(s)
s = strtrim(s);
if isempty(s) || strcmpi(s,'NA'), v = NaN; else, v = str2double(s); end
end

% =========================================================================
% statistics
% =========================================================================
function S = summarise(D)
isA = D.attack & ~isnan(D.resp);
isF = ~D.attack & ~isnan(D.resp);

S.nA = sum(isA);  S.kA = sum(D.resp(isA)==1);
S.nF = sum(isF);  S.kF = sum(D.resp(isF)==1);
S.TP = S.kA / S.nA;
S.q  = S.kF / S.nF;

t = D.tdet;      % attacks only
S.t = t;  S.nT = numel(t);
S.tmean = mean(t);
S.lam   = ig_shape(t);            % inverse-Gaussian MLE shape
S.stat  = S.lam / S.tmean;        % scale-free -- this is what identifies M

S.M = invert_M(S.stat);
S.alpha  = alphahat_from_rates(S.TP, S.q, S.M);
S.q1     = 1 - (1 - S.q)^(1/S.M);
S.theta1 = -log(S.q1);
S.Kmax   = log(1 - S.q1) / log(1 - S.q1^2);
end

function lam = ig_shape(t)
% MLE for the inverse-Gaussian shape at fixed mean: 1/lam = mean(1/t - 1/mu)
mu  = mean(t);
lam = 1 / mean(1./t - 1/mu);
end

function B = bootstrap_joint(D, opts)
%BOOTSTRAP_JOINT  Cluster bootstrap, stratified by event type so the thin
%   flyby stratum (26 clusters) cannot vanish. Propagates BOTH rate sampling
%   error and timing/M uncertainty -- alpha depends on M, so resampling them
%   independently would understate the interval.
rng(opts.seed);
[keys, ~, ci] = unique(D.clust);
isA = D.attack;
cA  = unique(ci(isA));  cF = unique(ci(~isA));
nrep = opts.nboot;
B.M = nan(nrep,1); B.alpha = nan(nrep,1); B.TP = nan(nrep,1); B.q = nan(nrep,1);

for b = 1:nrep
    pick = [cA(randi(numel(cA), numel(cA), 1)); ...
            cF(randi(numel(cF), numel(cF), 1))];
    idx = [];
    for j = 1:numel(pick), idx = [idx; find(ci == pick(j))]; end %#ok<AGROW>

    a  = D.attack(idx);  r = D.resp(idx);  t = D.tfirst(idx);
    okA = a & ~isnan(r); okF = ~a & ~isnan(r);
    if ~any(okA) || ~any(okF), continue; end
    TP = mean(r(okA)==1);  q = mean(r(okF)==1);
    tt = t(a & ~isnan(t));            % detection latencies only
    if numel(tt) < 20 || TP <= q, continue; end

    M = invert_M(ig_shape(tt)/mean(tt));
    B.TP(b) = TP; B.q(b) = q; B.M(b) = M;
    B.alpha(b) = alphahat_from_rates(TP, q, M);
end

ok = ~isnan(B.alpha);
B.Mci     = prctile(B.M(ok),     [2.5 97.5]);
B.alphaci = prctile(B.alpha(ok), [2.5 97.5]);
B.pOver   = mean(B.alpha(ok) > linf_ceiling(B.M(ok)));
B.nvalid  = sum(ok);
% Minimum excess over replicates. This is the number the paper reports as
% "a minimum excess of 0.20", NOT the point-estimate gap of 0.376, which is
% what panel C shades at Mhat. The two differ because the gap decreases in M
% and the bootstrap reaches well past Mhat.
B.minexcess = min(B.alpha(ok) - linf_ceiling(B.M(ok)));
fprintf(['  bootstrap: M in [%.2f, %.2f], alpha in [%.3f, %.3f]\n' ...
         '             %d/%d valid, P(alpha > L_inf) = %.4f, min excess = %.4f\n'], ...
        B.Mci(1), B.Mci(2), B.alphaci(1), B.alphaci(2), ...
        B.nvalid, nrep, B.pOver, B.minexcess);
end

% =========================================================================
% panels
% =========================================================================
function panelA(D, S, opts)
%PANEL A  response rate vs shoal area. This is where TP and q come from, and
%   it carries the point that the false-alarm rate is measured on flybys.
[f, ax] = newfig(opts, 'A');
plotbinned(ax, D, true,  opts.colAttack, 'o', opts);
plotbinned(ax, D, false, opts.colFlyby,  's', opts);
xg = linspace(min(D.area(~isnan(D.area))), max(D.area(~isnan(D.area))), 200)';
plot(ax, xg, logistic_curve(D, true,  xg), '-', 'Color', opts.colAttack, 'LineWidth', opts.LW);
plot(ax, xg, logistic_curve(D, false, xg), '-', 'Color', opts.colFlyby,  'LineWidth', opts.LW);
xlim(ax, [min(xg) max(xg)]); ylim(ax, [0 1]);
lab(ax, 'shoal area $A$ (m$^2$)', 'response probability', opts);
txt(ax, 0.58, 0.90, 'true positives (attacks)', opts.colAttack, opts);
txt(ax, 0.58, 0.18, 'false positives (flybys)', opts.colFlyby, opts);
finish(f, opts, 'fig5A');
end

function panelB(D, S, B, opts) %#ok<INUSL>
%PANEL B  first-response-time distribution, with the calibration that turns
%   it into M shown as an inset.
[f, ax] = newfig(opts, 'B');
edges = 0:1:ceil(max(S.t));
histogram(ax, S.t, edges, 'Normalization','pdf', ...
    'FaceColor', opts.colHist, 'EdgeColor','w', 'FaceAlpha', 1);
tt = linspace(0.2, max(S.t), 400);
plot(ax, tt, igpdf(tt, S.tmean, S.lam), '-', 'Color', opts.colAttack, 'LineWidth', opts.LW);
xlim(ax, [0 min(max(S.t), 20)]); ylim(ax, [0 1.15*max(igpdf(tt, S.tmean, S.lam))]);
lab(ax, 'first response time $t_{(1)}$ (s)', 'probability density', opts);
txt(ax, 0.62, 0.62, 'group minimum', opts.colAttack, opts, -2);

p  = get(ax, 'Position');
ai = axes('Parent', f, 'Position', [p(1)+0.52*p(3), p(2)+0.52*p(4), 0.40*p(3), 0.38*p(4)]);
hold(ai, 'on'); box(ai, 'on');
[Mg, sg] = calibration_grid();
plot(ai, Mg, sg, '-o', 'Color', [0.2 0.2 0.2], 'MarkerFaceColor','w', ...
     'LineWidth', opts.LW-1.5, 'MarkerSize', 6);
plot(ai, [Mg(1) S.M S.M], [S.stat S.stat 0], ':', 'Color', opts.colAttack, 'LineWidth', 2.5);
set(ai, 'XScale','log', 'XLim',[2 90], 'YLim',[0 1.15*max(sg)], ...
    'XTick',[2 5 10 20 50 90], 'FontSize', opts.fs-10, 'TickDir','out', ...
    'TickLabelInterpreter','latex', 'LineWidth', 1.5);
xlabel(ai, 'pooling count $M$', 'Interpreter','latex', 'FontSize', opts.fs-8);
ylabel(ai, 'shape $\lambda/\bar t$', 'Interpreter','latex', 'FontSize', opts.fs-8);
txt(ai, 0.52, 0.20, sprintf('$\\hat M = %.1f$', S.M), opts.colAttack, opts, -8);
finish(f, opts, 'fig5B');
end

function panelC(S, B, opts)
%PANEL C  the verdict. x is the pooling count M throughout -- alpha-hat is a
%   function of M alone, K enters only the benchmark. Since K <= M and L_inf
%   is increasing, the optimum at abscissa M lies anywhere in [0, L_inf(M)];
%   shading that region and putting alpha-hat above it says what the (M,K)
%   heatmap said, without a second axis and without the K/M conflation.
[f, ax] = newfig(opts, 'C');
Mg = linspace(2, opts.Mmax, 400)';
Lg = linf_ceiling(Mg);
Ag = alphahat_from_rates(S.TP, S.q, Mg);

% 95% interval on Mhat, drawn as a capped bar near the top rather than as a
% full-height shaded band. The band read as an unexplained region and had to
% be described in the caption; a bar with caps is self-evidently an interval.
yb = 1.035;  cap = 0.022;
plot(ax, B.Mci, [yb yb], '-', 'Color', opts.colEst, 'LineWidth', 3);
plot(ax, [B.Mci(1) B.Mci(1)], yb+[-cap cap], '-', 'Color', opts.colEst, 'LineWidth', 3);
plot(ax, [B.Mci(2) B.Mci(2)], yb+[-cap cap], '-', 'Color', opts.colEst, 'LineWidth', 3);
plot(ax, S.M, yb, 'o', 'MarkerFaceColor', opts.colEst, 'MarkerEdgeColor','w', ...
     'MarkerSize', 9, 'LineWidth', 1.5);
fillbetween(ax, Mg, Lg, Ag, opts.colGap, 0.22);
fillbetween(ax, Mg, zeros(size(Mg)), Lg, opts.colOpt, 0.16);
plot(ax, [2 opts.Mmax], [1 1], '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 2);
plot(ax, [2 opts.Mmax], [0 0], ':', 'Color', [0.65 0.65 0.65], 'LineWidth', 2);
plot(ax, Mg, Lg, '-',  'Color', opts.colOpt, 'LineWidth', opts.LW);
plot(ax, Mg, Ag, '-',  'Color', opts.colEst, 'LineWidth', opts.LW);
errorbar(ax, S.M, S.alpha, S.alpha - B.alphaci(1), B.alphaci(2) - S.alpha, ...
    'o', 'Color', opts.colEst, 'MarkerFaceColor', opts.colEst, ...
    'MarkerSize', 11, 'LineWidth', 2.5, 'CapSize', 10);
xlim(ax, [2 opts.Mmax]); ylim(ax, [-0.03 1.10]);
lab(ax, 'pooling count $M$', 'social discounting rate $\alpha$', opts);

txt(ax, 0.80, 0.955, '$\alpha = 1$',        [0.55 0.55 0.55], opts, -6);
txt(ax, 0.24, 0.045, 'naive $\alpha = 0$',  [0.55 0.55 0.55], opts, -6);
txt(ax, 0.79, 0.855, '$\hat\alpha(M)$',     opts.colEst,      opts);
txt(ax, 0.58, 0.66,  'excess discounting',  [0.35 0.35 0.35], opts, -2);
txt(ax, 0.63, 0.32,  'individually Bayesian', opts.colOpt,    opts, -2);
txt(ax, 0.63, 0.24,  '$L_\infty(K),\ K \le M$', opts.colOpt,  opts, -4);
fM = (B.Mci(2) + 1.5 - 2)/(opts.Mmax - 2);
txt(ax, fM, (1.035 + 0.03)/1.13, '$\hat M$, 95\%', opts.colEst, opts, -6);
fa = (S.M + 2.2 - 2)/(opts.Mmax - 2);
txt(ax, fa, (S.alpha - 0.09 + 0.03)/1.13, sprintf('$\\hat\\alpha = %.2f$', S.alpha), ...
    opts.colEst, opts, -2);
finish(f, opts, 'fig5C');
end

% =========================================================================
% model
% =========================================================================

function a = alphahat_from_rates(TP, q, M)
%ALPHAHAT_FROM_RATES  Eq. (13) inversion at pooling count M.
%   Per-unit rates from the shoal-level ones, then the Mobius form.
%   Reproduces 0.8836 at M = 7 and 0.9511 at M = 14. Kept explicit here
%   because the bootstrap calls it with RESAMPLED rates, which calib_pool's
%   cached rates would not reflect.
q1 = 1 - (1 - q).^(1./M);
m1 = (1 - TP).^(1./M);
r  = log(q1) ./ log(m1);
a  = (r - 1) ./ (r + 1);
end

function M = invert_M(stat)
%INVERT_M  map lam/mean onto the pooling count. Delegates to calib_pool so
%   every figure inverts against ONE curve, calibrated self-consistently at
%   theta_1(M) rather than at an arbitrary threshold. See calib_pool.m for
%   why the threshold is not free.
M = calib_pool('invert', stat);
end

function [Mg, sg] = calibration_grid()
[Mg, sg] = calib_pool('grid');
end


% =========================================================================
% helpers
% =========================================================================
function y = igpdf(t, mu, lam)
y = sqrt(lam./(2*pi*t.^3)) .* exp(-lam*(t-mu).^2 ./ (2*mu^2*t));
end


function plotbinned(ax, D, wantAttack, col, mk, opts)
sel = (D.attack == wantAttack) & ~isnan(D.resp) & ~isnan(D.area);
a = D.area(sel); r = D.resp(sel);
q = quantile_local(a, [0 .2 .4 .6 .8 1]);
for j = 1:5
    in = a >= q(j) & a <= q(j+1);
    if sum(in) < 3, continue; end
    k = sum(r(in)==1); n = sum(in);
    [lo, hi] = wilson(k, n);
    x = median(a(in));
    plot(ax, [x x], [lo hi], '-', 'Color', col, 'LineWidth', 2.0);
    plot(ax, x, k/n, mk, 'MarkerFaceColor', col, 'MarkerEdgeColor', col, ...
         'MarkerSize', 4 + 10*n/numel(a)*5, 'LineWidth', 2.0);
end
end

function y = logistic_curve(D, wantAttack, xg)
sel = (D.attack == wantAttack) & ~isnan(D.resp) & ~isnan(D.area);
X = [ones(sum(sel),1) D.area(sel)]; yv = double(D.resp(sel)==1);
b = irls_logit(X, yv);
y = 1 ./ (1 + exp(-[ones(numel(xg),1) xg] * b));
end

function b = irls_logit(X, y)
% plain Newton-Raphson logistic fit -- no Statistics Toolbox
b = zeros(size(X,2),1);
for it = 1:60
    p = 1./(1+exp(-X*b));
    W = p.*(1-p) + 1e-9;
    z = X*b + (y - p)./W;
    bn = (X' * (W .* X)) \ (X' * (W .* z));
    if norm(bn - b) < 1e-10, b = bn; break; end
    b = bn;
end
end

function q = quantile_local(x, p)
x = sort(x(:)); n = numel(x);
h = (n-1)*p(:) + 1; lo = floor(h); hi = ceil(h);
q = x(lo) + (h - lo).*(x(hi) - x(lo));
end

function fillbetween(ax, x, ylo, yhi, col, alpha)
patch(ax, [x; flipud(x)], [ylo; flipud(yhi)], col, 'EdgeColor','none', 'FaceAlpha', alpha);
end

function [f, ax] = newfig(opts, letter)
%NEWFIG  house style, as make_fig4: 22 x 19 cm, LaTeX throughout, ticks out,
%   no box, LW 3.5 on data, 2 on axes. Panel letters are added in assembly,
%   so opts.letters is false by default.
f  = figure('Units','centimeters','Position',[1 1 opts.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
set(ax,'TickDir','out','TickLabelInterpreter','latex', ...
       'FontSize',opts.fs,'LineWidth',opts.axLW,'Layer','top');
if ~isempty(letter) && opts.letters
    annotation(f,'textbox',[0.005 0.90 0.10 0.09],'String',letter, ...
        'FontWeight','bold','FontSize',opts.fs+12,'EdgeColor','none', ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
end
end

function lab(ax, xs, ys, opts)
xlabel(ax, xs, 'Interpreter','latex', 'FontSize', opts.fs+8);
ylabel(ax, ys, 'Interpreter','latex', 'FontSize', opts.fs+8);
end

function txt(ax, fx, fy, str, col, opts, dfs)
%TXT  place a label at fractional position (fx,fy) of the CURRENT axis
%   limits, converted to data units. Setting 'Units','normalized' in the
%   text() call itself is unreliable: Position is applied in data units
%   first and only reinterpreted afterwards. Requires xlim/ylim to be set
%   before the call.
if nargin < 7, dfs = 0; end
xl = get(ax,'XLim'); yl = get(ax,'YLim');
if strcmp(get(ax,'XScale'),'log')
    x = 10^(log10(xl(1)) + fx*(log10(xl(2))-log10(xl(1))));
else
    x = xl(1) + fx*diff(xl);
end
y = yl(1) + fy*diff(yl);
text(ax, x, y, str, 'Color', col, 'FontSize', opts.fs + dfs, ...
     'Interpreter','latex', 'HorizontalAlignment','center', ...
     'VerticalAlignment','middle', 'Clipping','off');
end

function finish(f, opts, name)
%FINISH  export both formats into opts.outdir, as make_fig4's export().
if exist(opts.outdir,'dir')
    exportgraphics(f, fullfile(opts.outdir,[name '.pdf']), 'ContentType','vector');
    exportgraphics(f, fullfile(opts.outdir,[name '.png']), 'Resolution',400);
    fprintf('           wrote %s.pdf / .png\n', name);
else
    fprintf('           outdir "%s" not found -- %s not exported\n', opts.outdir, name);
end
end

function opts = defaults(opts)
% Second argument is an options STRUCT here but an output DIRECTORY string in
% make_figS3--S6. Accept either, so make_fig5(csv,'output') does the obvious
% thing rather than erroring in fieldnames().
if ischar(opts) || isstring(opts), opts = struct('outdir', char(opts)); end
d = struct( ...
    'FIG', [22 19], ...                         % centimetres, as make_fig4
    'fs', 26, 'LW', 3.5, 'axLW', 2.0, ...       % P.FSZ, P.LW
    'letters', false, ...                       % panel letters added in assembly
    'nboot', 4000, 'seed', 20260811, ...
    'Mmax', 40, ...
    'colAttack', [0.20 0.55 0.32], ...          % P.greens(4,:)
    'colFlyby',  [0.08 0.24 0.48], ...          % P.blues(3,:)
    'colHist',   [0.75 0.90 0.78], ...          % P.greens(1,:)
    'colOpt',    [0.30 0.30 0.30], ...   % theory: achromatic
    'colEst',    [0.35 0.10 0.55], ...
    'colGap',    [0.35 0.10 0.55], ...   % excess, tinted to match colEst
    'outdir', 'output');
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = d.(fn{i}); end
end
end
