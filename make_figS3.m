function out = make_figS3(outdir)
%MAKE_FIGS3  Robustness of the pooling count $\hat M$.
%
%  Three standalone panels, fully self-contained (no data file needed).
%
%  A  the scale-free statistic lambda/Tbar against theta, at several M, with
%     the self-consistent locus theta_1(M) marked. The statistic is NOT
%     threshold-free at the M of interest: measured spread across theta in
%     [1,6] is 54% at M=5, 18% at M=10, 12% at M=14, 8% at M=20, 5% at
%     M=40. This panel motivates the self-consistent calibration rather than
%     asserting independence: the threshold is fixed by the false alarm rate
%     through theta_1(M) = -log(1-(1-q)^(1/M)), so the loop closes and
%     nothing is tuned. The markers show where each M is actually read.
%  B  right-censoring sensitivity. Simulate at a known M, censor the top x%
%     of latencies, and read the implied Mhat back off the calibration. The
%     observed tail quantiles are marked; the data show no censoring
%     signature, so this bounds how much a censoring reading could move M.
%  C  correlated accumulators. Repeat the calibration with equicorrelated
%     streams and plot Mhat against the correlation. This turns the
%     independence caveat into a number.
%
%  NEITHER B NOR C IS DRAWN OUTSIDE ITS VALID RANGE. The moment ratio
%  inverts against calib_pool's grid (M up to 90) with interp1 'extrap', so
%  anything past the top of that grid is extrapolation; the likelihood
%  searches a bounded grid, so a value sitting at the ceiling is the ceiling
%  rather than an estimate. Panel B now blanks both cases, and panel C's
%  correlation sweep stops at c = 0.45, the last point where Mhat stays
%  inside the calibrated range.
%
%  Usage:  out = make_figS3('output');
%
%  ZPK 2026

if nargin < 1 || isempty(outdir), outdir = 'output'; end
P.outdir = outdir;
P.LW = 3.5;  P.FSZ = 26;  P.FIG = [22 19];
P.greens = [0.75 0.90 0.78; 0.55 0.80 0.60; 0.36 0.68 0.44;
            0.20 0.55 0.32; 0.10 0.42 0.24; 0.03 0.28 0.15];
P.blue   = [0.08 0.24 0.48];
P.nsamp    = 2e5;
P.Mmle_max = 300;                               % likelihood search ceiling
P.stat_obs = 3.93;                              % lambda/Tbar in the molly data
P.Mhat     = calib_pool('invert', P.stat_obs);  % 13.53; rates live in calib_pool

rng(4);
out.A = panelA(P);
out.B = panelB(P);
out.C = panelC(P);
fprintf('\n');
end

%% ============== A: the statistic is near theta-free ====================
function A = panelA(P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
th_g = linspace(1, 6, 11);  Mg = [5 10 14 20 40];
S = nan(numel(Mg), numel(th_g));  Ssc = nan(size(Mg));
for i = 1:numel(Mg)
    a = calib_pool('alpha', Mg(i));
    for j = 1:numel(th_g)
        S(i,j) = poolstat(Mg(i), th_g(j), a, P.nsamp);
    end
    plot(ax, th_g, S(i,:), '-', 'Color', P.greens(min(i+1,6),:), 'LineWidth', P.LW);
    % where this M is actually read: the self-consistent threshold
    thsc = calib_pool('theta1', Mg(i));
    Ssc(i) = poolstat(Mg(i), thsc, a, P.nsamp);
    plot(ax, thsc, Ssc(i), 'o', 'MarkerFaceColor', P.greens(min(i+1,6),:), ...
         'MarkerEdgeColor','w', 'MarkerSize', 14, 'LineWidth', 2);
end
plot(ax, arrayfun(@(m) calib_pool('theta1',m), Mg), Ssc, '-', ...
     'Color', [0.35 0.10 0.55], 'LineWidth', P.LW-1);
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'threshold $\theta$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'shape $\lambda/\bar T_{(1)}$','Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[1 6]);
rel = (max(S,[],2)-min(S,[],2))./mean(S,2);
fprintf('  panel A: relative spread over theta ='); fprintf(' %.3f', rel); fprintf('\n');
fprintf('           self-consistent theta_1(M) =');
fprintf(' %.2f', arrayfun(@(m) calib_pool('theta1',m), Mg)); fprintf('\n');
%  POST: label curves M = 5, 10, 14, 20, 40 (light to dark), placing the
%        M = 5 label clear of the purple locus, which it currently collides
%        with. Markers are the self-consistent points theta = theta_1(M),
%        joined by that locus. State the spreads and say plainly that the
%        statistic is threshold-dependent at small M, which is WHY the
%        calibration is run along the locus rather than at a fixed theta.
export(f,'figs3A',P);
A = struct('theta',th_g,'M',Mg,'stat',S,'relspread',rel,'selfcon',Ssc);
end

%% ============== B: right-censoring sensitivity =========================
function B = panelB(P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
[Mgrid, scurve] = calib(P);
% Range stops at 30%. Beyond that the moment ratio leaves calib_pool's grid
% and the panel would be plotting extrapolation, the same thing panel C's
% sweep was truncated to avoid. Both numbers the caption quotes, 22 at 5%
% and 42 at 20%, sit inside this range, and the argument is about the LEFT
% edge in any case, since the molly latencies show no censoring signature.
frac = 0:0.05:0.30;  Mtrue = [10 round(P.Mhat) 20];
Mimp = nan(numel(Mtrue), numel(frac));  Mmle = nan(numel(Mtrue), numel(frac));
for i = 1:numel(Mtrue)
    a = calib_pool('alpha', Mtrue(i));
    % Threshold at the SELF-CONSISTENT locus theta_1(M), not a fixed value.
    % Panel A exists to show lambda/Tbar is threshold-dependent at these M
    % (12% spread over theta in [1,6] at M = 14), so simulating panel B at a
    % hardcoded theta = 2.0 read the shape off a different point of the very
    % curve panel A plots. Panel C already used theta_1(M); this matches it.
    x = poolsample(Mtrue(i), calib_pool('theta1', Mtrue(i)), a, P.nsamp);
    for j = 1:numel(frac)
        keep = x <= quantile_local(x, 1-frac(j));
        Mimp(i,j)  = interp1(scurve, Mgrid, shapeover(x(keep)), 'linear', 'extrap');
        Mmle(i,j)  = fit_pool_mle(x(keep), struct('ngrid',120,'Mmax',P.Mmle_max));
        % Neither estimator is drawn outside the range it is valid on. Past
        % the top of the calibration grid the moment ratio is extrapolating;
        % at the likelihood's search ceiling the value IS the ceiling, not an
        % estimate.
        if Mimp(i,j) > max(Mgrid),           Mimp(i,j) = NaN; end
        if Mmle(i,j) > 0.98*P.Mmle_max,      Mmle(i,j) = NaN; end
    end
    plot(ax, 100*frac, Mimp(i,:), '--', 'Color', P.greens(2*i,:), 'LineWidth', P.LW-1);
    plot(ax, 100*frac, Mmle(i,:), '-',  'Color', P.greens(2*i,:), 'LineWidth', P.LW);
    plot(ax, 0, Mtrue(i), 'o', 'MarkerFaceColor', P.greens(2*i,:), ...
         'MarkerEdgeColor', P.greens(2*i,:), 'MarkerSize', 10);
end
plot(ax, xlim(ax), P.Mhat*[1 1], ':', 'Color',[0.4 0.4 0.4], 'LineWidth', 2);
set(ax,'YScale','log','TickDir','out','TickLabelInterpreter','latex', ...
       'FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'latencies censored (\%)','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'implied $\hat M$','Interpreter','latex','FontSize',P.FSZ+8);
fprintf('  panel B: at  5%% censoring, true M=%d reads as %.1f (moment ratio) / %.1f (MLE)\n', ...
        Mtrue(2), interp1(frac, Mimp(2,:), 0.05), interp1(frac, Mmle(2,:), 0.05));
fprintf('           at 20%% censoring, true M=%d reads as %.1f (moment ratio) / %.1f (MLE)\n', ...
        Mtrue(2), interp1(frac, Mimp(2,:), 0.20), interp1(frac, Mmle(2,:), 0.20));
nb = sum(isnan(Mimp(:))) + sum(isnan(Mmle(:)));
fprintf('           %d point(s) blanked as out of range (grid top %.0f, MLE ceiling %d)\n', ...
        nb, max(Mgrid), P.Mmle_max);
%  POST: DASHED is the moment ratio (the headline estimator), SOLID the
%        scale-free likelihood cross-check; dotted line is the OBSERVED
%        Mhat = 13.5, which is a different object from panel C's dotted
%        line (the simulation's true M = 14) and should be labelled so.
%        Both estimators are sensitive to censoring, since both must read
%        the shape of a distribution whose scale is not identified. Label
%        the three curves inline in their own colours, as panels A and C do,
%        leaving the legend to carry only the solid/dashed key. Say in the
%        caption that the molly latencies show NO censoring signature (95th
%        pct 9.1 s, max 27.8 s, one event within 2 s of the max), so the
%        relevant point is the left edge; the curves bound how far a
%        censoring reading could move the estimate, not how far it does.
export(f,'figs3B',P);
B = struct('frac',frac,'Mtrue',Mtrue,'Mimplied',Mimp,'Mmle',Mmle);
end

%% ============== C: correlated accumulators =============================
function C = panelC(P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
[Mgrid, scurve] = calib(P);
% Sweep stops at 0.45, the last point at which the recovered Mhat is still
% inside the calibration grid. Past it Mhat falls below M = 2 by
% extrapolation, and a pooling count under two is meaningless.
c_g = 0:0.05:0.45;  Mtrue = round(P.Mhat);  a = calib_pool('alpha', Mtrue);
Mc = nan(size(c_g));
for j = 1:numel(c_g)
    x = poolsample_corr(Mtrue, calib_pool('theta1',Mtrue), a, P.nsamp, c_g(j));
    Mc(j) = interp1(scurve, Mgrid, shapeover(x), 'linear', 'extrap');
end
plot(ax, c_g, Mc, '-', 'Color', P.greens(4,:), 'LineWidth', P.LW);
plot(ax, c_g, Mtrue*ones(size(c_g)), ':', 'Color',[0.4 0.4 0.4], 'LineWidth', 2);
yyaxis(ax,'right');
plot(ax, c_g, arrayfun(@(m) calib_pool('alpha',m), Mc), '-', ...
     'Color', P.blue, 'LineWidth', P.LW);
plot(ax, c_g, arrayfun(@(m) linf_ceiling(m), Mc), '--', 'Color', P.blue, 'LineWidth', 2.5);
ylabel(ax,'$\hat\alpha$ and $L_\infty(\hat M)$','Interpreter','latex','FontSize',P.FSZ+8);
set(ax,'YColor',P.blue); ylim(ax,[0 1]);
yyaxis(ax,'left'); ylabel(ax,'recovered $\hat M$','Interpreter','latex','FontSize',P.FSZ+8);
set(ax,'YColor',P.greens(4,:));
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'accumulator correlation $c$','Interpreter','latex','FontSize',P.FSZ+8);
ga = arrayfun(@(m) calib_pool('alpha',m) - linf_ceiling(m), Mc);
fprintf('  panel C: c = 0 -> Mhat %.1f;  c = %.2f -> Mhat %.1f\n', Mc(1), c_g(end), Mc(end));
fprintf('           gap at c = 0: %+.3f;  at c = %.2f: %+.3f\n', ga(1), c_g(end), ga(end));
fprintf('           gap over the sweep: min %+.3f at c = %.2f, max %+.3f at c = %.2f\n', ...
        min(ga), c_g(find(ga == min(ga),1)), max(ga), c_g(find(ga == max(ga),1)));
fprintf('           gap stays positive: %d\n', all(ga > 0));
%  POST: correlation moves Mhat DOWN sharply -- at c = 0.2 the estimate is
%        already halved -- and alpha-hat with it. The gap against
%        L_inf(Mhat) is NOT monotone: it widens as Mhat falls through about
%        seven, peaking at 0.44 near c = 0.25, then turns over, and over the
%        whole calibrated range it varies only between 0.33 and 0.44. That
%        is the honest statement, not "narrows without closing". The dotted
%        line here is the SIMULATION's true M = 14, not panel B's observed
%        Mhat = 13.5; label the green curve Mhat, since it is recovered.
%        Report the point estimate as conditional on independence and the
%        SIGN as robust over the range tested.
export(f,'figs3C',P);
C = struct('c',c_g,'Mhat',Mc,'gap',ga);
end

%% ============================ helpers ==================================
function [Mg, sg] = calib(~)
[Mg, sg] = calib_pool('grid');
end

function s = poolstat(M, th, a, n)
s = shapeover(poolsample(M, th, a, n));
end

function x = poolsample(M, th, a, n)
M = round(M);
y = igrnd(th/(1-a), th^2/2, n);
k = floor(n/M);
x = min(reshape(y(1:k*M), k, M), [], 2);
end

function x = poolsample_corr(M, th, a, n, c, dt, Tmax)
%POOLSAMPLE_CORR  Minimum of M EQUICORRELATED first passages.
%   Each unit accumulates sqrt(c)*(shared Brownian) + sqrt(1-c)*(own
%   Brownian) at drift 1-a against threshold th, so c is the correlation
%   between accumulators WITHIN a trial. Simulated directly, since there is
%   no closed form for the minimum of correlated first passages.
if nargin < 6 || isempty(dt),   dt   = 0.01; end
if nargin < 7 || isempty(Tmax), Tmax = 80;   end
M   = round(M);  ntr = max(round(n/M), 2000);
mu  = 1 - a;  sig = sqrt(2);  sc = sig*sqrt(dt);
X   = zeros(ntr, M);  out = nan(ntr,1);  t = 0;
for k = 1:round(Tmax/dt)
    t = t + dt;
    zc = randn(ntr,1);  zi = randn(ntr,M);
    X  = X + mu*dt + sc*(sqrt(c)*zc + sqrt(1-c)*zi);
    hit = any(X >= th, 2) & isnan(out);
    if any(hit), out(hit) = t; X(hit,:) = -Inf; end
    if ~any(isnan(out)), break; end
end
x = out(~isnan(out));
end

function s = shapeover(t)
m = mean(t);  s = (1/(mean(1./t) - 1/m)) / m;
end



function q = quantile_local(x, p)
x = sort(x(:)); n = numel(x);
h = (n-1)*p(:) + 1; lo = floor(h); hi = ceil(h);
q = x(lo) + (h - lo).*(x(hi) - x(lo));
end

function export(f, name, P)
if exist(P.outdir,'dir')
    exportgraphics(f, fullfile(P.outdir,[name '.pdf']), 'ContentType','vector');
    exportgraphics(f, fullfile(P.outdir,[name '.png']), 'Resolution',400);
    fprintf('           wrote %s.pdf / .png\n', name);
else
    fprintf('           outdir "%s" not found -- %s not exported\n', P.outdir, name);
end
end
