function out = make_figS4(csvfile, outdir)
%MAKE_FIGS4  The over-discounting gap does not depend on the attended count.
%
%  This is the panel retired from Fig. 5, which belongs in the supplement:
%  Fig. 5C collapses the two-dimensional statement onto its worst case, and
%  this figure shows the full plane behind it.
%
%  A  the gap alpha-hat(M) - L_inf(K) over the admissible wedge K <= M, with
%     the diagonal, the admissibility ceiling K_max, and the Mhat interval
%  B  joint bootstrap density over (M, alpha) with L_inf(M) across it, and
%     P(alpha-hat > L_inf) printed
%  C  heterogeneity in the attended count: E[L_inf(K)] against the
%     coefficient of variation at fixed mean, showing the gap WIDENS by
%     concavity rather than closing
%
%  Usage:  out = make_figS4('sciadv_adt8600_data_s1.csv', 'output');
%
%  ZPK 2026

if nargin < 1 || isempty(csvfile), csvfile = 'sciadv_adt8600_data_s1.csv'; end
if nargin < 2 || isempty(outdir),  outdir  = 'output'; end
P.outdir = outdir;
P.LW = 3.5;  P.FSZ = 26;  P.FIG = [22 19];
P.est  = [0.35 0.10 0.55];
P.grey = [0.30 0.30 0.30];
P.greens = [0.75 0.90 0.78; 0.55 0.80 0.60; 0.36 0.68 0.44;
            0.20 0.55 0.32; 0.10 0.42 0.24; 0.03 0.28 0.15];
P.B    = 4000;
% SAME seed and replicate count as make_fig5's bootstrap_joint (opts.seed =
% 20260811, opts.nboot = 4000), so the Mhat interval quoted in the main text
% and drawn three times -- Fig. 5C as a capped bar, panel A as a shaded band,
% panel B as the cloud -- comes from one stream rather than three. The seed is
% applied immediately before jointboot, NOT here: rates() below calls
% calib_pool, which draws 2e6 samples per grid point on a cold cache and would
% otherwise leave the two resamples starting from different states.
P.seed = 20260811;
D = read_pacher(csvfile);
assert_key(D);
S = rates(D);
fprintf('\n  TP = %.4f, q = %.4f, stat = %.3f -> Mhat = %.1f, alpha = %.4f\n', ...
        S.TP, S.q, S.stat, S.M, S.alpha);
% bootstrap FIRST: panel A's Mhat band comes from it
rng(P.seed);                     % see the note above -- must be HERE
[Mb, ab] = jointboot(D, P.B);
ok = ~isnan(ab);  Mb = Mb(ok);  ab = ab(ok);
S.Mci = prctile(Mb, [2.5 97.5]);
S.aci = prctile(ab, [2.5 97.5]);
fprintf('  bootstrap: M in [%.1f, %.1f], alpha in [%.3f, %.3f], %d replicates\n', ...
        S.Mci(1), S.Mci(2), S.aci(1), S.aci(2), numel(ab));
out.A = panelA(S, P);
out.B = panelB(Mb, ab, S, P);
out.C = panelC(S, P);
fprintf('\n');
end

%% ================= A: the gap over the (M,K) wedge =====================
function A = panelA(S, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG(1)+2 P.FIG(2)],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
Mg = linspace(2, 40, 400); Kg = linspace(1, 40, 400);
[MM, KK] = meshgrid(Mg, Kg);
% alpha_hat(M, q, TP) rather than calib_pool('alpha', M). The two agree at the
% point estimate -- calib_pool caches exactly these rates -- but alpha_hat is
% elementwise arithmetic and so is guaranteed to vectorise over the 400x400
% grid, and it makes the rates the gap is computed at explicit rather than
% cached. S.q and S.TP are the same rates panel B's dashed benchmark uses.
G = alpha_hat(MM, S.q, S.TP) - linf_ceiling(KK);
G(KK > MM) = NaN;                                  % K <= M by counting
pc = pcolor(ax, Mg, Kg, G); set(pc,'EdgeColor','none','FaceColor','texturemap')
colormap(ax, greenmap()); caxis(ax,[0 1]);
lev = 0.3:0.1:0.9;
contour(ax, Mg, Kg, G, lev, 'LineColor',[1 1 1], 'LineWidth', 1.8);
% Label on a MASKED copy, so labels appear only where the wedge is wide
% enough to hold one. The contours converge in the lower-left corner, where
% K is small and the gap is largest, and labelling all seven levels there
% stacked three of them on top of each other. Lines come from G, labels from
% Gl, so every contour stays drawn while the text thins out.
Gl = G;  Gl(KK < 8 | MM < 12) = NaN;
[c1,h1] = contour(ax, Mg, Kg, Gl, lev, 'LineColor','none');
clabel(c1, h1, lev, 'Color',[1 1 1], 'FontSize', 18, 'LabelSpacing', 1000);
plot(ax, Mg, Mg, '-', 'Color',[0.1 0.1 0.1], 'LineWidth', P.LW);
plot(ax, xlim(ax), S.Kmax*[1 1], '--', 'Color',[0.1 0.1 0.1], 'LineWidth', 2.5);
patch(ax, [S.Mci(1) S.Mci(2) S.Mci(2) S.Mci(1)], [1 1 40 40], ...
      [1 1 1], 'EdgeColor','none', 'FaceAlpha', 0.30);
plot(ax, S.Mci(1)*[1 1], [1 40], ':', 'Color', P.est, 'LineWidth', 2.5);
plot(ax, S.Mci(2)*[1 1], [1 40], ':', 'Color', P.est, 'LineWidth', 2.5);
plot(ax, S.M*[1 1], [1 40], '-', 'Color', P.est, 'LineWidth', P.LW);
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'pooling count $M$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'attended neighbors $K$','Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[2 40]); ylim(ax,[1 40]);
cb = colorbar(ax);
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\hat{\alpha}(M) - L_\infty(K)$'; cb.Label.FontSize = P.FSZ;
cb.TickLabelInterpreter = 'latex';
fprintf('  panel A: positive fraction of the admissible wedge = %.4f\n', ...
        mean(G(~isnan(G)) > 0));
fprintf('           minimum gap on the wedge = %.3f\n', min(G(~isnan(G))));
%  POST labels: black diagonal "$K = M$", dashed "$K_{\max}$", and ONE label
%        reading "$\hat M$, 95\%" covering the solid purple line and both
%        dotted edges together, as in Fig. 5C. Contour labels are drawn from
%        a masked copy, so the lower-left corner is deliberately bare rather
%        than missing. State that the whole admissible wedge is positive.
export(f,'figs4A',P);
A = struct('M',Mg,'K',Kg,'gap',G);
end

%% ================= B: joint bootstrap over (M, alpha) ==================
function B = panelB(Mb, ab, S, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG(1)+2 P.FIG(2)],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
frac  = mean(ab > linf_ceiling(Mb));
% Minimum excess over replicates. The Fig. S4B caption and §pacherfitting both
% quote this as 0.20, and it is NOT the point-estimate gap of 0.376 -- the two
% differ because the excess decreases in M and the bootstrap reaches well past
% Mhat. Printed so neither number is ever quoted from memory.
minex = min(ab - linf_ceiling(Mb));

xe = linspace(log10(2), log10(60), 90);  ye = linspace(0.7, 1.0, 90);
% hist2 snaps out-of-range values into the edge bin ('nearest','extrap'), so
% replicates above M = 60 would pile into a spurious column at the frame edge.
% Exclude them from the IMAGE only; frac, minex and the printed quantiles
% below stay on the full Mb and ab, which is what the caption quotes.
in = Mb <= 60;
H  = smooth2(hist2(log10(Mb(in)), ab(in), xe, ye), 2.2);  H = H/max(H(:));
pc = pcolor(ax, 10.^xe, ye, H'); set(pc,'EdgeColor','none','FaceColor','texturemap')
colormap(ax, purplemap()); caxis(ax,[0 1]);
Mg = logspace(log10(2), log10(60), 300);
plot(ax, Mg, linf_ceiling(Mg), '-',  'Color','w', 'LineWidth', P.LW+1.5);
plot(ax, Mg, linf_ceiling(Mg), '--', 'Color', P.grey, 'LineWidth', P.LW-1);
set(ax,'XScale','log','TickDir','out','TickLabelInterpreter','latex', ...
       'FontSize',P.FSZ,'LineWidth',2);
set(ax,'XTick',[2 5 10 20 40 60],'XTickLabel',{'2','5','10','20','40','60'});
xlim(ax,[2 60]);
% ylim spans the BENCHMARK, not just the bootstrap cloud. At [0.7 1] the
% L_inf(M) curve only entered the frame past M ~ 33, so the panel showed a
% purple blob and a stray dashed segment rather than the gap it exists to
% display.
ylim(ax,[0 1]);
xlabel(ax,'pooling count $M$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'social discounting rate $\alpha$','Interpreter','latex','FontSize',P.FSZ+8);
cb = colorbar(ax);
cb.Label.Interpreter = 'latex';
cb.Label.String = 'bootstrap density (scaled)'; cb.Label.FontSize = P.FSZ-2;
cb.TickLabelInterpreter = 'latex';
qM = prctile(Mb,[2.5 50 97.5]); qa = prctile(ab,[2.5 50 97.5]);
fprintf('  panel B: M = %.1f [%.1f, %.1f], alpha = %.3f [%.3f, %.3f]\n', ...
        qM(2), qM(1), qM(3), qa(2), qa(1), qa(3));
fprintf('           P(alpha > L_inf(M)) = %.4f over %d replicates\n', frac, numel(ab));
fprintf('           minimum excess alpha - L_inf(M) over replicates = %.4f\n', minex);
fprintf('           %d replicate(s) above the M = 60 axis top, excluded from the image\n', ...
        sum(~in));
%  POST: dashed curve is $L_\infty(M)$, the largest value the individually
%        Bayesian benchmark can take at that pooling count. Annotate the
%        fraction of the density lying above it.
export(f,'figs4B',P);
B = struct('M',Mb,'alpha',ab,'frac_above',frac,'min_excess',minex, ...
           'ciM',qM,'cia',qa);
end

%% ================= C: heterogeneity in K ===============================
function C = panelC(S, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
cv = linspace(0, 1.2, 40);  Kbar = [5 7 14];
EL = nan(numel(Kbar), numel(cv));
for i = 1:numel(Kbar)
    for j = 1:numel(cv)
        EL(i,j) = mean_Linf_gamma(Kbar(i), cv(j));
    end
    plot(ax, cv, EL(i,:), '-', 'Color', P.greens(2*i,:), 'LineWidth', P.LW);
end
plot(ax, xlim(ax), S.alpha*[1 1], '-', 'Color', P.est, 'LineWidth', P.LW);
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'coefficient of variation of $K$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'$\langle L_\infty(K)\rangle$','Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[0 1.2]); ylim(ax,[0 1]);
fprintf('  panel C: Kbar = 7, E[L_inf] = %.4f, %.4f, %.4f at CV = 0, 0.6, 1.0\n', ...
        mean_Linf_gamma(7,0), mean_Linf_gamma(7,0.6), mean_Linf_gamma(7,1.0));
%  POST: purple line is $\hat\alpha$. Curves are $\bar K = 5, 7, 14$. The
%        point is that heterogeneity moves the benchmark DOWN, since
%        $L_\infty$ is concave, so the gap widens.
export(f,'figs4C',P);
C = struct('cv',cv,'Kbar',Kbar,'EL',EL);
end

%% ============================ helpers ==================================
function S = rates(D)
isA = D.attack & ~isnan(D.resp);  isF = ~D.attack & ~isnan(D.resp);
S.TP = mean(D.resp(isA)==1);  S.q = mean(D.resp(isF)==1);
t = D.tdet;      % attacks only
S.stat = (1/(mean(1./t) - 1/mean(t)))/mean(t);
S.M = calib_pool('invert', S.stat);
S.alpha = calib_pool('alpha', S.M);
q1 = 1-(1-S.q)^(1/S.M);
S.Kmax = log(1-q1)/log(1-q1^2);
% Mci is filled by the joint bootstrap below rather than hardcoded, so
% panel A's band always matches panel B's density.
S.Mci = [NaN NaN];
end

function [Mb, ab] = jointboot(D, B)
[~, ~, ci] = unique(D.clust);
cA = unique(ci(D.attack)); cF = unique(ci(~D.attack));
Mb = nan(B,1); ab = nan(B,1);
for b = 1:B
    pick = [cA(randi(numel(cA),numel(cA),1)); cF(randi(numel(cF),numel(cF),1))];
    idx = [];
    for j = 1:numel(pick), idx = [idx; find(ci == pick(j))]; end %#ok<AGROW>
    a = D.attack(idx); r = D.resp(idx); t = D.tfirst(idx);
    okA = a & ~isnan(r); okF = ~a & ~isnan(r);
    tt = t(a & ~isnan(t));            % detection latencies only
    if ~any(okA) || ~any(okF) || numel(tt) < 20, continue; end
    TP = mean(r(okA)==1); q = mean(r(okF)==1);
    if TP <= q, continue; end
    s = (1/(mean(1./tt) - 1/mean(tt)))/mean(tt);
    Mb(b) = calib_pool('invert', s);
    ab(b) = alpha_hat(Mb(b), q, TP);
end
end

function E = mean_Linf_gamma(Kbar, cv)
if cv <= 0, E = linf_ceiling(Kbar); return; end
sh = 1/cv^2; sc = Kbar/sh;
k = max(gamrnd_local(sh, sc, 2e5), 1);
E = mean(linf_ceiling(k));
end

function x = gamrnd_local(a, b, n)
% Marsaglia-Tsang, no Statistics Toolbox
if a < 1
    x = gamrnd_local(a+1, b, n) .* rand(n,1).^(1/a); return
end
d = a - 1/3; c = 1/sqrt(9*d); x = nan(n,1); m = 0;
while m < n
    z = randn(n,1); v = (1 + c*z).^3; u = rand(n,1);
    ok = v > 0 & log(u) < 0.5*z.^2 + d - d*v + d*log(v);
    take = min(sum(ok), n-m);
    if take > 0
        vv = d*v(ok); x(m+1:m+take) = vv(1:take); m = m + take;
    end
end
x = b*x;
end


function H = hist2(x, y, xe, ye)
H = zeros(numel(xe), numel(ye));
ix = interp1(xe, 1:numel(xe), x, 'nearest', 'extrap');
iy = interp1(ye, 1:numel(ye), y, 'nearest', 'extrap');
for k = 1:numel(ix), H(ix(k),iy(k)) = H(ix(k),iy(k)) + 1; end
end

function S = smooth2(H, sd)
r = ceil(3*sd); [gx,gy] = meshgrid(-r:r,-r:r);
g = exp(-(gx.^2+gy.^2)/(2*sd^2)); g = g/sum(g(:));
S = conv2(H, g, 'same');
end

function m = greenmap()
m = [linspace(1,0.03,64)', linspace(1,0.28,64)', linspace(1,0.15,64)'];
end

function m = purplemap()
m = [linspace(1,0.35,64)', linspace(1,0.10,64)', linspace(1,0.55,64)'];
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
