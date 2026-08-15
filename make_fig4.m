function out = make_fig4(bayes, outdir)
%MAKE_FIG4  Fig. 4 as FOUR STANDALONE FIGURES, one per panel.
%
%  Writes fig4A, fig4B, fig4C, fig4D as .pdf and .png,
%  each sized and typeset for standalone
%  use, to be assembled externally.  Nothing is tiled.
%
%  A  total survival drift (N-1)lambda^(N)(t) vs rescaled time, inset of the
%     extrapolated tails against the closed form L_inf(N)
%     -- REQUIRES the solver output from solve_bayes
%  B  group first-departure density, two stacked tiles in one figure:
%     naive over Bayesian at the self-consistent alpha = L_inf(N)
%  C  group speed-accuracy Pareto frontier
%  D  discounting rate REQUIRED to hold the false alarm rate, over (q_1, K)
%
%  Panels B, C, D are analytic and self-contained.  Panel A takes
%      bayes(i).N  .t  .lam  .thetaN       (lam is PER-NEIGHBOUR lambda^(N))
%  and is skipped if bayes = [].
%
%  THRESHOLD (K = 1 anchor).  theta_1 = -log(q_1) is the SOLITARY threshold,
%  so the solitary false alarm rate is exp(-theta_1) with no alpha in it.
%  Holding the group rate there gives
%      theta_N = (theta_1 + log N)/(1 + alpha),
%  not theta_1 + log(N)/(1+alpha), which pinned to exp(-(1+alpha)*theta_1),
%  the rate of an agent discounting neighbours it does not have.
%
%  q_1 IS PER-UNIT, NOT THE GROUP RATE.  The observed q_N = 0.321 is the
%  SHOAL-level false alarm rate, pooled over the M responders whose earliest
%  crossing registers as a group response. The solitary rate is
%      q_1 = 1 - (1-q_N)^(1/M) = 0.0282   at Mhat = 13.52,
%  giving theta_1 = 3.5677. Earlier versions of this file used q_N directly,
%  which set theta_1 = 1.1363 and made panels B, C, D the rate of a lone fish
%  false-alarming a third of the time. Panel A was unaffected, since it reads
%  the threshold from the solver struct rather than rebuilding it.
%
%  Consequence for panel D: at q_1 = 0.0282 the observed point is
%  alpha_req(7) = 0.542 against L_inf(7) = 0.451 and K_max = 35.9. Under the
%  old q_N anchor it read 1.572 and 3.56, i.e. inadmissible. The reversal is
%  the point of the panel, so check it after any change here.
%
%  NO FIGURE, BUT REPORTED: report_scaling prints the threshold-scaling
%  worked example quoted in the main text (Sec. threshold_scaling), so that
%  those numbers come off the same primitives as the panels rather than
%  being recomputed by hand. Two of the three errors caught in this project
%  came from separately recomputed quantities.
%
%  Usage:  bayes = solve_bayes([2 5 10 20 40 100], 3.5678);
%          out   = make_fig4(bayes, 'output');
%
%  ZPK 2026

if nargin < 1, bayes = []; end
if nargin < 2 || isempty(outdir), outdir = 'output'; end

P = struct();
P.qN_obs = 0.3210;                       % SHOAL-level rate, M-pooled
P.Mhat   = 13.5241;                      % pooling count, from make_fig5
P.q1_obs = 1 - (1-P.qN_obs)^(1/P.Mhat);  % 0.028220, PER-UNIT
% Propagated from the Mhat bootstrap interval [8.8249, 29.4032] (corrected
% file|bout cluster key, 4000 replicates). The earlier [9.9, 26.0] came from
% the pre-fix key and is superseded; q_1 is DECREASING in M, so the upper M
% gives the lower q_1 and the endpoints enter reversed.
P.q1_CI  = [1-(1-P.qN_obs)^(1/29.4032), 1-(1-P.qN_obs)^(1/8.8249)];
P.theta1 = -log(P.q1_obs);               % 3.5677, the K = 1 anchor
P.K_obs  = 7;
P.alpha_hat = 0.9489;                    % from make_fig5, for report_scaling
P.outdir = outdir;
P.LW     = 3.5;  P.FSZ = 26;  P.FIG = [22 19];
P.greens = [0.75 0.90 0.78; 0.55 0.80 0.60; 0.36 0.68 0.44;
            0.20 0.55 0.32; 0.10 0.42 0.24; 0.03 0.28 0.15];
P.blues  = [0.72 0.86 0.94; 0.35 0.60 0.82; 0.08 0.24 0.48];

fprintf(['\n  q_N = %.4f (shoal) -> q_1 = %.4f (per unit) at Mhat = %.2f\n' ...
         '  theta_1 = -log(q_1) = %.4f   (K = 1 anchor)\n'], ...
        P.qN_obs, P.q1_obs, P.Mhat, P.theta1);
out = struct('theta1', P.theta1);
[out.N_A, out.tail_A]   = panelA(bayes, P);
out.T1_B                = panelB(P);
[out.frontier, out.N_C] = panelC(P);
out.D                   = panelD(P);
out.scaling             = report_scaling(P);
fprintf('\n');
end

%% ======================= A: total survival drift ========================
function [Ns, tails] = panelA(bayes, P)
Ns = []; tails = [];
if isempty(bayes)
    fprintf('  panel A skipped (no solver output passed)\n'); return
end
f  = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
nb = numel(bayes); Ns = nan(1,nb); tails = nan(1,nb);
for i = 1:nb
    N = bayes(i).N; t = bayes(i).t(:); lam = bayes(i).lam(:); thN = bayes(i).thetaN;
    L = (N-1)*lam;  x = t/thN;
    plot(ax, x, L, '-', 'Color', P.greens(min(i,size(P.greens,1)),:), 'LineWidth', P.LW);
    % Tail by EXTRAPOLATION on a window common in rescaled time, not a
    % window mean: fit L ~ A + B*(t/theta)^(-2/3), take the intercept.
    w  = x >= 3 & x <= 6.8 & isfinite(L);
    cf = [ones(nnz(w),1), x(w).^(-2/3)] \ L(w);
    Ns(i) = N; tails(i) = cf(1);
end
plot(ax, xlim(ax), [1 1], ':', 'Color',[0.4 0.4 0.4], 'LineWidth', 2);
xlabel(ax,'rescaled time $t/\theta_N$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'survival correction $\lambda^{(N)}_{\rm tot}(t)$', ...
       'Interpreter','latex','FontSize',P.FSZ+8);
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
pos = get(ax,'Position');
axI = axes(f,'Position',[pos(1)+0.50*pos(3), pos(2)+0.52*pos(4), 0.40*pos(3), 0.36*pos(4)]);
hold(axI,'on'); box(axI,'off')
Ng = logspace(log10(2), log10(200), 300);
plot(axI, Ng, linf(Ng), '--', 'Color',[0.2 0.2 0.2], 'LineWidth', 2.5);
plot(axI, Ns, tails, 'o', 'MarkerSize', 10, 'LineWidth', 2.5, ...
     'MarkerEdgeColor',[0.03 0.28 0.15], 'MarkerFaceColor','w');
set(axI,'XScale','log','TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ-10);
xlabel(axI,'$N$','Interpreter','latex','FontSize',P.FSZ-4);
ylabel(axI,'tail $L_\infty$','Interpreter','latex','FontSize',P.FSZ-4);
ylim(axI,[0 1.05]);
fprintf('  panel A: tails / L_inf ='); fprintf(' %.3f', tails./linf(Ns)); fprintf('\n');
%  POST: tails approach L_inf FROM BELOW, tightening with N. The draft's
%        "from above, exceeding by ~20%%" was a window-mean artifact.
export(f, 'fig4A', P);
end

%% ==================== B: first-departure density ========================
function T1 = panelB(P)
% One figure, two stacked axes sharing an x-axis: top naive, bottom Bayesian
% at the self-consistent alpha = L_inf(N). Independent y-limits (the naive
% densities are taller); a single y-label spans both tiles.
f  = figure('Units','centimeters','Position',[1 1 P.FIG(1) P.FIG(2)*1.10],'Color','w');
tl = tiledlayout(f, 2, 1, 'TileSpacing','tight','Padding','compact');

t  = linspace(1e-3, 1.3*P.theta1, 6000)';   % plotting grid
% The mean is a separate, LONGER integral. At theta_1 = 3.57 the N = 1
% survival still carries ~5% of its mass past t = 8, so computing T1 by
% trapz on the plotting grid truncates it (3.376 against the exact 3.5677 =
% theta by Wald). N > 1 is unaffected, since the minimum concentrates at
% small t -- which is why this stayed hidden at the old theta_1 = 1.136.
tI = linspace(1e-4, 60*P.theta1, 200000)';  % integration grid
NB = [1 10 100];  T1 = nan(2,numel(NB));
axB = gobjects(1,2);
for j = 1:2
    axB(j) = nexttile(tl); hold(axB(j),'on'); box(axB(j),'off')
    for i = 1:numel(NB)
        N = NB(i);
        % Lower tile uses the model's own self-consistent rate at each N,
        % alpha = L_inf(N), not a single fitted value: panel A establishes
        % that the Bayesian rate varies with N, and alpha_hat is empirical,
        % belonging to Fig. 5. L_inf(1) = 0, so N = 1 coincides in both tiles.
        if j == 1, a = 0; else, a = linf(N); end
        Sv = surv(t, 1-a, P.theta1);  fv = igd(t, 1-a, P.theta1);
        plot(axB(j), t, N*fv.*Sv.^(N-1), '-', 'Color', P.blues(i,:), 'LineWidth', P.LW);
        T1(j,i) = trapz(tI, surv(tI, 1-a, P.theta1).^N);
    end
    set(axB(j),'TickDir','out','TickLabelInterpreter','latex', ...
               'FontSize',P.FSZ-2,'LineWidth',2);
    % x-limit scales with the threshold: at theta_1 = 3.57 the N = 1 mean is
    % 3.57, so the old [0 1] window cut off the naive tile entirely.
    axis(axB(j),[0 1.2*P.theta1 0 2.2]);
end
set(axB(1),'XTickLabel',[]);                       % shared x-axis
% MODEL UNITS, not seconds. The model carries no timescale of its own: the
% only place seconds appear is Fig. 5, where the Pacher latencies set one.
% The caption states the axis is truncated at 1.2*theta_1, so the label must
% not promise a physical unit the abscissa does not carry.
xlabel(tl,'time $t$','Interpreter','latex','FontSize',P.FSZ+6);
ylabel(tl,'first-departure density, $H=1$','Interpreter','latex','FontSize',P.FSZ+6);

%  POST: top tile is alpha = 0 (naive), bottom is alpha = L_inf(N)
%        (Bayesian); line darkness is N = 1, 10, 100. No legend by design --
%        annotate in assembly.
fprintf('  panel B: mean T_(1) naive    ='); fprintf(' %7.3f', T1(1,:)); fprintf('\n');
fprintf('           mean T_(1) Bayesian ='); fprintf(' %7.3f', T1(2,:)); fprintf('\n');
export(f, 'fig4B', P);
end

%% ======================= C: Pareto frontier =============================
function [front, NC] = panelC(P)
f  = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
th_g = linspace(0.2, 3.5*P.theta1, 200);  a_g = linspace(0, 0.95, 35);  NC = [1 5 20];
tC = linspace(1e-4, 400, 150000)';  front = cell(1,numel(NC));
for i = 1:numel(NC)
    N = NC(i); Pts = zeros(0,2);
    for a = a_g
        if N == 1 && a > 0, continue, end     % alpha is inert for a lone agent
        for th = th_g
            Pts(end+1,:) = [trapz(tC, surv(tC,1-a,th).^N), ...
                            1-(1-exp(-(1+a)*th))^N]; %#ok<AGROW>
        end
    end
    Pts = sortrows(Pts,1); keep = false(size(Pts,1),1); best = inf;
    for r = 1:size(Pts,1)
        if Pts(r,2) < best, keep(r) = true; best = Pts(r,2); end
    end
    front{i} = Pts(keep,:);
    plot(ax, front{i}(:,1), front{i}(:,2), '-', 'Color', P.blues(i,:), 'LineWidth', P.LW);
end
set(ax,'YScale','log','TickDir','out','TickLabelInterpreter','latex', ...
       'FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'detection time $\bar T_{(1)}$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'group false alarm rate $q_N$','Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[0 2.5*P.theta1]); ylim(ax,[1e-4 1]);
legend(ax, arrayfun(@(n) sprintf('$N=%d$',n), NC, 'uni',0), ...
       'Interpreter','latex','Box','off','Location','northeast','FontSize',P.FSZ-6);
fprintf('  panel C: frontier points ='); fprintf(' %d', cellfun(@(c) size(c,1), front)); fprintf('\n');
export(f, 'fig4C', P);
end

%% =================== D: required discounting rate =======================
function D = panelD(P)
f  = figure('Units','centimeters','Position',[1 1 P.FIG(1)+2 P.FIG(2)],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
q_g = logspace(log10(1e-3), log10(0.5), 420);
K_g = logspace(log10(1),    log10(100), 380);
[Q,K] = meshgrid(q_g, K_g);
A = areq(K,Q); Aplot = A; Aplot(A > 1 | A < 0) = NaN;
pc = pcolor(ax, q_g, K_g, Aplot);
set(pc,'EdgeColor','none','FaceColor','texturemap')
colormap(ax, plasma()); caxis(ax,[0 1])
Kmax = log(1-q_g)./log(1-q_g.^2);                    % exact alpha = 1 boundary
patch(ax, [q_g fliplr(q_g)], [Kmax 1e4*ones(size(q_g))], ...
      [0.88 0.88 0.88], 'EdgeColor','none','FaceAlpha',0.92)
plot(ax, q_g, Kmax, '-', 'Color',[0.15 0.15 0.15], 'LineWidth', P.LW);
[c1,h1] = contour(ax, q_g, K_g, A, 0.1:0.1:0.9, 'LineColor',[1 1 1], 'LineWidth', 1.6);
clabel(c1, h1, 'Color',[1 1 1], 'FontSize', 18, 'LabelSpacing', 600)
q_linf = arrayfun(@(k) q_at_alpha(k, linf(k)), K_g);
plot(ax, q_linf, K_g, '-', 'Color',[0.10 0.45 0.85], 'LineWidth', P.LW);
% observed operating point, with the Mhat-propagated interval on q_1
plot(ax, P.q1_CI, [P.K_obs P.K_obs], '-', 'Color',[0.35 0.10 0.55], 'LineWidth', P.LW);
plot(ax, P.q1_obs, P.K_obs, 'o', 'MarkerFaceColor',[0.35 0.10 0.55], ...
     'MarkerEdgeColor','w', 'MarkerSize', 16, 'LineWidth', 2.5);
set(ax,'XScale','log','YScale','log','TickDir','out', ...
       'TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlim(ax,[q_g(1) q_g(end)]); ylim(ax,[1 100]);
set(ax,'XTick',[1e-3 1e-2 1e-1 0.5],'XTickLabel',{'0.001','0.01','0.1','0.5'});
set(ax,'YTick',[1 2 5 10 20 50 100],'YTickLabel',{'1','2','5','10','20','50','100'});
xlabel(ax,'false alarm rate held constant, $q_1$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'attended neighbors $K$','Interpreter','latex','FontSize',P.FSZ+8);
cb = colorbar(ax); cb.Label.String = 'required $\alpha$';
cb.Label.Interpreter = 'latex'; cb.Label.FontSize = P.FSZ; cb.TickLabelInterpreter = 'latex';
%  POST labels:
%    grey wedge  "discounting alone insufficient: threshold must rise"
%                (NOT "inadmissible" -- the rate is attainable, just not at
%                fixed theta)
%    black curve "$K_{\max}\simeq 1/q_1$"
%    blue curve  "$\alpha_{\rm req} = L_\infty(K)$", with the reading
%                glossed either side: left, individually Bayesian discounting
%                already suffices; right, it does not.
%    purple dot  the sulphur molly operating point, $(q_1, K) = (0.028, 7)$,
%                bar the interval propagated from Mhat. It sits to the RIGHT
%                of the blue curve, so the required rate exceeds what
%                individually Bayesian updating supplies -- the same
%                conclusion the empirical fit reaches, from the false alarm
%                rate alone.
fprintf(['  panel D: alpha_req(K=%d, q_1=%.4f) = %.3f  vs  L_inf(%d) = %.3f\n' ...
         '           K_max(obs) = %.2f\n'], ...
        P.K_obs, P.q1_obs, areq(P.K_obs,P.q1_obs), P.K_obs, linf(P.K_obs), ...
        log(1-P.q1_obs)/log(1-P.q1_obs^2));
% the two endpoints of the drawn interval bar, for the caption
fprintf('           q_1 interval from Mhat = [%.4f, %.4f]\n', P.q1_CI(1), P.q1_CI(2));
export(f, 'fig4D', P);
D = struct('q',q_g,'K',K_g,'alpha_req',A,'Kmax_of_q',Kmax,'q_at_Linf',q_linf);
end

%% ============ threshold scaling worked example (no figure) ==============
function S = report_scaling(P)
% Sec. threshold_scaling quotes: at K = 7 the naive rule needs theta_K =
% 5.51 where a discounter at alpha_hat needs 2.83, and the group detection
% time at MATCHED false alarm rate falls from 1.79 to 0.96 at N = 20. The
% thresholds were already right; the two times had no console line behind
% them, which is what this block supplies. Same primitives as panel B, same
% long integration grid, so no separate recomputation.
K = P.K_obs; a = P.alpha_hat; Nrep = [7 20];
thK = [(P.theta1 + log(K))/1, (P.theta1 + log(K))/(1 + a)];   % naive, discounted
tI  = linspace(1e-4, 4000, 400000)';
S = struct('K',K,'alpha',a,'theta_K',thK,'N',Nrep,'T1',nan(2,numel(Nrep)));
for j = 1:2
    for i = 1:numel(Nrep)
        S.T1(j,i) = trapz(tI, surv(tI, 1-(j-1)*a, thK(j)).^Nrep(i));
    end
end
fprintf('  scaling: K = %d, alpha_hat = %.4f\n', K, a);
fprintf('           theta_K naive = %.4f   discounted = %.4f\n', thK(1), thK(2));
for i = 1:numel(Nrep)
    fprintf('           N = %3d:  T_(1) naive = %.4f   discounted = %.4f   (%.2fx faster)\n', ...
            Nrep(i), S.T1(1,i), S.T1(2,i), S.T1(1,i)/S.T1(2,i));
end
%  POST: these are at MATCHED group false alarm rate, since theta_K is set by
%        the pinning relation in each case. Contrast with panel B, which
%        holds theta fixed instead and therefore shows discounting COSTING
%        speed. Both statements are in the text and they are not in conflict.
end

%% ============================ helpers ===================================
function export(f, name, P)
if exist(P.outdir,'dir')
    exportgraphics(f, fullfile(P.outdir,[name '.pdf']), 'ContentType','vector');
    exportgraphics(f, fullfile(P.outdir,[name '.png']), 'Resolution',400);
    fprintf('           wrote %s.pdf / .png\n', name);
else
    fprintf('           outdir "%s" not found -- %s not exported\n', P.outdir, name);
end
end

function f = igd(t, mu, th)
t = t(:); f = zeros(size(t)); ok = t > 0; tt = t(ok);
f(ok) = th ./ sqrt(4*pi*tt.^3) .* exp(-(th - mu*tt).^2 ./ (4*tt));
end

function S = surv(t, mu, th)
t = t(:); S = ones(size(t)); ok = t > 0; tt = t(ok);
S(ok) = ncdf((th - mu*tt)./sqrt(2*tt)) - exp(mu*th).*ncdf((-th - mu*tt)./sqrt(2*tt));
S = min(max(S,0),1);
end

function p = ncdf(x)
% Local standard normal CDF. Every other script in the repo carries one of
% these; this file was the last caller of the Statistics Toolbox normcdf,
% which meant the figure pipeline was not toolbox-free as the README claims.
p = 0.5*erfc(-x./sqrt(2));
end

function L = linf(N)
k = N - 1;
L = ((k + 2) - 2*sqrt(k + 1)) ./ max(k, eps);
L(N <= 1) = 0;
end

function a = areq(K, q)
a = log(1 ./ (1 - (1-q).^(1./K))) ./ log(1./q) - 1;
a(K == 1) = 0;
end

function q = q_at_alpha(K, a)
if K <= 1 || a <= 0, q = NaN; return, end
g = @(lq) areq(K, exp(lq)) - a;
lo = log(1e-14); hi = log(0.5 - 1e-9);
if g(lo)*g(hi) > 0, q = NaN; return, end
q = exp(fzero(g, [lo hi]));
end

function m = plasma()
m = [0.050 0.030 0.528; 0.179 0.016 0.592; 0.281 0.012 0.632;
     0.367 0.012 0.659; 0.445 0.026 0.672; 0.519 0.052 0.674;
     0.586 0.089 0.661; 0.647 0.130 0.638; 0.700 0.171 0.609;
     0.749 0.215 0.576; 0.792 0.260 0.540; 0.832 0.308 0.502;
     0.866 0.358 0.462; 0.897 0.411 0.422; 0.922 0.467 0.381;
     0.943 0.526 0.341; 0.958 0.589 0.302; 0.968 0.654 0.266;
     0.973 0.722 0.234; 0.975 0.793 0.210; 0.973 0.864 0.199;
     0.970 0.936 0.212];
end