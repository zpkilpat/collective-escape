%% MAKE_FIG3  The cost of ignoring a neighbour's hesitation (dyad).
%
%  A  a false-alarm cascade under H=0, naive vs discounting on shared noise.
%     Uses an ILLUSTRATIVE alphaA and thA, not the fitted values -- see below.
%     NOTE: the heuristic kick theta + alpha*t exactly cancels the accumulated
%     discount alpha*t, so both rules land the survivor in the SAME post-kick
%     position. What differs is the TRIGGER: under H=0 the discounting agent
%     drifts at -(1+alpha) and false-alarms later. That delay, not a weaker
%     kick, is where the accuracy gain comes from.
%  B  cascade probability ratio q_dyad/q^2 vs kick size, theta = 1,2,3.
%  C  cascade delay vs kick size.
%  D  Bayesian survival correction lambda(t) against the uncoupled single-agent
%     hazard gap, with the closed-form ceiling L_inf(2) the heuristic adopts.
%  E  survival functions, naive / heuristic / Bayesian.
%  F  speed-accuracy frontier: naive, the lower envelope of the heuristic
%     family swept over alpha, and the exact Bayesian dyad by Monte Carlo.
%
%  PANEL SELECTION. Set PANELS below to any subset of 'ABCDEF'. The fixed
%  point is solved only when D or E is requested; F is independent of it and
%  runs its own per-theta solves. Rerunning 'DEF' leaves fig3A/B/C untouched
%  on disk.
%
%  WHICH PANELS USE theta1. Only D and E. Panel A runs at an illustrative
%  thA = 1.8 chosen for legibility, B and C at theta = 1, 2, 3, and F sweeps
%  theta as its abscissa. So a change of anchor requires rerunning 'DE' only.
%  The caption must not declare theta_1 globally: it applies to D and E, and
%  A's illustrative threshold has to be stated separately.
%
%  PANEL D HORIZON. Earlier versions stopped at t=8, which is t/theta = 2.4,
%  i.e. BEFORE the plateau window [3.0,6.8]*theta in which the tail constant
%  is measured, so the relaxation read as an unresolved transient. The axis
%  now runs to TSHOW = 7*theta1, the top of that window, inside the solve
%  horizon 8*theta1. The excess lambda - L_inf falls only as roughly 1/t, so
%  a visually converged panel would need t ~ 30-50 theta and would be mostly
%  flat line; the approach is reported as "from above" and left unquantified.
%  The tail exponent is printed as a console diagnostic only: over the 0.36
%  decades available it is not resolved well enough to publish.
%
%  PANEL E HORIZON. E keeps the original TSHOW_E = 8. It is a match-quality
%  panel and its content lives where survival runs from 1 to about 0.1;
%  stretching it to 7*theta1 crushes all the structure into the left third
%  and leaves half the axis at zero. D and E deliberately differ.
%
%  PANEL F ENVELOPE. The old version flattened the (alpha,theta) family into
%  one vector, binned along T_(2), and plotted each bin minimum at its own
%  jittered T. With about four points per bin the vertices were unevenly
%  spaced, so the polyline zigzagged, and bins where a single branch still had
%  support produced the drop near T_(2)=3.9. The family is now kept as a
%  (theta,alpha) matrix, each branch is interpolated in log10(q) onto a common
%  T grid over its own support only, the envelope is the pointwise minimum
%  where at least MINBRANCH branches are defined, and a cumulative minimum
%  enforces the Pareto property.
%
%  PANEL F RESULT. The tuned envelope improves on naive by at most 0.57% across
%  the whole range, and over about a quarter of it the minimising member of the
%  family is alpha = 0 itself. The three rules are therefore on essentially the
%  same frontier at N = 2, but say NEGLIGIBLE rather than UNRESOLVED: the exact
%  Bayesian circles sit several MC standard errors off the naive curve, so the
%  separation is statistically resolved and merely small. The draw order puts
%  naive ON TOP as a thick dashed underlay with the envelope beneath it, so the
%  panel reads as agreement rather than as one curve hiding two others.
%
%  Self-contained. ZPK 2026

clear; close all;

PANELS      = 'F';     % any subset of 'ABCDEF'
USE_CACHE_F = true;      % reuse output/fig3F_mc.mat for the Bayesian MC

want = @(c) any(PANELS == c);

st = paper_style();  outdir = 'output';
if ~exist(outdir,'dir'), mkdir(outdir); end

theta1 = 3.5678;
dt = 0.005;  Tmax = 8*theta1;
t  = (0:dt:Tmax)';  nt = numel(t);

TSHOW   = 7*theta1;      % D: top of the plateau window
TSHOW_E = 8;             % E: decision-relevant window, NOT tied to TSHOW

% The heuristic constant is the closed-form ceiling, not a window mean of
% lambda(t). validate_alpha_choice scores L_inf at 2.7% on the dyad
% first-departure time; a decision-window mean is 8.4% and worst on survival
% error. L_inf is also theta-independent, which Panel F exploits. The window
% mean is horizon-dependent besides: at theta1 = 3.5678 it reads 0.2153 over
% [3,6.8]*theta while lambda is still falling at 7*theta. That window mean sits
% 25% ABOVE L_inf and is the artifact retired from Fig. 4A -- never quote it as
% a tail value. The extrapolated tail (solve_bayes) is 0.1636, 4.7% BELOW.
% The 2.7% / 8.4% pair is produced by validate_alpha_choice.m, NOT here; it
% must be rerun at this anchor before those numbers are quoted.
alpha_dyad = linf_ceiling(2);

if want('D') || want('E')
    fprintf('Fig 3: solving Bayesian dyad fixed point...\n');
    [~, lambda, info] = solve_lambda_fixedpoint(2, theta1, dt, Tmax);
    fprintf('  alpha_dyad = L_inf(2) = %.4f  (numerical tail %.4f)\n', ...
            alpha_dyad, info.tail);
    fprintf('  solve horizon %.2f (= %.2f theta); D to %.2f, E to %.2f\n', ...
            Tmax, Tmax/theta1, TSHOW, TSHOW_E);
end

%% ---- A: a false-alarm cascade, naive vs discounting --------------------
if want('A')
thA = 1.8;  dtA = 0.004;  TA = 12;  tA = (0:dtA:TA)';  nA = numel(tA);
% ILLUSTRATIVE alpha, not alpha_dyad. The whole visual content of this panel is
% the within-pair gap alpha*t between a rule and its discounted counterpart. At
% alpha_dyad = 0.17 and a crossing near t = 2.8 that gap is 0.48, against a
% y-range set by the H=0 excursion, so it renders as two lines on top of each
% other. Raising alpha widens the gap but also makes a crossing rarer, since
% the naive path must now reach thA + alphaA*t under a drift of -1; alphaA much
% above 0.4 empties the search. Fig. 1C already declares an illustrative
% alpha = 0.8, so the caption convention exists. State alphaA there.
alphaA = 0.4;
best = [];  bestscore = -inf;
for seed = 1:200000
    rng(seed);
    e1 = sqrt(2*dtA)*randn(nA,1);  e2 = sqrt(2*dtA)*randn(nA,1);
    x1n = cumsum(-dtA + e1);  x2n = cumsum(-dtA + e2);
    x1b = x1n - alphaA*tA;  x2b = x2n - alphaA*tA;      % same noise, discounted

    kn = find(x1n >= thA, 1);   kb = find(x1b >= thA, 1);
    if isempty(kn) || isempty(kb) || kb <= kn, continue; end
    if any(x2n(1:kn) >= thA) || any(x2b(1:kb) >= thA), continue; end
    if x2n(kn) < 0.15, continue; end                 % naive kick must clear
    delay = tA(kb) - tA(kn);
    if delay < 0.5 || delay > 4.0 || tA(kn) < 1.0 || tA(kb) > 9, continue; end

    % TWO separations exist here and they compete. Between AGENTS (x1 against
    % x2) is incidental; within a PAIR (a rule against its discounted twin) is
    % the point. The old score MAXIMISED the between-agent term, which is what
    % drove agent 2 down to -2.7, set a y-range of 5.7 units, and left the
    % within-pair gap at under 9% of the axis. Now the between-agent spread is
    % targeted at a moderate value rather than maximised, and the depth of the
    % excursion is penalised outright, since it is what sets the y-range.
    sep   = mean(abs(x1n(1:kn) - x2n(1:kn)));
    depth = -min([x1b(1:kb); x2b(1:kb)]);
    if depth > 2.0, continue; end
    score = -abs(delay - 1.5) - 0.3*abs(sep - 1.0) - 0.6*depth;
    if score > bestscore
        bestscore = score;
        best = struct('kn',kn,'kb',kb,'x1n',x1n,'x2n',x2n,'x1b',x1b,'x2b',x2b);
    end
end
if isempty(best)
    error(['Panel A: no legible realisation found. alphaA = %.2f may be too ' ...
           'large -- the naive path must reach thA + alphaA*t under drift -1. ' ...
           'Lower alphaA, raise the seed count, or relax the depth cap.'], alphaA);
end

kn = best.kn;  kb = best.kb;
tn = tA(kn);   tb = tA(kb);
kapN = thA;                      % naive kick
kapB = thA + alphaA*tb;          % heuristic kick, grows with departure time
xEnd = min(TA, tb + 1.2);

fA = figure('Units','centimeters','Position',[1 1 st.FIG],'Color','w');
axA = axes(fA); hold(axA,'on');
yline(axA, thA, '-', 'Color',[0 0 0], 'LineWidth',st.LWt);

% Shade each rule against its discounted twin. A filled area reads at a glance
% where two nearly parallel lines do not, and it is the same device Fig. 1C
% uses for the discounted silence.
fill(axA, [tA(1:kn); flipud(tA(1:kn))], ...
          [best.x1b(1:kn); flipud(best.x1n(1:kn))], ...
     [0.55 0.72 0.55], 'FaceAlpha', 0.30, 'EdgeColor','none');
fill(axA, [tA(1:kn); flipud(tA(1:kn))], ...
          [best.x2b(1:kn); flipud(best.x2n(1:kn))], ...
     [0.55 0.72 0.55], 'FaceAlpha', 0.30, 'EdgeColor','none');

plot(axA, tA(1:kn), best.x1n(1:kn), '-', 'Color',st.col_naiv, 'LineWidth',st.LW);
plot(axA, tA(1:kn), best.x2n(1:kn), '-', 'Color',[0.78 0.78 0.78], 'LineWidth',st.LW);
plot(axA, [tn tn], [best.x2n(kn) best.x2n(kn)+kapN], ':', ...
     'Color',[0.78 0.78 0.78], 'LineWidth',st.LWt);
plot(axA, tn, thA, 'o', 'MarkerSize',12, 'MarkerFaceColor',st.col_naiv, ...
     'MarkerEdgeColor',st.col_naiv);

plot(axA, tA(1:kb), best.x1b(1:kb), '-', 'Color',st.col_bay, 'LineWidth',st.LW);
plot(axA, tA(1:kb), best.x2b(1:kb), '-', 'Color',[0.55 0.80 0.62], 'LineWidth',st.LW);
plot(axA, [tb tb], [best.x2b(kb) best.x2b(kb)+kapB], ':', ...
     'Color',[0.55 0.80 0.62], 'LineWidth',st.LWt);
plot(axA, tb, thA, 'o', 'MarkerSize',12, 'MarkerFaceColor',st.col_bay, ...
     'MarkerEdgeColor',st.col_bay);

yb = thA + 0.32;
plot(axA, [tn tb], [yb yb], '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);
plot(axA, [tn tn], [thA yb], ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1);
plot(axA, [tb tb], [thA yb], ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1);
text(axA, 0.5*(tn+tb), yb+0.16, 'delay', 'Interpreter','latex', ...
     'FontSize',st.FSZ-10, 'Color',[0.35 0.35 0.35], 'HorizontalAlignment','center');
text(axA, 0.15, thA-0.45, 'naive',       'Interpreter','latex', ...
     'FontSize',st.FSZ-10, 'Color',st.col_naiv);
text(axA, 0.15, thA-1.05, 'discounting', 'Interpreter','latex', ...
     'FontSize',st.FSZ-10, 'Color',st.col_bay);

apply_axes(axA, 'time $t$', 'belief $\xi(t)$');
xlim(axA,[0 xEnd]);
ylim(axA,[min([best.x1b(1:kb); best.x2b(1:kb)])-0.4, thA+0.9]);
export_panel(fA, 'fig3A', outdir);
yr = diff(ylim(axA));
fprintf(['  Panel A: alphaA = %.2f, naive trigger t=%.2f, discounting trigger ' ...
         't=%.2f (delay %.2f)\n'], alphaA, tn, tb, tb-tn);
fprintf('           within-pair gap at crossing %.2f, y-range %.2f (%.0f%% of axis)\n', ...
        alphaA*tb, yr, 100*alphaA*tb/yr);
end

%% ---- B, C: naive cascade observables vs kick size ----------------------
if want('B') || want('C')
theta_BC = [1 2 3];
shades   = [0.70 0.70 0.70; 0.40 0.40 0.40; 0 0 0];
kap_rel  = linspace(0.02, 4, 60);
kap_mc   = [0.25 0.75 1.5 2.3 3.2];

fB = figure('Units','centimeters','Position',[3 1 st.FIG],'Color','w');
axB = axes(fB); hold(axB,'on');
fC = figure('Units','centimeters','Position',[5 1 st.FIG],'Color','w');
axC = axes(fC); hold(axC,'on');

for j = 1:numel(theta_BC)
    th = theta_BC(j);  q = exp(-th);
    ratio = zeros(size(kap_rel));  delay = zeros(size(kap_rel));
    for i = 1:numel(kap_rel)
        [qc, ct] = dyad_cascade_analytic(th, kap_rel(i)*th);
        ratio(i) = qc/q^2;   delay(i) = ct/th;
    end
    plot(axB, kap_rel, ratio, '-', 'Color',shades(j,:), 'LineWidth',st.LW);
    plot(axC, kap_rel, delay, '-', 'Color',shades(j,:), 'LineWidth',st.LW);

    rmc = zeros(size(kap_mc));  dmc = zeros(size(kap_mc));
    for i = 1:numel(kap_mc)
        [rmc(i), dmc(i)] = dyad_mc(th, kap_mc(i)*th, 1e4, 1e-3, 30*th);
    end
    plot(axB, kap_mc, rmc/q^2, 'o', 'Color',shades(j,:), 'MarkerFaceColor','w', ...
         'MarkerSize',9, 'LineWidth',2);
    plot(axC, kap_mc, dmc/th,  'o', 'Color',shades(j,:), 'MarkerFaceColor','w', ...
         'MarkerSize',9, 'LineWidth',2);
end
yline(axB, 1, ':', 'Color',[0.85 0.3 0.3], 'LineWidth',st.LWt);
text(axB, 2.4, 1.25, 'independent', 'Interpreter','latex', ...
     'FontSize',st.FSZ-12,'Color',[0.85 0.3 0.3]);
set(axB,'YScale','log');
apply_axes(axB, 'jump $\kappa/\theta$', 'cascade ratio $q_{\rm dyad}/q^2$');
xlim(axB,[0 4]);
export_panel(fB, 'fig3B', outdir);
apply_axes(axC, 'jump $\kappa/\theta$', '$[\bar{T}_{(2)}-\bar{T}_{(1)}]/\theta$');
xlim(axC,[0 4]);  ylim(axC,[0 0.6]);
export_panel(fC, 'fig3C', outdir);
end

%% ---- D: the survival correction ----------------------------------------
if want('D')
h1s = ddm_ig(t,  1, theta1).h;
h0s = ddm_ig(t, -1, theta1).h;
gap = max(h1s - h0s, 0);

show = t <= TSHOW;

fD = figure('Units','centimeters','Position',[7 1 st.FIG],'Color','w');
axD = axes(fD); hold(axD,'on');
plot(axD, t(show), gap(show),    '--', 'Color',st.col_ref, 'LineWidth',st.LWt);
plot(axD, t(show), lambda(show), '-',  'Color',st.col_bay, 'LineWidth',st.LW);
yline(axD, alpha_dyad, ':', 'Color',st.col_bay, 'LineWidth',st.LWt);
text(axD, 0.62*TSHOW, alpha_dyad-0.035, sprintf('$L_\\infty(2) = %.2f$', alpha_dyad), ...
     'Interpreter','latex','FontSize',st.FSZ-10,'Color',st.col_bay);
text(axD, 0.10*TSHOW, 0.55, 'single agent', 'Interpreter','latex', ...
     'FontSize',st.FSZ-12,'Color',st.col_ref);
apply_axes(axD, 'time $t$', 'survival correction $\lambda(t)$');
% ylim must contain the single-agent gap; an earlier version clipped it
xlim(axD,[0 TSHOW]);  ylim(axD,[0 max(1.15*max(gap(t > 0.3)), 0.45)]);
export_panel(fD, 'fig3D', outdir);

% console diagnostic only: the tail is not resolved well enough to plot
tail = t >= 3.0*theta1 & t <= TSHOW;
dl   = lambda(tail) - alpha_dyad;
if all(dl > 0)
    p = polyfit(log10(t(tail)), log10(dl), 1);
    fprintf('  Panel D: excess above L_inf falls from %.4f at t=%.1f to %.4f at t=%.1f\n', ...
            dl(1), t(find(tail,1,'first')), dl(end), t(find(tail,1,'last')));
    fprintf('           apparent exponent %.3f over %.2f decades (NOT publishable)\n', ...
            p(1), log10(t(find(tail,1,'last'))/t(find(tail,1,'first'))));
else
    fprintf(2, '  Panel D: lambda reaches or crosses L_inf inside the window.\n');
    fprintf(2, '           min(lambda - L_inf) = %.3e. Check the survival floor.\n', min(dl));
end
end

%% ---- E: survival functions ---------------------------------------------
if want('E')
S1_bay = survival_drift( 1, lambda,     theta1, dt, nt);
S0_bay = survival_drift(-1, lambda,     theta1, dt, nt);
S1_heu = survival_drift( 1, alpha_dyad, theta1, dt, nt);
S0_heu = survival_drift(-1, alpha_dyad, theta1, dt, nt);
S1_nai = ddm_ig(t,  1, theta1).S;
S0_nai = ddm_ig(t, -1, theta1).S;

showE = t <= TSHOW_E;

fE = figure('Units','centimeters','Position',[9 1 st.FIG],'Color','w');
axE = axes(fE); hold(axE,'on');
plot(axE, t(showE), S1_nai(showE), '-',  'Color',st.col_naiv, 'LineWidth',st.LW);
plot(axE, t(showE), S0_nai(showE), '--', 'Color',st.col_naiv, 'LineWidth',st.LW);
plot(axE, t(showE), S1_heu(showE), '-',  'Color',[0.45 0.75 0.55], 'LineWidth',st.LW);
plot(axE, t(showE), S0_heu(showE), '--', 'Color',[0.45 0.75 0.55], 'LineWidth',st.LW);
plot(axE, t(showE), S1_bay(showE), '-',  'Color',st.col_bay, 'LineWidth',st.LWt);
plot(axE, t(showE), S0_bay(showE), '--', 'Color',st.col_bay, 'LineWidth',st.LWt);
apply_axes(axE, 'time $t$', 'survival $S(t\,|\,H)$');
xlim(axE,[0 TSHOW_E]);  ylim(axE,[0 1]);
export_panel(fE, 'fig3E', outdir);
fprintf('  max |S_bayes - S_heur|: H=1 %.3e, H=0 %.3e\n', ...
        max(abs(S1_bay-S1_heu)), max(abs(S0_bay-S0_heu)));
end

%% ---- F: speed-accuracy frontier ----------------------------------------
%  Three objects, as the caption describes:
%    naive     alpha = 0, kick kappa = theta
%    heuristic lower envelope of the (alpha, theta) family
%    Bayesian  exact fixed point, by direct Monte Carlo -- the analytic
%              killed-density formula assumes CONSTANT drift and cannot take
%              the time-varying lambda(t), so MC is the honest route here.
if want('F')
DTF = 0.005;  TFAC = 12;  NXF = 2000;
% Upper bound must exceed the largest theta_B: the Bayesian dyad is SLOWER at
% matched threshold (the anti-drift delays the second departure), so T_(2) at
% theta_B = 4.2 overshoots the naive curve's endpoint at theta_N = 4.4 and the
% offset diagnostic returns NaN for that circle.
theta_N = unique([linspace(0.8, 5.2, 30), theta1]);
alpha_N = unique([linspace(0, 0.75, 21), alpha_dyad]);
NTG = 400;  MINBRANCH = 6;      % 6 kills the right-edge ripple past T ~ 4.2
NSHOW = 6;                      % branches drawn as the background cloud
XMAXF = 4.3;

fprintf('  Panel F: naive sweep...\n');
[qn, Tn] = deal(nan(size(theta_N)));
for j = 1:numel(theta_N)
    [qn(j), Tn(j)] = dyad_frontier_heuristic(theta_N(j), 0, DTF, TFAC, NXF);
end

fprintf('  Panel F: heuristic (alpha,theta) family, %d alpha x %d theta...\n', ...
        numel(alpha_N), numel(theta_N));
Hq = nan(numel(theta_N), numel(alpha_N));
HT = nan(numel(theta_N), numel(alpha_N));
for ia = 1:numel(alpha_N)
    for j = 1:numel(theta_N)
        [Hq(j,ia), HT(j,ia)] = dyad_frontier_heuristic(theta_N(j), alpha_N(ia), ...
                                                       DTF, TFAC, NXF);
    end
    fprintf('    alpha=%.4f done\n', alpha_N(ia));
end

% envelope: interpolate log10(q) per branch on its own support, then min
T2g = linspace(min(HT(:)), max(HT(:)), NTG);
Q   = nan(NTG, numel(alpha_N));
for ia = 1:numel(alpha_N)
    [t2, o] = sort(HT(:,ia));  lq = log10(Hq(o,ia));
    g = isfinite(t2) & isfinite(lq);
    if nnz(g) < 4, continue; end
    in = T2g >= min(t2(g)) & T2g <= max(t2(g));
    Q(in,ia) = interp1(t2(g), lq(g), T2g(in), 'pchip');
end
keep = sum(~isnan(Q), 2) >= MINBRANCH;
eT   = T2g(keep);
eQ   = 10.^cummin(min(Q(keep,:), [], 2, 'omitnan'));   % Pareto, left to right

fprintf('  Panel F: envelope over T_(2) in [%.2f, %.2f], %d points, >= %d branches\n', ...
        eT(1), eT(end), numel(eT), MINBRANCH);

% which alpha attains the envelope, and by how much it beats naive
[~, iwin] = min(Q(keep,:), [], 2, 'omitnan');
fprintf('  Panel F: minimising alpha in [%.4f, %.4f]; alpha=0 optimal on %.0f%% of the range\n', ...
        min(alpha_N(iwin)), max(alpha_N(iwin)), 100*mean(alpha_N(iwin) == 0));
qn_at_eT = interp1(Tn, qn, eT, 'pchip', NaN);
gapfrac  = eQ(:)./qn_at_eT(:);
fin = isfinite(gapfrac);
fprintf('  Panel F: envelope/naive ratio in [%.4f, %.4f] -- tuning gains at most %.2f%%\n', ...
        min(gapfrac(fin)), max(gapfrac(fin)), 100*(1-min(gapfrac(fin))));

fprintf('  Panel F: exact Bayesian dyad by Monte Carlo...\n');
theta_B = unique([1.2 1.8 2.4 3.0 3.6 4.2, theta1]);
cacheF  = fullfile(outdir, 'fig3F_mc.mat');
haveF   = false;
if USE_CACHE_F && exist(cacheF,'file')
    C = load(cacheF);
    % 'bgk' in the cache marks the Broadie-Glasserman-Kou-corrected run. An
    % older cache lacks it and is rejected, since the uncorrected q_MC is
    % biased low against the analytic naive curve it is plotted with.
    haveF = isfield(C,'theta_B') && isequal(C.theta_B, theta_B) && ...
            isfield(C,'dt') && C.dt == dt && isfield(C,'bgk') && C.bgk;
    if haveF
        qB = C.qB;  TB = C.TB;
        fprintf('    loaded cached MC from %s\n', cacheF);
    end
end
if ~haveF
    [qB, TB] = deal(nan(size(theta_B)));
    for j = 1:numel(theta_B)
        th  = theta_B(j);
        Tm  = 10*th;
        [~, lamj] = solve_lambda_fixedpoint(2, th, dt, Tm, struct('tol',3e-4));
        % Length comes FROM the solver. round(Tm/dt)+1 disagrees by one whenever
        % 10*theta is not a whole number of dt steps, which no entry in the
        % original theta_B triggered but theta_1 = 3.5678 does: 35.678/0.005 is
        % 7135.6, so round gives 7137 against the solver's 7136 and
        % survival_drift rejects the mismatch.
        ntB = numel(lamj);
        S1j = max(survival_drift( 1, lamj, th, dt, ntB), 1e-12);
        S0j = max(survival_drift(-1, lamj, th, dt, ntB), 1e-12);
        Jtj = th - log(S1j./S0j);              % exact Bayesian kick, Eq. (S12)
        [qB(j), TB(j)] = bayes_dyad_mc(th, lamj, Jtj, dt, ntB, 4e4, 100+j);
        fprintf('    theta=%.2f: q_MC=%.3e  T2_MC=%.3f\n', th, qB(j), TB(j));
    end
    bgk = true;  save(cacheF, 'qB', 'TB', 'theta_B', 'dt', 'bgk');
end

% how far the MC circles sit from the naive curve, in units of MC error
qn_at_TB = interp1(Tn, qn, TB, 'pchip', NaN);
se_q     = sqrt(qB.*(1-qB)/4e4);
off      = (qB - qn_at_TB)./se_q;
fprintf('  Panel F: circle offset from naive, in MC standard errors:');
fprintf(' %.1f', off);  fprintf('\n');
fprintf('           same offset as a fraction of q:                ');
fprintf(' %+.1f%%', 100*(qB - qn_at_TB)./qn_at_TB);  fprintf('\n');
% Sigma units and relative size answer different questions. At q ~ 1e-2 with
% 4e4 trials, five standard errors is a ~25% difference in cascade rate: the
% separation is statistically RESOLVED. Whether it is NEGLIGIBLE is the
% relative line, and that is the one the caption should quote. A positive
% offset means the exact Bayesian dyad is WORSE than naive at matched response
% time, which is the paper's own point that individually optimal inference does
% not optimize the collective order-statistic objective -- but state it only if
% the sign survives the BGK correction, which pushes q_MC up and so widens it.

% MATCHED-THRESHOLD comparison, which is a different claim from the frontier.
% Sec. dyad says the Bayesian rule "cuts q_dyad by roughly a third" at fixed
% theta; the frontier compares at matched T_(2) instead. This is the number
% behind that sentence, and it was previously never computed here.
jB = find(abs(theta_B - theta1) < 1e-9, 1);
jN = find(abs(theta_N - theta1) < 1e-9, 1);
if ~isempty(jB) && ~isempty(jN)
    fprintf(['  Panel F: at MATCHED theta = %.4f, naive q_dyad = %.4e, ' ...
             'Bayesian %.4e\n           ratio %.3f (a reduction of %.1f%%)\n'], ...
            theta1, qn(jN), qB(jB), qB(jB)/qn(jN), 100*(1 - qB(jB)/qn(jN)));
end

% display cloud: pick branches FROM the grid (an independent list misses it)
is = unique(round(linspace(1, numel(alpha_N), NSHOW)));

fF = figure('Units','centimeters','Position',[11 1 st.FIG],'Color','w');
axF = axes(fF); hold(axF,'on');
plot(axF, HT(:,is), Hq(:,is), '.', 'Color',[0.90 0.94 0.91], 'MarkerSize',6);
hH = plot(axF, eT, eQ, '-',  'Color',[0.45 0.75 0.55], 'LineWidth',st.LW);
hN = plot(axF, Tn, qn, '--', 'Color',st.col_naiv, 'LineWidth',st.LW+2);
hB = plot(axF, TB, qB, 'o',  'Color',st.col_bay, 'MarkerFaceColor','w', ...
          'MarkerSize',11, 'LineWidth',2.5);
set(axF,'YScale','log');
apply_axes(axF, 'collective response time $\bar{T}_{(2)}$', ...
                'cascade rate $q_{\rm dyad}$');
xlim(axF,[0 XMAXF]);  ylim(axF,[1e-3 1]);
legend(axF, [hN hH hB], {'naive','heuristic','Bayesian'}, ...
       'Interpreter','latex', 'FontSize',st.FSZ-12, ...
       'Location','northeast', 'Box','off');
export_panel(fF, 'fig3F', outdir);
end

fprintf('Fig 3 done (panels %s). alpha_dyad = L_inf(2) = %.4f\n', PANELS, alpha_dyad);

%% ---- local: naive dyad Monte Carlo (panels B, C) -----------------------
function [qc, ct] = dyad_mc(theta, kappa, ntr, dtmc, Tmax)
    ns = round(Tmax/dtmc);  sd = sqrt(2*dtmc);
    out = [0 0];
    for H = [0 1]
        mu = 2*H - 1;
        x = zeros(ntr,2);  hit = nan(ntr,2);  kicked = false(ntr,1);
        for k = 1:ns
            live = isnan(hit);
            if ~any(live(:)), break; end
            x(live) = x(live) + mu*dtmc + sd*randn(nnz(live),1);
            just = live & (x >= theta);
            hit(just) = k*dtmc;
            one = (sum(~isnan(hit),2) == 1) & ~kicked;
            if any(one)
                s1 = one & isnan(hit(:,1));  x(s1,1) = x(s1,1) + kappa;
                s2 = one & isnan(hit(:,2));  x(s2,2) = x(s2,2) + kappa;
                kicked(one) = true;
            end
        end
        if H == 0
            out(1) = mean(all(~isnan(hit),2));
        else
            d = max(hit,[],2) - min(hit,[],2);
            out(2) = mean(d(~isnan(d)));
        end
    end
    qc = out(1);  ct = out(2);
end

%% ---- local: exact Bayesian dyad Monte Carlo (panel F) ------------------
function [q_casc, T2] = bayes_dyad_mc(theta, lambda, Jt, dt, nt, ntr, seed)
% Each agent drifts at mu_H - lambda(t) WHILE ITS NEIGHBOUR IS SILENT. On the
% neighbour's departure at t_k it receives the exact Bayesian kick
% Jt(k) = theta - Lambda^surv(t_k) and its anti-drift switches OFF. Shares no
% code with the Volterra quadrature, so agreement tests the solver.
%
% DISCRETE-BARRIER BIAS. Euler-Maruyama monitors the barrier only at grid
% points and misses excursions between them, so it UNDERSTATES the crossing
% rate. Uncorrected, this MC is compared against an analytic naive curve that
% has no such bias, which is exactly the error found in make_figS2. The
% Broadie-Glasserman-Kou shift moves the simulated barrier inward by
% 0.5826*sigma*sqrt(dt) (0.058 at dt = 5e-3, about 1.6% of theta_1). The
% correction pushes q_MC UP, so any excess of the Bayesian circles over the
% naive curve is a LOWER bound before it is applied.
    rng(seed);
    lambda = lambda(:);  Jt = Jt(:);
    sq  = sqrt(2*dt);
    thr = theta - 0.5826*sqrt(2)*sqrt(dt);   % BGK-corrected crossing level
    out = [0 0];
    for H = [0 1]
        mu = 2*H - 1;
        x1 = zeros(ntr,1);  x2 = zeros(ntr,1);
        d1 = false(ntr,1);  d2 = false(ntr,1);
        k1 = false(ntr,1);  k2 = false(ntr,1);
        t1 = nan(ntr,1);    t2 = nan(ntr,1);
        for k = 1:nt-1
            lam = lambda(k);
            a1 = ~d1;  a2 = ~d2;
            if any(a1)
                dr = mu - lam*double(~d2(a1) & ~k1(a1));
                x1(a1) = x1(a1) + dr*dt + sq*randn(nnz(a1),1);
            end
            if any(a2)
                dr = mu - lam*double(~d1(a2) & ~k2(a2));
                x2(a2) = x2(a2) + dr*dt + sq*randn(nnz(a2),1);
            end
            c1 = ~d1 & x1 >= thr;   c2 = ~d2 & x2 >= thr;
            kk1 = c1 & ~d2 & ~k2;  x2(kk1) = x2(kk1) + Jt(k);  k2(kk1) = true;
            kk2 = c2 & ~d1 & ~k1;  x1(kk2) = x1(kk2) + Jt(k);  k1(kk2) = true;
            n1 = ~d1 & x1 >= thr;  n2 = ~d2 & x2 >= thr;
            t1(n1) = k*dt;  t2(n2) = k*dt;
            d1 = d1 | x1 >= thr;   d2 = d2 | x2 >= thr;
            if all(d1 & d2), break; end
        end
        if H == 0
            out(1) = mean(d1 & d2);
        else
            both = ~isnan(t1) & ~isnan(t2);
            out(2) = mean(max(t1(both), t2(both)));
        end
    end
    q_casc = out(1);  T2 = out(2);
end