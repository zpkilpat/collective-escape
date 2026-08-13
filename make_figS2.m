function out = make_figS2(bayes, outdir)
%MAKE_FIGS2  Where the constant-alpha heuristic departs from the full
%            Bayesian fixed point, and what that costs.
%
%  A  the pointwise discrepancy S_Bayes(t) - S_heur(t) against rescaled time,
%     one curve per N, under H = 1. Shows WHERE the constant-alpha rule
%     departs: the gap opens through the early transient, peaks near the
%     first-passage mode, and decays as the correction relaxes to L_inf.
%     This is the mechanism Sec. S6 describes -- the constant is applied from
%     t = 0 while the true correction first rises through the hump -- and
%     panel B is its summary.
%
%     An earlier version plotted the survival curves themselves, solid against
%     dashed, for both hypotheses. That wasted a third of the panel on the
%     H = 0 pair, which overlaps into a single band near 1 (max|dS| ~ 1.5e-3),
%     and the visible solid/dashed gap in the H = 1 curves read as though it
%     contradicted the claim rather than quantifying it.
%  B  error in the DERIVED OBSERVABLE against N: the relative discrepancy in
%     the group detection time Tbar_(1) = int S(t|H=1)^N dt, computed under
%     the fixed point and under the heuristic. The max pointwise survival gap
%     is shown alongside for contrast.
%
%     WHY NOT max|dS| ALONE. That metric is the wrong summary for a steep
%     monotone curve: a small shift in TIME produces a large vertical gap
%     wherever the survival is falling fast, so it is dominated by the
%     transient and grows with N (2.6e-2, 1.1e-1, 1.9e-1, 3.0e-1, 4.0e-1 at
%     N = 2, 5, 10, 20, 40). It is
%     also not the quantity the heuristic is used for -- nothing in the paper
%     reads S(t) pointwise; everything reads integrals of it. The growth is
%     expected in any case, since the constant L_inf(N) is applied from t = 0
%     while the true correction first rises through the early-time hump,
%     which diverges as N/(log N)^2.
%
%  Both are evaluated at the SAME threshold with the same primitive, so the
%  comparison isolates the constant-versus-time-varying correction alone. The
%  heuristic survival is analytic, the inverse-Gaussian survival at the
%  shifted drifts mu_H - alpha, so only the Bayesian side needs the solver.
%
%  RESULT (theta_1 = 3.568, N = 2, 5, 10, 20, 40). The detection-time error
%  grows as log N, 1.2%% to 6.9%%, comparable to the saddle-point
%  approximation of Eq. (14) which the paper already accepts at ~10%% for
%  N = 100, and like it never used for a reported number: all detection times
%  come from the exact integral. The pointwise gap over the same range runs
%  2.6%% to 40%%.
%
%  INPUT IS OPTIONAL. Called with no argument the fixed point is solved here
%  via solve_bayes over N = 2, 5, 10, 20, 40 at theta_1 = 3.568. Pass
%  run_bayes_nagent output to override, with fields N, t, thetaN, S1, S0.
%
%  Usage:  out = make_figS2;                    % solves everything itself
%          out = make_figS2(bayes, 'output');   % or pass solver output
%
%  ZPK 2026

if nargin < 2 || isempty(outdir), outdir = 'output'; end
if nargin < 1 || isempty(bayes)
    fprintf('  no solver output passed -- solving the fixed point here\n');
    th1 = calib_pool('theta1', calib_pool('invert', 3.93));
    bayes = solve_bayes([2 5 10 20 40], th1);
end
P.outdir = outdir;
P.LW = 3.5;  P.FSZ = 26;  P.FIG = [22 19];
P.greens = [0.75 0.90 0.78; 0.55 0.80 0.60; 0.36 0.68 0.44;
            0.20 0.55 0.32; 0.10 0.42 0.24; 0.03 0.28 0.15];
P.grey = [0.45 0.45 0.45];

out.A = panelA(bayes, P);
out.B = panelB(bayes, P);
fprintf('\n');
end

%% ============== A: survival curves, both rules ==========================
function A = panelA(bayes, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
n = numel(bayes);
D = cell(1,n); Ns = nan(1,n); pk = nan(1,n); tpk = nan(1,n);
for i = 1:n
    b = bayes(i);  a = linf_ceiling(b.N);  t = b.t(:);  x = t/b.thetaN;
    d = b.S1(:) - surv(t, 1-a, b.thetaN);         % H = 1 only
    c = P.greens(min(i+1, size(P.greens,1)), :);
    plot(ax, x, d, '-', 'Color', c, 'LineWidth', P.LW);
    [pk(i), k] = max(abs(d));  tpk(i) = x(k);
    Ns(i) = b.N;  D{i} = [x, d];
end
plot(ax, xlim(ax), [0 0], ':', 'Color',[0.45 0.45 0.45], 'LineWidth', 2);
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'rescaled time $t/\theta_N$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'$S_{\mathrm{Bayes}} - S_{\mathrm{heur}}$', ...
       'Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[0 8]);
fprintf('  panel A: peak |dS| under H=1'); fprintf(' %.3f', pk); fprintf('\n');
fprintf('           at t/theta_N =       '); fprintf(' %.2f', tpk); fprintf('\n');
%  POST: label the curves by N (light to dark). The H = 0 discrepancies are
%        omitted -- flat at ~1.5e-3 for every N, since the survivor drifts
%        away from the barrier and the anti-drift is nearly inert. State that
%        in the caption rather than plotting a null result.
export(f,'figs2A',P);
A = struct('N',Ns,'diff',{D},'peak',pk,'t_peak',tpk);
end

%% ============== B: discrepancy against group size =======================
function B = panelB(bayes, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
n = numel(bayes);
Ns = nan(1,n); e1 = nan(1,n); e0 = nan(1,n); eT = nan(1,n);
for i = 1:n
    b = bayes(i); a = linf_ceiling(b.N); t = b.t(:);
    Sh1 = surv(t,  1-a, b.thetaN);          % heuristic, H = 1
    Sh0 = surv(t, -1-a, b.thetaN);          % heuristic, H = 0
    Ns(i) = b.N;
    e1(i) = max(abs(b.S1(:) - Sh1));        % pointwise, for contrast
    e0(i) = max(abs(b.S0(:) - Sh0));
    % the quantity the model actually reports: group detection time
    TB = trapz(t, b.S1(:).^b.N);
    TH = trapz(t, Sh1.^b.N);
    eT(i) = abs(TH - TB)/TB;
end
plot(ax, Ns, eT, '-o', 'Color', P.greens(5,:), 'MarkerFaceColor', P.greens(5,:), ...
     'LineWidth', P.LW, 'MarkerSize', 11);
plot(ax, Ns, e1, '--s', 'Color', P.grey, 'MarkerFaceColor','w', ...
     'LineWidth', P.LW-1.2, 'MarkerSize', 9);
set(ax,'XScale','log','YScale','log','TickDir','out', ...
       'TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'group size $N$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'relative error','Interpreter','latex','FontSize',P.FSZ+8);
fprintf('  panel B: detection time  rel err'); fprintf(' %.2e', eT); fprintf('\n');
fprintf('           max|dS| H=1             '); fprintf(' %.2e', e1); fprintf('\n');
fprintf('           max|dS| H=0             '); fprintf(' %.2e', e0); fprintf('\n');
%  POST: SOLID green circles are the error in the group detection time
%        Tbar_(1), the quantity the heuristic is used for. DASHED grey
%        squares are the max pointwise survival gap under H=1, shown for
%        contrast: it grows with N because the constant is applied through
%        the transient, which is exactly what Sec. S6 predicts. If the solid
%        series stays flat and small, the heuristic is fit for purpose and
%        the pointwise growth is a metric artifact; if it does not, that is
%        a real limitation and belongs in the text.
export(f,'figs2B',P);
B = struct('N',Ns,'err_T1',eT,'err_H1',e1,'err_H0',e0);
end

%% ============================ helpers ==================================
function S = surv(t, nu, th)
t = t(:); S = ones(size(t)); ok = t > 0; tt = t(ok);
S(ok) = ncdf((th - nu*tt)./sqrt(2*tt)) - exp(nu*th).*ncdf((-th - nu*tt)./sqrt(2*tt));
S = min(max(S,0),1);
end

function p = ncdf(z), p = 0.5*(1+erf(z/sqrt(2))); end


function export(f, name, P)
if exist(P.outdir,'dir')
    exportgraphics(f, fullfile(P.outdir,[name '.pdf']), 'ContentType','vector');
    exportgraphics(f, fullfile(P.outdir,[name '.png']), 'Resolution',400);
    fprintf('           wrote %s.pdf / .png\n', name);
else
    fprintf('           outdir "%s" not found -- %s not exported\n', P.outdir, name);
end
end
