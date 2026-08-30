function out = make_figS2(bayes, outdir)
%MAKE_FIGS2  Cost of the constant-alpha heuristic against the Bayesian fixed point.
%
%   OUT = MAKE_FIGS2() solves the Bayesian fixed point at N = 2, 5, 10, 20,
%   40, 100 and writes two standalone panels to ./output.
%
%   OUT = MAKE_FIGS2(BAYES) uses a precomputed solver struct array instead,
%   with fields N, t, thetaN, S1, S0.
%
%   OUT = MAKE_FIGS2(BAYES, OUTDIR) writes to OUTDIR. Pass [] for BAYES to
%   solve here and still choose the directory.
%
%   The heuristic replaces the time-varying survival correction lambda^(N)(t)
%   by the constant alpha = L_inf(N). Both sides are evaluated at the same
%   threshold, and the heuristic survival is analytic (inverse Gaussian at the
%   shifted drifts mu_H - alpha), so only the Bayesian side needs the solver
%   and the comparison isolates the constant-versus-time-varying correction.
%
%   Panels
%     A   pointwise discrepancy S_Bayes(t) - S_heur(t) under H = 1 against
%         rescaled time t/theta_N, one curve per N. The gap opens through the
%         early transient, where the constant is imposed while the true
%         correction is still rising, and closes as the correction relaxes.
%     B   relative error against N in the group detection time
%         Tbar_(1) = int S(t|H=1)^N dt, the quantity the heuristic is used
%         for, with the maximum pointwise survival gap shown for contrast.
%         The pointwise gap is the wrong summary for a steep monotone curve,
%         since a small shift in time gives a large vertical gap wherever the
%         survival is falling fast, and no result in the paper reads S(t)
%         pointwise.
%
%   Printed diagnostics
%     Peak |dS| and its location per N, flagged when the peak falls outside
%     the plotted window; the detection-time error, max|dS| under each
%     hypothesis, and the tail-resolution check (min S1 and the fraction of
%     the window below the resolution floor, since at large N the H = 1
%     survivor pool underflows late and lambda is held at its last resolved
%     value). The caption's quoted figures come from this printout.
%
%   Outputs
%     figs2A.pdf/.png, figs2B.pdf/.png in OUTDIR, and OUT.A / OUT.B holding
%     the plotted series. Both panels are regenerated on every call, so a
%     composite assembled from stale panels will disagree with OUT.
%
%   Requires  solve_bayes, calib_pool, linf_ceiling
%
%   ZPK 2026

if nargin < 2 || isempty(outdir), outdir = 'output'; end
if nargin < 1 || isempty(bayes)
    fprintf('  no solver output passed -- solving the fixed point here\n');
    th1 = calib_pool('theta1', calib_pool('invert', 3.93));
    bayes = solve_bayes([2 5 10 20 40 100], th1);
end

P.outdir = outdir;
P.LW  = 3.5;   P.FSZ = 26;   P.FIG = [22 19];
P.xmax = 8;                 % plotted window in rescaled time
P.Sfloor = 1e-8;            % below this the H=1 hazard is at quadrature noise
P.greens = [0.75 0.90 0.78; 0.55 0.80 0.60; 0.36 0.68 0.44;
            0.20 0.55 0.32; 0.10 0.42 0.24; 0.03 0.28 0.15];
P.grey = [0.45 0.45 0.45];

fprintf('  run at N ='); fprintf(' %d', [bayes.N]); fprintf('\n');
out.A = panelA(bayes, P);
out.B = panelB(bayes, P);
fprintf('\n');
end

%% ---------- A: pointwise discrepancy under H = 1 ------------------------
function A = panelA(bayes, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
n = numel(bayes);
shade = shades(n, P.greens);
D = cell(1,n); Ns = nan(1,n); pk = nan(1,n); tpk = nan(1,n); inwin = false(1,n);
for i = 1:n
    b = bayes(i);  a = linf_ceiling(b.N);  t = b.t(:);  x = t/b.thetaN;
    d = b.S1(:) - surv(t, 1-a, b.thetaN);
    plot(ax, x, d, '-', 'Color', shade(i,:), 'LineWidth', P.LW);
    [pk(i), k] = max(abs(d));  tpk(i) = x(k);
    inwin(i) = tpk(i) < P.xmax;      % is the peak inside the plotted window
    Ns(i) = b.N;  D{i} = [x, d];
end
plot(ax, [0 P.xmax], [0 0], ':', 'Color', P.grey, 'LineWidth', 2);
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'rescaled time $t/\theta_N$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'$S_{\mathrm{Bayes}} - S_{\mathrm{heur}}$', ...
       'Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[0 P.xmax]);
fprintf('  panel A: peak |dS| under H=1  '); fprintf(' %8.3f', pk);  fprintf('\n');
fprintf('           at t/theta_N =       '); fprintf(' %8.2f', tpk); fprintf('\n');
if ~all(inwin)
    fprintf('           NOTE: peak outside the window at N =');
    fprintf(' %d', Ns(~inwin));
    fprintf(' -- those are lower bounds, do not quote them\n');
end
%  POST: label the curves by N (light to dark). H = 0 is omitted, flat at
%        ~1.5e-3 for every N; state that in the caption rather than plotting
%        a null result. The y-axis is left on auto so no curve is clipped.
export(f,'figS2A',P);
A = struct('N',Ns,'diff',{D},'peak',pk,'t_peak',tpk,'peak_in_window',inwin);
end

%% ---------- B: error against group size ---------------------------------
function B = panelB(bayes, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
n = numel(bayes);
Ns = nan(1,n); e1 = nan(1,n); e0 = nan(1,n); eT = nan(1,n);
Smin = nan(1,n); frac = nan(1,n);
for i = 1:n
    b = bayes(i); a = linf_ceiling(b.N); t = b.t(:);
    Sh1 = surv(t,  1-a, b.thetaN);
    Sh0 = surv(t, -1-a, b.thetaN);
    TB = trapz(t, b.S1(:).^b.N);
    Ns(i) = b.N;
    e1(i) = max(abs(b.S1(:) - Sh1));
    e0(i) = max(abs(b.S0(:) - Sh0));
    eT(i) = abs(trapz(t, Sh1.^b.N) - TB)/TB;
    Smin(i) = min(b.S1(:));
    frac(i) = mean(b.S1(:) < P.Sfloor);
end
plot(ax, Ns, eT, '-o', 'Color', P.greens(5,:), 'MarkerFaceColor', P.greens(5,:), ...
     'LineWidth', P.LW, 'MarkerSize', 11);
plot(ax, Ns, e1, '--s', 'Color', P.grey, 'MarkerFaceColor','w', ...
     'LineWidth', P.LW-1.2, 'MarkerSize', 9);
set(ax,'XScale','log','YScale','log','TickDir','out', ...
       'TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'group size $N$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'relative error','Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[min(Ns)/1.3 max(Ns)*1.3]);
fprintf('  panel B: N                    '); fprintf(' %8d',   Ns);   fprintf('\n');
fprintf('           detection-time error  '); fprintf(' %8.2e', eT);   fprintf('\n');
fprintf('           max|dS| H=1           '); fprintf(' %8.2e', e1);   fprintf('\n');
fprintf('           max|dS| H=0           '); fprintf(' %8.2e', e0);   fprintf('\n');
fprintf('           min S1                '); fprintf(' %8.1e', Smin); fprintf('\n');
fprintf('           frac below floor      '); fprintf(' %8.2f', frac); fprintf('\n');
if any(frac > 0.5)
    fprintf('           WARNING: unresolved window at N =');
    fprintf(' %d', Ns(frac > 0.5));
    fprintf(' -- restate the caption at the largest resolved N\n');
end
%  POST: SOLID green circles are the detection-time error, the quantity the
%        heuristic is used for. DASHED grey squares are the max pointwise
%        survival gap under H=1, shown for contrast; its growth with N is
%        expected, since the constant is applied through the transient.
export(f,'figS2B',P);
B = struct('N',Ns,'err_T1',eT,'err_H1',e1,'err_H0',e0, ...
           'min_S1',Smin,'frac_unresolved',frac);
end

%% ---------- helpers -----------------------------------------------------
function C = shades(n, pal)
%SHADES  One distinct shade per curve, light to dark, for any n.
if n < 2, C = pal(end,:); return; end
C = pal(round(linspace(1, size(pal,1), n)), :);
end

function S = surv(t, nu, th)
%SURV  Inverse-Gaussian survival at drift NU, threshold TH, sigma^2 = 2.
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
