%% MAKE_FIGS3  Closed-form validation of the saturation ceiling and hump scaling.
%
%  A  1 - L_inf(N) vs N-1 on log-log axes. Circles are the exact closed form
%     (Eq. S18); the dashed line is the leading asymptote 2/sqrt(N-1), a line
%     of slope -1/2. The +2/k correction lifts the small-N points above it.
%  B  early-time hump height max_t (N-1)lambda^(N)(t) vs N. Exact matched-
%     asymptotics form (N-1)G(beta_N)/theta_N^2 with beta_N = (1-L_inf)theta_N;
%     the dashed line is the driftless leading term with C0 = 1.0207. Unlike
%     the tail, the hump DIVERGES, as N/(log N)^2.
%
%  Fully closed-form: no solver, no data file.  ZPK 2026

clear; close all; st = paper_style();  outdir = 'output';
theta1 = 3.5678;
col_main = [0.10 0.45 0.25];  col_ref = [0.45 0.75 0.55];

%% ---- A ------------------------------------------------------------------
Ns = unique(round(logspace(log10(2), log10(2000), 30)));
kk = Ns - 1;
one_minus = 1 - linf_ceiling(Ns)';

fA = figure('Units','centimeters','Position',[1 1 st.FIG],'Color','w');
axA = axes(fA); hold(axA,'on');
loglog(axA, kk, 2./sqrt(kk), '--', 'Color',col_ref, 'LineWidth',st.LW);
loglog(axA, kk, one_minus, 'o', 'Color',col_main, 'MarkerFaceColor','w', ...
       'MarkerSize',9, 'LineWidth',2.5);
set(axA,'XScale','log','YScale','log');
apply_axes(axA, '$N-1$', '$1 - L_\infty(N)$');
xlim(axA,[1 2200]);  ylim(axA,[0.03 1.1]);
text(axA, 30, 0.42, '$\sim 2/\sqrt{N\!-\!1}$', 'Interpreter','latex', ...
     'FontSize',st.FSZ-8,'Color',col_ref);
export_panel(fA, 'figS3A', outdir);

%% ---- B ------------------------------------------------------------------
% C0 = max_u Phi(u), u* the root of dlog(Phi)/du = 0  (Eqs. S20-S22)
dlogPhi = @(u) 1.5./u - 1 - (exp(-u)./sqrt(u)/sqrt(pi))./erf(sqrt(u));
ustar   = fzero(dlogPhi, 1.3);
Phi     = @(u) (4*u).^1.5 .* exp(-u) ./ (sqrt(4*pi).*erf(sqrt(u)));
C0      = Phi(ustar);

Nh   = [2 5 10 20 40 100 200 400 1000 2000];
Lh   = linf_ceiling(Nh);  Lh = Lh(:).';        % row, to match Nh
thN  = (theta1 + log(Nh))./(1 + Lh);           % Eq. (threshold_scaling) at alpha = L_inf(N)
muT  = 1 - Lh;

hump_exact = arrayfun(@(i) (Nh(i)-1)*peak_hazard(muT(i), thN(i))/thN(i)^2, 1:numel(Nh));
hump_drift = (Nh-1).*C0./thN.^2;

fB = figure('Units','centimeters','Position',[3 1 st.FIG],'Color','w');
axB = axes(fB); hold(axB,'on');
loglog(axB, Nh, hump_drift, '--', 'Color',col_ref, 'LineWidth',st.LW);
loglog(axB, Nh, hump_exact, '-',  'Color',col_main, 'LineWidth',st.LWt);
loglog(axB, Nh, hump_exact, 'o',  'Color',col_main, 'MarkerFaceColor','w', ...
       'MarkerSize',9, 'LineWidth',2.5);
plot(axB, [Nh(1) Nh(end)], [1 1], ':', 'Color',st.col_ref, 'LineWidth',2);
set(axB,'XScale','log','YScale','log');
apply_axes(axB, '$N$', '$\max_t\,(N\!-\!1)\lambda^{(N)}(t)$');
xlim(axB,[1.7 2200]);  ylim(axB,[0.2 150]);
text(axB, 3, 1.35, 'tail ceiling $=1$', 'Interpreter','latex', ...
     'FontSize',st.FSZ-12,'Color',st.col_ref);
export_panel(fB, 'figS3B', outdir);

fprintf('Fig S1: C0 = %.4f at u* = %.4f  (paper quotes 1.0207, 1.30437)\n', C0, ustar);
fprintf('  L_inf(2) = %.6f  vs  3-2*sqrt(2) = %.6f\n', linf_ceiling(2), 3-2*sqrt(2));

%% ---- local -------------------------------------------------------------
function G = peak_hazard(nu, th)
% theta^2 * peak first-passage hazard at drift nu, threshold th, sigma^2 = 2.
% Universal in beta = nu*th; this is G(beta) of Eq. (S24).
    t = linspace(1e-4, 30*th, 4e5)';
    o = ddm_ig(t, nu, th);
    G = th^2 * max(o.h);
end

