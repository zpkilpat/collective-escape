%% MAKE_FIG2  Evidence accumulation and the single-agent speed-accuracy tradeoff.
%
%  A  first-passage densities f(t|H): the H=1 detection density and the
%     defective H=0 false-alarm density (area = q = exp(-theta)).
%  B  survival log-likelihood ratio Lambda^surv(t) for theta = 1,2,3.
%  C  the single-agent frontier (T,q) = (theta, exp(-theta)).
%
%  Fully closed-form; no solver, no .mat.  ZPK 2026

clear; close all; st = paper_style();  outdir = 'output';

%% ---- A ------------------------------------------------------------------
theta = 2.0;
t  = linspace(0.005, 8.0, 3000)';
f1 = ddm_ig(t,  1, theta).f;
f0 = ddm_ig(t, -1, theta).f;
q  = exp(-theta);

fA = figure('Units','centimeters','Position',[1 1 st.FIG],'Color','w');
axA = axes(fA); hold(axA,'on');
fill(axA, [t; flipud(t)], [f1; zeros(size(f1))], st.col_H1, 'FaceAlpha',0.16, 'EdgeColor','none');
plot(axA, t, f1, '-', 'Color',st.col_H1, 'LineWidth',st.LW);
fill(axA, [t; flipud(t)], [f0; zeros(size(f0))], st.col_H1, 'FaceAlpha',0.85, 'EdgeColor','none');
plot(axA, t, f0, '-', 'Color',st.col_H1, 'LineWidth',st.LW);
text(axA, 2.6, max(f1)*0.62, {'true detection','$f(t\,|\,H{=}1)$'}, ...
     'Interpreter','latex','FontSize',st.FSZ-8,'Color',st.col_H1);
text(axA, 4.0, max(f0)*1.8, {'false alarm','$f(t\,|\,H{=}0)$'}, ...
     'Interpreter','latex','FontSize',st.FSZ-8,'Color',st.col_H1);
text(axA, 4.0, max(f0)*0.9, sprintf('area $= q = %.2f$', q), ...
     'Interpreter','latex','FontSize',st.FSZ-12,'Color',st.col_H1);
apply_axes(axA, 'decision time $t$', 'first-passage density $f(t\,|\,H)$');
xlim(axA,[0 8]);  ylim(axA,[0 max(f1)*1.15]);
export_panel(fA, 'fig2A', outdir);

%% ---- B ------------------------------------------------------------------
theta_vals = [1 2 3];
shades = [0.70 0.70 0.70; 0.40 0.40 0.40; 0 0 0];
tb = (0.01:0.005:6)';

fB = figure('Units','centimeters','Position',[3 1 st.FIG],'Color','w');
axB = axes(fB); hold(axB,'on');
yline(axB, 0, ':', 'Color',st.col_ref, 'LineWidth',st.LWt);
for k = 1:numel(theta_vals)
    S1 = ddm_ig(tb,  1, theta_vals(k)).S;
    S0 = ddm_ig(tb, -1, theta_vals(k)).S;
    plot(axB, tb, log(max(S1,1e-12)./max(S0,1e-12)), '-', ...
         'Color',shades(k,:), 'LineWidth',st.LW);
end
apply_axes(axB, 'time $t$', 'survival LLR $\Lambda^{\rm surv}(t)$');
xlim(axB,[0 6]);  ylim(axB,[-8 0.5]);
export_panel(fB, 'fig2B', outdir);

%% ---- C ------------------------------------------------------------------
theta_F = linspace(0.2, 8.5, 400)';
fC = figure('Units','centimeters','Position',[5 1 st.FIG],'Color','w');
axC = axes(fC); hold(axC,'on');
plot(axC, theta_F, exp(-theta_F), '-', 'Color',[0.1 0.1 0.1], 'LineWidth',st.LW);
set(axC,'YScale','log');
apply_axes(axC, 'detection time $\bar{T}$', 'false alarm probability $q$');
xlim(axC,[0 8.5]);  ylim(axC,[1e-4 1]);
export_panel(fC, 'fig2C', outdir);

fprintf('Fig 2: theta=%.2f -> q=%.4f; frontier (T,q)=(theta, e^{-theta})\n', theta, q);
