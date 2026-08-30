%% MAKE_FIGS1  Single-agent first-passage statistics and the speed-accuracy
%             frontier.
%
%  A  first-passage densities f(t|H): the H=1 detection density and the
%     defective H=0 false-alarm density (area = q = exp(-theta)).
%  B  the single-agent frontier (T,q) = (theta, exp(-theta)).
%
%  Both panels were Fig. 2A and 2C in the Science Advances version; they moved
%  to the SI when the main text went to four figures. The third panel of that
%  old figure, the survival LLR, is now Fig. 1D.
%
%  Fully closed-form; no solver, no .mat.  ZPK 2026

clear; close all; st = paper_style();  outdir = 'output';

%% ---- A: first-passage densities ----------------------------------------
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
export_panel(fA, 'figS1A', outdir);

%% ---- B: frontier -------------------------------------------------------
theta_F = linspace(0.2, 8.5, 400)';
fB = figure('Units','centimeters','Position',[5 1 st.FIG],'Color','w');
axB = axes(fC); hold(axB,'on');
plot(axB, theta_F, exp(-theta_F), '-', 'Color',[0.1 0.1 0.1], 'LineWidth',st.LW);
set(axB,'YScale','log');
apply_axes(axB, 'detection time $\bar{T}$', 'false alarm probability $q$');
xlim(axB,[0 8.5]);  ylim(axB,[1e-4 1]);
export_panel(fB, 'figS1B', outdir);

fprintf('Fig S1: theta=%.2f -> q=%.4f; frontier (T,q)=(theta, e^{-theta})\n', theta, q);

