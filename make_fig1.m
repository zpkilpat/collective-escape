%% MAKE_FIG1  Combining private evidence and social information.
%
%  Panel A (bird/shoal cartoon) is drawn separately in Illustrator.
%
%  B  single-agent DDM sample paths: one H=1 realisation absorbed at theta,
%     one H=0 realisation drifting away.
%  C  the mechanism, naive vs discounting on shared axes. Both beliefs are the
%     SAME noise realisation under a real threat; they differ only by the
%     discounting rate alpha, so the vertical gap is exactly alpha*t. One
%     neighbour departs at T1 and both receive the SAME kick: the naive belief
%     is carried across theta (false alarm), the discounting belief is not.
%
%  Self-contained: no .mat dependency.  ZPK 2026

clear; close all; st = paper_style();  outdir = 'output';
theta = 3.5678;  dt = 0.01;

%% ---- B: sample paths ---------------------------------------------------
rng(3);
[tb1, xb1] = simpath(+1, theta, dt, 40);
[tb0, xb0] = simpath(-1, theta, dt,  8);

fB = figure('Units','centimeters','Position',[1 1 st.FIG],'Color','w');
axB = axes(fB); hold(axB,'on');
yline(axB, theta, '--', 'Color',st.col_ref, 'LineWidth',st.LWt);
plot(axB, tb0, xb0, '-', 'Color',st.col_H0, 'LineWidth',st.LW);
plot(axB, tb1, xb1, '-', 'Color',st.col_H1, 'LineWidth',st.LW);
plot(axB, tb1(end), theta, 'o', 'MarkerSize',12, 'MarkerFaceColor',st.col_H1, ...
     'MarkerEdgeColor','k', 'LineWidth',1.5);
text(axB, 0.3, theta+0.4, '$\theta$', 'Interpreter','latex', ...
     'FontSize',st.FSZ, 'Color',st.col_ref);
apply_axes(axB, 'time $t$', 'private belief $\xi(t)$');
xlim(axB,[0 max(tb1(end),tb0(end))*1.05]);  ylim(axB,[-4 theta+1]);
export_panel(fB, 'fig1B', outdir);

%% ---- C: naive vs discounting ------------------------------------------
xmax = 4.0;  T1 = 2.0;  stub = 0.5;
alpha_naive = 0.0;  alpha_disc = 0.8;
t2 = (0:dt:xmax)';

rng(5);
base = cumsum(1.0*dt + sqrt(dt)*0.65*randn(numel(t2),1));  base = base - base(1);
priv = base - alpha_naive*t2;
bel  = base - alpha_disc *t2;

k      = find(t2 >= T1, 1);
xpre_n = priv(k);  xpre_b = bel(k);
kappa  = theta - xpre_n + 0.18;        % shared kick; naive just clears theta

fC = figure('Units','centimeters','Position',[3 1 st.FIG],'Color','w');
axC = axes(fC); hold(axC,'on');
yline(axC, theta, '--', 'Color',st.col_ref, 'LineWidth',st.LWt);
fill(axC, [t2(1:k); flipud(t2(1:k))], [bel(1:k); flipud(priv(1:k))], ...
     [0.55 0.72 0.55], 'FaceAlpha',0.40, 'EdgeColor','none');
plot(axC, t2(1:k), priv(1:k), '-', 'Color',st.col_naiv, 'LineWidth',st.LW);
plot(axC, t2(1:k), bel(1:k),  '-', 'Color',st.col_bay,  'LineWidth',st.LW);
draw_kick(axC, T1, xpre_n, xpre_n+kappa, dt, stub, st.col_naiv, st.col_ref, st.LW, true);
draw_kick(axC, T1, xpre_b, xpre_b+kappa, dt, stub, st.col_bay,  st.col_ref, st.LW, false);

text(axC, T1+stub+0.08, xpre_n+kappa, 'naive ($\alpha{=}0$) crosses $\theta$', ...
     'Interpreter','latex','FontSize',st.FSZ-10,'Color',st.col_naiv,'VerticalAlignment','middle');
text(axC, T1+stub+0.08, xpre_b+kappa, 'discounting ($\alpha{=}0.8$) stays below', ...
     'Interpreter','latex','FontSize',st.FSZ-10,'Color',st.col_bay,'VerticalAlignment','middle');
mid = find(t2 >= 1.35, 1);
text(axC, 0.45, 2.45, 'discounted silence', 'Interpreter','latex', ...
     'FontSize',st.FSZ-11,'Color',st.col_bay*0.8);
text(axC, 0.45, 2.10, '($-\alpha t$)', 'Interpreter','latex', ...
     'FontSize',st.FSZ-11,'Color',st.col_bay*0.8);
plot(axC, [0.95 1.30], [2.05 0.5*(priv(mid)+bel(mid))], '-', ...
     'Color',st.col_bay*0.8, 'LineWidth',1.4);
text(axC, 0.2, theta+0.4, '$\theta$', 'Interpreter','latex','FontSize',st.FSZ,'Color',st.col_ref);
text(axC, T1-0.55, -0.62, 'neighbor departs', 'Interpreter','latex', ...
     'FontSize',st.FSZ-11,'HorizontalAlignment','center','Color',st.col_ref);
apply_axes(axC, 'time $t$', 'private belief $\xi_2(t)$');
xlim(axC,[0 xmax]);  ylim(axC,[-0.8 theta+1.6]);
export_panel(fC, 'fig1C', outdir);

fprintf('Fig 1: naive clears theta by %.2f; discounting short by %.2f\n', ...
        xpre_n+kappa-theta, theta-(xpre_b+kappa));

%% ---- local functions ---------------------------------------------------
function [t, x] = simpath(mu, theta, dt, Tmax)
    n = round(Tmax/dt);  x = zeros(n,1);  t = (0:n-1)'*dt;
    for i = 2:n
        x(i) = x(i-1) + mu*dt + sqrt(2*dt)*randn;   % sigma^2 = 2
        if x(i) >= theta, x(i) = theta; x = x(1:i); t = t(1:i); return; end
    end
end

function draw_kick(ax, Tj, xpre, xtop, dt, stub, col, colm, LW, cross_marker)
    plot(ax, [Tj Tj], [xpre xtop], '-', 'Color',col, 'LineWidth',LW);
    ts = (Tj:dt:Tj+stub)';
    plot(ax, ts, xtop*ones(size(ts)), '-', 'Color',col, 'LineWidth',LW);
    plot(ax, Tj, xpre, 'v', 'MarkerSize',11, 'MarkerFaceColor',colm, ...
         'MarkerEdgeColor','k', 'LineWidth',1.1);
    if cross_marker
        plot(ax, Tj, xtop, 'o', 'MarkerSize',11, 'MarkerFaceColor',col, ...
             'MarkerEdgeColor','k', 'LineWidth',1.4);
    end
end
