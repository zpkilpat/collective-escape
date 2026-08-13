function out = make_figS5(csvfile, outdir)
%MAKE_FIGS5  Rate estimation, clustering, and what the intervals rest on.
%
%  A  response rate against shoal area for both event types, logistic fits
%     with CI ribbons, plus the interaction test. Attacks rise with area and
%     flybys are flat, but the interaction does not reach significance, so
%     the honest statement is that the two slopes are not resolved apart.
%  B  iid against file|bout cluster bootstrap intervals for every rate-derived
%     quantity. The clustered intervals are the ones to report.
%  C  the cluster structure itself. bout_id restarts within each recording,
%     so the key must be file|bout rather than bout alone, with singleton
%     labels for rows carrying no bout. 73 clusters, 47 carrying attacks,
%     26 carrying flybys.
%
%     The 83/47/17 triple quoted in earlier drafts came from keying on
%     location|bout instead of file|bout, a column-index error. So did the
%     claim that bout_id is NA on every flyby row: flybys carry their own
%     labels, so the singleton branch almost never fires. Stratifying the
%     resample by event type is still required, but for the ordinary reason
%     that an unstratified draw can empty the flyby stratum.
%
%     Note the coincidence: there are 26 flyby-carrying clusters AND 26
%     flyby responses out of 81. Unrelated quantities, same number.
%
%  Usage:  out = make_figS5('sciadv_adt8600_data_s1.csv', 'output');
%
%  ZPK 2026

if nargin < 1 || isempty(csvfile), csvfile = 'sciadv_adt8600_data_s1.csv'; end
if nargin < 2 || isempty(outdir),  outdir  = 'output'; end
P.outdir = outdir;
P.LW = 3.5;  P.FSZ = 26;  P.FIG = [22 19];
P.green = [0.20 0.55 0.32];  P.blue = [0.08 0.24 0.48];
P.grey  = [0.45 0.45 0.45];
P.B = 4000;
% SAME seed as make_fig5's bootstrap_joint and make_figS4 (opts.seed =
% 20260811). boot(...,true) below is structurally identical to that routine --
% same stratified cluster draw, same tt selection, same alpha_hat -- so with a
% matched seed the clustered alpha interval here reproduces the one quoted in
% the main text EXACTLY rather than approximately. It is the same interval
% drawn in three places, so it must come from one stream.
P.seed = 20260811;

rng(P.seed);
D = read_pacher(csvfile);
assert_key(D);
out.A = panelA(D, P);
out.B = panelB(D, P);
out.C = panelC(D, P);
fprintf('\n');
end

%% ============== A: rates against area, with the interaction ============
function A = panelA(D, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
xg = linspace(min(D.area(~isnan(D.area))), max(D.area(~isnan(D.area))), 200)';
[bA, sA] = fitone(D, true,  xg, ax, P.green, 'o', P);
[bF, sF] = fitone(D, false, xg, ax, P.blue,  's', P);
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'shoal area $A$ (m$^2$)','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'response probability','Interpreter','latex','FontSize',P.FSZ+8);
ylim(ax,[0 1]);

% pooled model with an event-type interaction
sel = ~isnan(D.resp) & ~isnan(D.area);
a = D.area(sel); y = double(D.resp(sel)==1); g = double(D.attack(sel));
X = [ones(sum(sel),1), a, g, a.*g];
[b, se] = irls_logit(X, y);
z = b(4)/se(4); p = 2*(1-normcdf_local(abs(z)));
fprintf('  panel A: attacks slope %+.4f (SE %.4f, P = %.3f)\n', bA(2), sA(2), ...
        2*(1-normcdf_local(abs(bA(2)/sA(2)))));
fprintf('           flybys  slope %+.4f (SE %.4f, P = %.3f)\n', bF(2), sF(2), ...
        2*(1-normcdf_local(abs(bF(2)/sF(2)))));
fprintf('           interaction %.4f, z = %.2f, P = %.3f\n', b(4), z, p);
%  POST: green circles attacks, blue squares flybys, quintile means with
%        Wilson intervals, sized by bin count; ribbons are 95% CI on the
%        logistic fit. State the interaction test in the caption and say
%        plainly that the slopes are not resolved apart.
export(f,'figs5A',P);
A = struct('beta_attack',bA,'se_attack',sA,'beta_flyby',bF,'se_flyby',sF, ...
           'interaction',b(4),'z',z,'p',p);
end

%% ============== B: iid against clustered intervals =====================
function B = panelB(D, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
% Reseed before EACH resample. Without this the clustered draw starts from a
% stream already advanced by 4000 iid replicates and cannot match make_fig5,
% which reseeds immediately before its own clustered draw.
rng(P.seed);  Ri = boot(D, P.B, false);   % iid resample of events
rng(P.seed);  Rc = boot(D, P.B, true);    % file|bout cluster resample, stratified
% The area slopes belong here. Sec. scaling says of the logistic standard
% errors that "the clustered versions are reported in fig. S5B" -- before this
% they were not: the panel carried only the rate-derived quantities. Slopes are
% O(0.03), so they are shown x10 to share the axis, as theta_1 and K_max are
% shown rescaled.
names = {'$q$','$\mathrm{TP}$','$\hat\alpha$','$q_1$','$\theta_1/10$', ...
         '$K_{\max}/100$','$10\beta_{\rm atk}$','$10\beta_{\rm fly}$'};
fi = {@(r) r.q, @(r) r.TP, @(r) r.alpha, @(r) r.q1, @(r) r.th1/10, ...
      @(r) r.Kmax/100, @(r) 10*r.slopeA, @(r) 10*r.slopeF};
for k = 1:numel(fi)
    vi = fi{k}(Ri); vc = fi{k}(Rc);
    vi = vi(~isnan(vi)); vc = vc(~isnan(vc));
    ci_i = prctile(vi,[2.5 97.5]);  ci_c = prctile(vc,[2.5 97.5]);
    plot(ax, k-0.14+[0 0], ci_i, '-', 'Color', P.grey,  'LineWidth', P.LW);
    plot(ax, k+0.14+[0 0], ci_c, '-', 'Color', P.green, 'LineWidth', P.LW);
    plot(ax, k-0.14, median(vi), 'o', 'MarkerFaceColor', P.grey,  'MarkerEdgeColor', P.grey,  'MarkerSize', 9);
    plot(ax, k+0.14, median(vc), 'o', 'MarkerFaceColor', P.green, 'MarkerEdgeColor', P.green, 'MarkerSize', 9);
    fprintf('  %-14s iid [%.3f, %.3f]   clustered [%.3f, %.3f]  (x%.2f wider)\n', ...
            names{k}, ci_i(1), ci_i(2), ci_c(1), ci_c(2), diff(ci_c)/diff(ci_i));
end
set(ax,'XTick',1:numel(fi),'XTickLabel',names,'TickLabelInterpreter','latex', ...
       'TickDir','out','FontSize',P.FSZ,'LineWidth',2);
ylabel(ax,'estimate and 95\% interval','Interpreter','latex','FontSize',P.FSZ+8);
xlim(ax,[0.5 numel(fi)+0.5]);
%  POST: grey is the iid event resample, green the file|bout cluster
%        resample stratified by event type. Report the clustered intervals
%        throughout; note the rescaling of theta_1 and K_max on the axis.
export(f,'figs5B',P);
B = struct('iid',Ri,'clustered',Rc);
end

%% ============== C: the cluster structure ===============================
function C = panelC(D, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
[keys, ~, ci] = unique(D.clust);
nc = numel(keys);
nA = accumarray(ci, double(D.attack), [nc 1]);
nF = accumarray(ci, double(~D.attack), [nc 1]);
[~, ord] = sort(nA + nF, 'descend');
barh(ax, 1:nc, nA(ord), 'FaceColor', P.green, 'EdgeColor','none');
barh(ax, 1:nc, -nF(ord), 'FaceColor', P.blue,  'EdgeColor','none');
plot(ax, [0 0], [0 nc+1], '-', 'Color','w', 'LineWidth', 1.5);
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'events per cluster','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax,'recording $\times$ bout cluster','Interpreter','latex','FontSize',P.FSZ+8);
ylim(ax,[0 nc+1]);
fprintf('  panel C: %d clusters, %d carrying attacks, %d carrying flybys\n', ...
        nc, sum(nA>0), sum(nF>0));
%  POST: green to the right is attacks, blue to the left flybys. The
%        twenty-six flyby-carrying clusters are what every q-derived
%        interval rests on, which is the reason the clustered intervals in
%        panel B are wider than the raw count of 81 flybys suggests. Read
%        the count off the console line rather than transcribing it here.
export(f,'figs5C',P);
C = struct('nclust',nc,'nAttackClust',sum(nA>0),'nFlybyClust',sum(nF>0));
end

%% ============================ helpers ==================================
function R = boot(D, B, clustered)
[~, ~, ci] = unique(D.clust);
cA = unique(ci(D.attack)); cF = unique(ci(~D.attack));
iA = find(D.attack & ~isnan(D.resp)); iF = find(~D.attack & ~isnan(D.resp));
R.q = nan(B,1); R.TP = nan(B,1); R.alpha = nan(B,1);
R.q1 = nan(B,1); R.th1 = nan(B,1); R.Kmax = nan(B,1);
R.slopeA = nan(B,1); R.slopeF = nan(B,1);
for b = 1:B
    if clustered
        pick = [cA(randi(numel(cA),numel(cA),1)); cF(randi(numel(cF),numel(cF),1))];
        idx = [];
        for j = 1:numel(pick), idx = [idx; find(ci == pick(j))]; end %#ok<AGROW>
    else
        idx = [iA(randi(numel(iA),numel(iA),1)); iF(randi(numel(iF),numel(iF),1))];
    end
    a = D.attack(idx); r = D.resp(idx); t = D.tfirst(idx);
    okA = a & ~isnan(r); okF = ~a & ~isnan(r);
    tt = t(a & ~isnan(t));            % detection latencies only
    if ~any(okA) || ~any(okF) || numel(tt) < 20, continue; end
    TP = mean(r(okA)==1); q = mean(r(okF)==1);
    if TP <= q, continue; end
    s = (1/(mean(1./tt) - 1/mean(tt)))/mean(tt);
    M = calib_pool('invert', s);
    q1 = 1-(1-q)^(1/M);
    R.q(b) = q; R.TP(b) = TP; R.alpha(b) = alpha_hat(M, q, TP);
    R.q1(b) = q1; R.th1(b) = -log(q1);
    R.Kmax(b) = log(1-q1)/log(1-q1^2);
    % logistic slope in shoal area, per event type, on this same resample
    ar = D.area(idx);
    for e = [1 0]
        m = (a == logical(e)) & ~isnan(r) & ~isnan(ar);
        if sum(m) < 10 || numel(unique(ar(m))) < 5, continue; end
        bb = irls_logit([ones(sum(m),1) ar(m)], double(r(m)==1));
        if e, R.slopeA(b) = bb(2); else, R.slopeF(b) = bb(2); end
    end
end
end

function [b, se] = fitone(D, wantAttack, xg, ax, col, mk, P)
sel = (D.attack == wantAttack) & ~isnan(D.resp) & ~isnan(D.area);
a = D.area(sel); y = double(D.resp(sel)==1);
[b, se, V] = irls_logit([ones(numel(a),1) a], y);
X = [ones(numel(xg),1) xg];
% Prediction variance needs the FULL covariance, x' V x, not just the
% diagonal. With area running 30-75 rather than centred, the intercept and
% slope are strongly negatively correlated, and dropping the off-diagonal
% term inflates the ribbon until it spans the panel.
eta = X*b; s = sqrt(sum((X*V).*X, 2));
lo = 1./(1+exp(-(eta-1.96*s))); hi = 1./(1+exp(-(eta+1.96*s)));
patch(ax, [xg; flipud(xg)], [lo; flipud(hi)], col, 'FaceAlpha',0.18,'EdgeColor','none');
plot(ax, xg, 1./(1+exp(-eta)), '-', 'Color', col, 'LineWidth', P.LW);
q = quantile_local(a, [0 .2 .4 .6 .8 1]);
for j = 1:5
    in = a >= q(j) & a <= q(j+1);
    if sum(in) < 3, continue; end
    k = sum(y(in)); n = sum(in); [l, h] = wilson(k, n);
    plot(ax, [median(a(in)) median(a(in))], [l h], '-', 'Color', col, 'LineWidth', 2);
    plot(ax, median(a(in)), k/n, mk, 'MarkerFaceColor', col, ...
         'MarkerEdgeColor', col, 'MarkerSize', 5+12*n/numel(a));
end
end

function [b, se, V] = irls_logit(X, y)
b = zeros(size(X,2),1);
for it = 1:100
    p = 1./(1+exp(-X*b)); W = max(p.*(1-p),1e-9);
    z = X*b + (y(:)-p)./W;
    bn = (X'*(W.*X)) \ (X'*(W.*z));
    if max(abs(bn-b)) < 1e-10, b = bn; break; end
    b = bn;
end
p = 1./(1+exp(-X*b)); W = max(p.*(1-p),1e-9);
V  = inv(X'*(W.*X));          % full parameter covariance
se = sqrt(diag(V));
end


function q = quantile_local(x, p)
x = sort(x(:)); n = numel(x);
h = (n-1)*p(:) + 1; lo = floor(h); hi = ceil(h);
q = x(lo) + (h - lo).*(x(hi) - x(lo));
end

function p = normcdf_local(z), p = 0.5*(1+erf(z/sqrt(2))); end



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
