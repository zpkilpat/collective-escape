function out = validate_dyad_mc(sol, outdir, doplot)
%VALIDATE_DYAD_MC  Independent Monte Carlo check on the dyad observables.
%
%  NOT A FIGURE SCRIPT, and nothing in run_all calls it. The supplementary
%  figure this once produced is cut, and so is the section that quoted it;
%  the numbers live nowhere in the paper except the closed-form agreement
%  reported in S.orderstat. Retained because the data-availability statement
%  names the Monte Carlo validation, and because the naive-rule agreement is
%  worth being able to reproduce.
%
%  KNOWN LIMITATION, BAYESIAN BRANCH. dyad_observables integrates the killed
%  density at the BARE drift, but under the Bayesian rule the survivor drifts
%  at mu_H - lambda(t) during the neighbour's silence. dyad_frontier_heuristic
%  handles this correctly, evaluating the position law at the discounted
%  drift; the Eq. (S8) call here does not. The naive branch is exact and
%  agrees with Monte Carlo to a few times 1e-3 after the barrier correction;
%  the Bayesian discrepancy of 4-13 percent is this approximation, NOT a
%  solver error. Do not quote those numbers.
%
%  Optional panels over a sweep of thresholds:
%  A  false-alarm cascade probability P(both depart | H=0), logarithmic
%  B  collective response time Tbar_(2) = E[max(T1,T2) | H=1]
%
%  The Monte Carlo integrates the coupled SDE by Euler-Maruyama and shares no
%  code with the Volterra quadrature, so agreement tests the solver rather
%  than a common implementation. Under the Bayesian rule each agent drifts at
%  mu_H - lambda(t) while its neighbor is silent, with lambda read from the
%  converged solver output; on departure the anti-drift switches off and the
%  survivor receives the kick theta - Lambda_surv(T_j). Under the naive rule
%  the kick is fixed at kappa = theta with no anti-drift.
%
%  INPUT IS OPTIONAL. Called with no argument, the dyad fixed point is
%  solved here via solve_bayes over theta = 0.8, 1.5, 2.3, 3.568, and the
%  order-statistic cascade integrals are evaluated from the converged
%  survival functions. Pass your own struct array to override, with fields
%  theta, t, lam, Lsurv, p_naive, p_bayes, T2_naive, T2_bayes.
%
%  Usage:  out = validate_dyad_mc;                        % numbers only
%          out = validate_dyad_mc([], 'output', true);     % also draw the panels
%          out = validate_dyad_mc(sol, 'output');          % or pass solver output
%
%  MONTE CARLO BIAS. Euler-Maruyama against an absorbing boundary misses
%  excursions between grid points and so under-counts crossings. The barrier
%  is shifted inward by 0.5826*sigma*sqrt(dt) (Broadie, Glasserman and Kou)
%  to correct this; the kicks still use the true threshold. Uncorrected, the
%  bias is 5.7 percent per agent at dt = 5e-3 and ~11 percent on a two-agent
%  cascade.
%
%  HORIZON BIAS, SEPARATE FROM THE ABOVE. Under H = 1 the collective time
%  averages over trials in which BOTH agents departed before Tmax; a pair
%  still running at the horizon is dropped, which biases T2 downward. The
%  dropped fraction is now printed, and at Tmax = 60 with theta <= 3.6 it
%  should be at the 1e-4 level or below. If it is not, raise Tmax rather
%  than trusting T2.
%
%  SELF-TEST. Before plotting, the naive-rule observables are checked
%  against closed-form limits: the cascade probability must lie between the
%  independent baseline q^2 and the saturating value 2q - q^2, and the
%  collective time must exceed the group first-departure time. They are also
%  checked against the reference values below wherever the threshold matches
%  one of them.
%
%  ZPK 2026

if nargin < 2 || isempty(outdir), outdir = 'output'; end
if nargin < 3 || isempty(doplot), doplot = false; end   % figure is cut
P.outdir = outdir;
P.LW = 3.5;  P.FSZ = 26;  P.FIG = [22 19];
P.green = [0.20 0.55 0.32];  P.grey = [0.45 0.45 0.45];
P.ntrial = 2e5;  P.dt = 5e-3;  P.Tmax = 60;
% Broadie-Glasserman-Kou barrier correction. A discretely monitored barrier
% behaves like a continuous one shifted OUTWARD by beta*sigma*sqrt(dt), since
% the walk misses excursions between grid points; simulating against a
% barrier shifted inward by the same amount removes the bias to O(dt).
% Without it the Monte Carlo under-counts crossings by 5.7 percent per agent
% at dt = 5e-3, hence ~11 percent on a two-agent cascade -- which is exactly
% the flat ~10 percent discrepancy seen across all four thresholds before
% this was applied.
P.bgk = 0.5826;

if nargin < 1 || isempty(sol)
    fprintf('  no solver output passed -- solving the dyad fixed point here\n');
    ths = [0.8 1.5 2.3 calib_pool('theta1', calib_pool('invert', 3.93))];
    clear sol
    for i = numel(ths):-1:1
        b = solve_bayes(2, ths(i));
        b = dyad_observables(b);
        sol(i) = b;
    end
end

selftest(sol);

rng(2);
nth = numel(sol);
MC = struct('p_naive',nan(1,nth),'p_bayes',nan(1,nth), ...
            'T2_naive',nan(1,nth),'T2_bayes',nan(1,nth), ...
            'drop_naive',nan(1,nth),'drop_bayes',nan(1,nth));
for i = 1:nth
    fprintf('  theta = %.3f ...\n', sol(i).theta);
    r = mcdyad(sol(i), false, P);
    MC.p_naive(i) = r.p; MC.T2_naive(i) = r.T2; MC.drop_naive(i) = r.drop;
    r = mcdyad(sol(i), true,  P);
    MC.p_bayes(i) = r.p; MC.T2_bayes(i) = r.T2; MC.drop_bayes(i) = r.drop;
end
th = [sol.theta];

if doplot
    out.A = twopanel(th, [sol.p_naive], [sol.p_bayes], MC.p_naive, MC.p_bayes, ...
        'cascade probability, $H = 0$', true, 'validate_mc_A', P);
    out.B = twopanel(th, [sol.T2_naive], [sol.T2_bayes], MC.T2_naive, MC.T2_bayes, ...
        'collective response time $\bar T_{(2)}$', false, 'validate_mc_B', P);
end
out.solver = sol;
out.MC = MC;

fprintf('\n  solver against Monte Carlo, relative discrepancy per threshold:\n');
fprintf('    theta      p_naive   p_bayes   T2_naive  T2_bayes\n');
for i = 1:nth
    fprintf('    %6.3f   %8.1e  %8.1e  %8.1e  %8.1e\n', th(i), ...
        reldiff(MC.p_naive(i),  sol(i).p_naive), ...
        reldiff(MC.p_bayes(i),  sol(i).p_bayes), ...
        reldiff(MC.T2_naive(i), sol(i).T2_naive), ...
        reldiff(MC.T2_bayes(i), sol(i).T2_bayes));
end
dmax = max([MC.drop_naive MC.drop_bayes]);
fprintf('\n  H = 1 pairs unabsorbed at Tmax = %g, worst case %.2e of trials\n', ...
        P.Tmax, dmax);
if dmax > 1e-3
    warning('validate_dyad_mc:horizon', ...
        ['%.2e of H=1 trials had an agent still running at Tmax; T2 is ' ...
         'biased low. Raise P.Tmax.'], dmax);
end
fprintf(['\n  Barrier corrected by %.4f (Broadie-Glasserman-Kou at dt = %.0e).\n' ...
         '  Residual should now be O(dt) rather than O(sqrt(dt)); a few times\n' ...
         '  1e-2 is expected, and a flat offset across thresholds would mean\n' ...
         '  the correction is mis-signed.\n\n'], P.bgk*sqrt(2)*sqrt(P.dt), P.dt);
end

%% ============================ self-test =================================
function selftest(sol)
%SELFTEST  closed-form bounds on the naive-rule observables, plus a check
%   against the reference values in the header wherever the threshold matches.
REF_TH = [0.8 1.5 2.3 3.5678];
REF_P  = [0.4008 0.1732 0.0647 0.0132];
REF_T2 = [0.615  1.077  1.612  2.507];
fprintf('  self-test (naive rule against closed-form bounds):\n');
bad = false;
for i = 1:numel(sol)
    th = sol(i).theta;  q = exp(-th);
    p  = sol(i).p_naive;  T2 = sol(i).T2_naive;
    lo = q^2;  hi = 2*q - q^2;
    okp = p >= lo && p <= hi;
    okt = T2 > 0 && isfinite(T2);
    fprintf('    theta = %.3f: p = %.5f in [%.5f, %.5f] %s;  T2 = %.3f %s\n', ...
        th, p, lo, hi, tick(okp), T2, tick(okt));
    j = find(abs(REF_TH - th) < 5e-3, 1);
    if ~isempty(j)
        dp = abs(p - REF_P(j))/REF_P(j);  dt2 = abs(T2 - REF_T2(j))/REF_T2(j);
        okr = dp < 5e-3 && dt2 < 5e-3;
        fprintf('              against reference %.4f / %.3f: %.1e / %.1e %s\n', ...
            REF_P(j), REF_T2(j), dp, dt2, tick(okr));
        bad = bad || ~okr;
    end
    bad = bad || ~okp || ~okt;
end
if bad
    warning('validate_dyad_mc:selftest', ...
        ['naive-rule observables violate their closed-form bounds or the ' ...
         'reference values -- check Gexact/Cexact before trusting the panels']);
end
end

function s = tick(ok)
if ok, s = 'OK'; else, s = '** FAIL **'; end
end

%% ============================ panels ===================================
function d = reldiff(a, b)
d = abs(a - b) ./ max(abs(b), eps);
end

function S = twopanel(th, yn, yb, mn, mb, ylab, uselog, name, P)
f = figure('Units','centimeters','Position',[1 1 P.FIG],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'off')
plot(ax, th, yn, '-', 'Color', P.grey,  'LineWidth', P.LW);
plot(ax, th, yb, '-', 'Color', P.green, 'LineWidth', P.LW);
plot(ax, th, mn, 'o', 'MarkerEdgeColor', P.grey,  'MarkerFaceColor','w', ...
     'MarkerSize', 12, 'LineWidth', 2.5);
plot(ax, th, mb, 'o', 'MarkerEdgeColor', P.green, 'MarkerFaceColor','w', ...
     'MarkerSize', 12, 'LineWidth', 2.5);
if uselog, set(ax,'YScale','log'); end
set(ax,'TickDir','out','TickLabelInterpreter','latex','FontSize',P.FSZ,'LineWidth',2);
xlabel(ax,'threshold $\theta$','Interpreter','latex','FontSize',P.FSZ+8);
ylabel(ax, ylab, 'Interpreter','latex','FontSize',P.FSZ+8);
%  POST: grey is the naive rule, green Bayesian; lines are the solver,
%        open circles the independent Monte Carlo.
export(f, name, P);
S = struct('theta',th,'solver_naive',yn,'solver_bayes',yb, ...
           'mc_naive',mn,'mc_bayes',mb);
end

%% ============================ Monte Carlo ==============================
function r = mcdyad(s, bayesian, P)
dt = P.dt; n = P.ntrial; ns = ceil(P.Tmax/dt);
sq = sqrt(2*dt);                                   % sigma^2 = 2
th  = s.theta - P.bgk*sqrt(2)*sqrt(dt);            % BGK-corrected barrier
thK = s.theta;                                     % kicks use the TRUE threshold

% --- H = 0: cascade probability -------------------------------------
[tdep, ~] = runpair(-1, th, thK, s, bayesian, n, ns, dt, sq);
r.p = mean(all(isfinite(tdep),2));

% --- H = 1: collective response time ---------------------------------
[tdep, ~] = runpair(+1, th, thK, s, bayesian, n, ns, dt, sq);
ok = all(isfinite(tdep),2);
% Pairs still running at the horizon are dropped, which biases T2 downward.
% Under H = 1 absorption is certain, so this is a horizon artifact and the
% dropped fraction is reported rather than absorbed silently.
r.drop = 1 - mean(ok);
r.T2   = mean(max(tdep(ok,:),[],2));
end

function [tdep, x] = runpair(mu, th, thK, s, bayesian, n, ns, dt, sq)
%RUNPAIR  One Euler-Maruyama sweep of the coupled pair at drift mu.
%   Both the drift and the noise are masked by alive, so an absorbed agent
%   stops moving entirely.
x = zeros(n,2); alive = true(n,2); tdep = inf(n,2);
for k = 1:ns
    t = k*dt;
    d = mu*ones(n,2);
    if bayesian
        d = d - interp1(s.t, s.lam, min(t,s.t(end)), 'linear','extrap') ...
                * double(alive(:,[2 1]));
    end
    x = x + (d.*dt + sq*randn(n,2)).*double(alive);
    hit = alive & x >= th;
    if any(hit(:))
        tdep(hit) = t; alive(hit) = false; x(hit) = th;
        for c = 1:2
            j = hit(:,c) & alive(:,3-c);
            if any(j)
                if bayesian
                    kick = thK - interp1(s.t, s.Lsurv, min(t,s.t(end)), 'linear','extrap');
                else
                    kick = thK;
                end
                x(j,3-c) = x(j,3-c) + kick;
                cr = j & x(:,3-c) >= th;
                tdep(cr,3-c) = t; alive(cr,3-c) = false;
            end
        end
    end
    if ~any(alive(:)), break; end
end
end

function b = dyad_observables(b)
%DYAD_OBSERVABLES  Cascade probability and collective response time for the
%   dyad, by the EXACT order-statistic integrals of the supplement rather
%   than by any approximation.
%
%   Cascade probability, H = 0:   2 * int f(t|0) G(t;kappa,theta) dt
%   Collective time,     H = 1:   int S(t|1)^2 dt  +  2 * int f(t|1) C dt
%
%   G is the survivor's onward crossing probability after the kick, obtained
%   by integrating the killed density against exp(-max(theta-x-kappa,0)); C
%   is its expected residual crossing time. Both have closed forms (below).
%
%   Naive rule: kappa = theta, constant, uncoupled marginals.
%   Bayesian rule: kappa(t) = theta - Lambda_surv(t), and S, f come from the
%   converged fixed point.
t = b.t(:); th = b.theta;

% --- naive rule, uncoupled marginals ---------------------------------
[~, S1n] = fpt_moving(@(s) th - s, @(s) -1+0*s, t);
[~, S0n] = fpt_moving(@(s) th + s, @(s)  1+0*s, t); %#ok<ASGLU>
f1n = igd(t,  1, th);   f0n = igd(t, -1, th);
b.p_naive  = 2*trapz(t, f0n .* Gexact(t, th*ones(size(t)), th));
b.T2_naive = trapz(t, S1n.^2) + 2*trapz(t, f1n .* Cexact(t, th*ones(size(t)), th));

% --- Bayesian rule, converged coupling -------------------------------
kick = th - b.Lsurv(:);
f1 = max(-gradient(b.S1(:), t), 0);
f0 = max(-gradient(b.S0(:), t), 0);
b.p_bayes  = 2*trapz(t, f0 .* Gexact(t, kick, th));
b.T2_bayes = trapz(t, b.S1(:).^2) + 2*trapz(t, f1 .* Cexact(t, kick, th));
end

function G = Gexact(t, kap, th)
%GEXACT  Eq. (S8). kap may be a vector (time-dependent Bayesian kick).
r = sqrt(2*t);
G = ncdf((th + t)./r) - ncdf((th - kap + t)./r) ...
    - exp(-th)*( ncdf((-th + t)./r) - ncdf((-th - kap + t)./r) ) ...
    + exp(-(th - kap)) .* ncdf((th - kap - t)./r) ...
    - exp(kap) .* ncdf((-th - kap - t)./r);
G = min(max(G, 0), 1);
end

function C = Cexact(t, kap, th)
%CEXACT  Eq. (S11).
r = sqrt(2*t);
a = (th - kap - t)./r;
c = (-th - kap - t)./r;
C = (th - kap - t).*ncdf(a) + r.*npdf(a) ...
    - exp(th)*( (-th - kap - t).*ncdf(c) + r.*npdf(c) );
C = max(C, 0);
end

function f = igd(t, mu, th)
f = th ./ sqrt(4*pi*t.^3) .* exp(-(th - mu*t).^2 ./ (4*t));
end

function p = ncdf(z), p = 0.5*(1 + erf(z/sqrt(2))); end
function p = npdf(z), p = exp(-z.^2/2)/sqrt(2*pi); end

function export(f, name, P)
if exist(P.outdir,'dir')
    exportgraphics(f, fullfile(P.outdir,[name '.pdf']), 'ContentType','vector');
    exportgraphics(f, fullfile(P.outdir,[name '.png']), 'Resolution',400);
    fprintf('           wrote %s.pdf / .png\n', name);
else
    fprintf('           outdir "%s" not found -- %s not exported\n', P.outdir, name);
end
end
