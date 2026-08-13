function sol = solve_bayes(N, theta, opts)
%SOLVE_BAYES  Bayesian survival-correction fixed point for an N-agent group.
%
%  Each agent's effective drift while its N-1 neighbors are silent is
%  mu_H - (N-1)*lambda(t); the hazards h_H are the first-passage hazards at
%  that drift; and the correction closes lambda = h_1 - h_0. Solved by damped
%  fixed-point iteration over the moving-boundary primitive in fpt_moving.
%
%  Returns everything make_figS2 and make_figS3 need, so both run without an
%  external solver:
%     sol.N, sol.theta, sol.t
%     sol.lam       converged per-neighbor correction
%     sol.Ltot      (N-1)*lam
%     sol.Linf      closed-form ceiling for comparison
%     sol.S1, sol.S0, sol.h1, sol.h0
%     sol.Lsurv     log[S1/S0] <= 0, the survival log-likelihood ratio
%
%  Usage:  sol = solve_bayes(2, 3.602);
%          sol = solve_bayes([2 5 10 20], 3.602);   % struct array
%
%  ZPK 2026

if nargin < 3, opts = struct; end
if ~isfield(opts,'nt'),   opts.nt   = 4000;  end
if ~isfield(opts,'tol'),  opts.tol  = 1e-4;  end
if ~isfield(opts,'maxit'),opts.maxit= 200;   end
if ~isfield(opts,'quiet'),opts.quiet= false; end

if numel(N) > 1
    sol = solve_bayes(N(1), theta, opts);
    for i = 2:numel(N)
        o = opts; o.warm = sol(i-1);
        sol(i) = solve_bayes(N(i), theta, o);
    end
    return
end

Tmax = max(8*theta, 3*theta^2/max(log(max(N,2)),0.5));
t    = linspace(Tmax/opts.nt, Tmax, opts.nt)';

if isfield(opts,'warm') && ~isempty(opts.warm)
    lam = interp1(opts.warm.t, opts.warm.lam, t, 'linear', 'extrap');
    lam = max(lam, 0);
else
    [~,S1,h1] = fpt_moving(@(s) theta - s, @(s) -1+0*s, t);
    [~,S0,h0] = fpt_moving(@(s) theta + s, @(s)  1+0*s, t);
    lam = max(h1 - h0, 0);
end

for it = 1:opts.maxit
    L = (N-1)*lam;
    IL = cumtrapz(t, L);
    [~,S1,h1] = fpt_moving(@(s) theta - s + interp1(t,IL,s,'linear','extrap'), ...
                           @(s) -1 + interp1(t,L,s,'linear','extrap'), t);
    [~,S0,h0] = fpt_moving(@(s) theta + s + interp1(t,IL,s,'linear','extrap'), ...
                           @(s)  1 + interp1(t,L,s,'linear','extrap'), t);
    new = max(h1 - h0, 0);

    % hold lambda past the point where the H=1 survivor pool underflows,
    % rather than propagating quadrature noise into the fixed point
    bad = S1 < 1e-10;
    if any(bad)
        k = find(~bad, 1, 'last');
        if ~isempty(k), new(bad) = new(k); end
    end

    if     it <= 2,        w = min(0.5, 1.5/N);
    else,                  w = 0.5;
    end
    upd = w*new + (1-w)*lam;
    rel = max(abs(upd - lam)) / max(max(abs(upd)), eps);
    lam = upd;
    if rel < opts.tol, break; end
end

if ~opts.quiet
    fprintf('  N = %3d, theta = %.3f: %3d iters, rel %.1e, Ltail %.4f (Linf %.4f)\n', ...
            N, theta, it, rel, tailval(t, (N-1)*lam, theta), linf_ceiling(N));
end

sol = struct('N',N,'theta',theta,'thetaN',theta,'t',t,'lam',lam, ...
             'Ltot',(N-1)*lam,'Linf',linf_ceiling(N), ...
             'S1',S1,'S0',S0,'h1',h1,'h0',h0, ...
             'Lsurv',log(max(S1,eps)./max(S0,eps)), ...
             'iters',it,'resid',rel);
end

function v = tailval(t, L, th)
w = t/th >= 3 & t/th <= 6.8 & isfinite(L);
if nnz(w) < 5, v = L(end); return; end
x = t(w)/th;  cf = [ones(nnz(w),1), x.^(-2/3)] \ L(w);  v = cf(1);
end
