function [t, lambda, info] = solve_lambda_fixedpoint(N, theta, dt, Tmax, opts)
% SOLVE_LAMBDA_FIXEDPOINT  Bayesian survival correction lambda^(N)(t).
%
%   Self-consistent closure: each agent's effective drift while its N-1
%   neighbours stay silent is mu_H - (N-1)*lambda(t); the hazards h_H are the
%   first-passage hazards at that drift; and the correction must satisfy
%
%       lambda(t) = max( h1[lambda](t) - h0[lambda](t), 0 )
%
%   solved by damped fixed-point iteration.  Rectification enforces that a
%   neighbour's silence is never positive evidence for a threat.
%
%   HORIZON WARNING.  The late-time tail sets L_inf(N) and is the quantity
%   identified with alpha.  It converges slowly, so Tmax must resolve the H=0
%   tail: the default max(8*theta, 3*theta^2/log N) follows the detection
%   scale.  Verify the tail is flat before trusting it -- info.tail_slope
%   reports the drift over the last decade of the grid.
%
%   ZPK 2026

if nargin < 3 || isempty(dt),   dt   = 0.005; end
if nargin < 4 || isempty(Tmax), Tmax = max(8*theta, 3*theta^2/log(max(N,2))); end
if nargin < 5, opts = struct(); end
if ~isfield(opts,'iters'), opts.iters = 200; end
if ~isfield(opts,'tol'),   opts.tol   = 1e-4; end
if ~isfield(opts,'lam0'),  opts.lam0  = []; end

t  = (0:dt:Tmax)';  nt = numel(t);  k = N - 1;

% seed: uncoupled single-agent hazard gap (exact, no solver)
if isempty(opts.lam0)
    h1 = ddm_ig(t, 1, theta).h;
    h0 = ddm_ig(t,-1, theta).h;
    lambda = max(h1 - h0, 0);
else
    lambda = interp1(linspace(0,1,numel(opts.lam0))', opts.lam0(:), ...
                     linspace(0,1,nt)', 'linear', 'extrap');
end

% adaptive relaxation: gentler for large N where the coupling (N-1)lambda is stiff
w = min(0.5, 1.5/N);

% Survival floor. The hazard is f/S; once the H=1 survivor pool underflows,
% BOTH are at the quadrature noise level and their ratio is meaningless. The
% late-time tail is exactly where this bites, so lambda is held at its last
% resolved value beyond that point rather than propagating noise into the
% fixed point. (Without this, N=2 converges to lambda ~ 0.92, which fails its
% own self-consistency condition lambda = (1-lambda)^2/4 by a factor of 500.)
SMIN = 1e-8;

for it = 1:opts.iters
    [S1,~,h1] = survival_drift( 1, k*lambda, theta, dt, nt);
    [~, ~,h0] = survival_drift(-1, k*lambda, theta, dt, nt);
    lam_raw   = max(h1 - h0, 0);

    ok    = S1 > SMIN;
    ilast = find(ok, 1, 'last');
    if isempty(ilast), ilast = 1; end
    lam_raw(ilast+1:end) = lam_raw(ilast);

    lam_new = w*lam_raw + (1-w)*lambda;
    err     = max(abs(lam_new - lambda)) / (max(abs(lam_new)) + 1e-12);
    lambda  = lam_new;
    if err < opts.tol, break; end
end

S1  = survival_drift(1, k*lambda, theta, dt, nt);
ok  = S1 > SMIN;

% Read the tail on a CLEAN WINDOW: past the early-time hump (t > 3*theta) and
% inside the resolved region. lambda(end) is horizon-noisy and must not be used.
lo  = 3*theta;
hi  = 0.85*t(end);
if ~isempty(find(ok,1,'last')), hi = min(hi, t(find(ok,1,'last'))); end
if lo >= hi, lo = 0.5*hi; end
win = (t >= lo) & (t <= hi);
if nnz(win) < 5, win = ok & (t > 0.5*hi); end

tw = t(win);  lw = lambda(win);
info.iters      = it;
info.err        = err;
info.tail       = k*mean(lw);
info.tail_slope = k*(lw(end) - lw(1)) / max(tw(end) - tw(1), eps);
info.window     = [lo hi];
info.frac_resol = nnz(ok)/nt;
info.Linf       = linf_ceiling(N);
info.theta      = theta;
info.Tmax       = Tmax;

fprintf(['  N=%3d theta=%.3f: %3d iters, err=%.1e, (N-1)lam plateau=%.4f ', ...
         '(L_inf=%.4f), slope=%+.2e, window t/th=[%.1f,%.1f], resolved=%.0f%%\n'], ...
        N, theta, it, err, info.tail, info.Linf, info.tail_slope, ...
        info.window(1)/theta, info.window(2)/theta, 100*info.frac_resol);
end
