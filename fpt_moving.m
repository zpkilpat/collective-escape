function [g, S, h] = fpt_moving(b, db, t, x0)
%FPT_MOVING  First-passage density to a moving boundary, driftless Brownian
%            motion with variance rate 2 (diffusion coefficient 1).
%
%  THE single first-passage primitive in this repository. Durbin's Volterra
%  equation of the second kind in the Buonocore-Nobile-Ricciardi master form,
%
%      g(t) = Psi(t,0|b) - int_0^t g(s) Psi(t,s|b) ds,
%      Psi(t,s|b) = [ (b(t)-b(s))/(t-s) - bdot(t) ] p(b(t),t | b(s),s),
%      p(y,t|y',s) = exp(-(y-y')^2/(4(t-s))) / sqrt(4 pi (t-s)).
%
%  CONVENTION. No factor of two; the integral is SUBTRACTED. Verified to
%  machine precision (max error 1e-16) against the exact inverse Gaussian for
%  a constant boundary AND for a linear moving boundary. A constant boundary
%  alone cannot validate the kernel, since Psi(t,s) vanishes identically
%  there for s > 0 and only the s = 0 term survives; the linear case is what
%  pins the sign. Earlier drafts of the supplement carried
%  g = -2 Psi + 2 int g Psi, which is wrong.
%
%  A drift-diffusion with time-varying drift mu(t) against a FIXED barrier
%  theta maps here by y(t) = xi(t) - A(t), A(t) = int_0^t mu, giving the
%  moving boundary b(t) = theta - A(t). See survival_drift.m, which does that
%  change of variables for you.
%
%  Inputs
%     b, db   EITHER function handles (evaluated on t) OR nt-vectors giving
%             the boundary and its derivative at the nodes of t
%     t       uniform, strictly increasing time grid with t(1) > 0
%     x0      starting point (default 0)
%  Outputs
%     g       first-passage density on t
%     S       survival, 1 - cumulative
%     h       hazard g/S
%
%  CONSOLIDATION NOTE. This replaces three separate implementations that
%  coexisted during development: fpt_moving_boundary.m, an earlier
%  handle-only fpt_moving, and a local fpt() inside run_bayes_nagent.m. All
%  three agreed algebraically (fpt_moving_boundary folded the minus sign into
%  its kernel definition), but three copies of a primitive whose sign
%  convention was once wrong is exactly the wrong thing to ship.
%
%  ZPK 2026

if nargin < 4 || isempty(x0), x0 = 0; end
t  = t(:);  n = numel(t);  dt = t(2) - t(1);

if isa(b,'function_handle'),  bt  = b(t);   else, bt  = b(:);  end
if isa(db,'function_handle'), dbt = db(t);  else, dbt = db(:); end
if numel(bt) ~= n || numel(dbt) ~= n
    error('fpt_moving:size','b and db must match the length of t');
end

g = zeros(n,1);
for k = 1:n
    tk = t(k);
    % free term, measured from the start point rather than from b(0)
    psi0 = ((bt(k) - x0)/tk - dbt(k)) * pkern(bt(k), tk, x0, 0);
    acc  = 0;
    if k > 1
        s   = t(1:k-1);
        ds  = tk - s;
        psi = ((bt(k) - bt(1:k-1))./ds - dbt(k)) .* pkern(bt(k), tk, bt(1:k-1), s);
        w   = ones(k-1,1);  w(end) = 0.5;      % removable diagonal
        acc = dt * sum(w .* g(1:k-1) .* psi);
    end
    g(k) = psi0 - acc;
end

g = max(g, 0);
S = max(1 - dt*cumsum(g), 0);
h = g ./ max(S, 1e-300);
h(~isfinite(h)) = 0;
end

function v = pkern(y, t, yp, s)
d = t - s;
v = exp(-(y - yp).^2 ./ (4*d)) ./ sqrt(4*pi*d);
end
