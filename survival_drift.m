function [S, f, h] = survival_drift(mu, lambda_tot, theta, dt, nt)
%SURVIVAL_DRIFT  Survival, density and hazard for a DDM with time-varying
%                social anti-drift against a FIXED threshold.
%
%   dxi = (mu - lambda_tot(t)) dt + sigma dW,  sigma^2 = 2,  absorb at theta.
%
%   The deterministic shift y(t) = xi(t) - int_0^t (mu - lambda_tot) ds turns
%   this into a driftless process against the moving boundary
%       b(t) = theta - int_0^t (mu - lambda_tot) ds,
%   which fpt_moving solves. Both hypotheses integrate FORWARD in the same
%   physical-time direction -- the property the earlier backward-Kolmogorov
%   implementation lacked, and the reason the saturation ceiling is a hazard
%   balance rather than a discretization artifact.
%
%   lambda_tot may be a scalar (constant discount, the heuristic) or an
%   nt-vector (the Bayesian fixed point).
%
%   ZPK 2026

if isscalar(lambda_tot), lambda_tot = lambda_tot*ones(nt,1); end
lambda_tot = lambda_tot(:);

t     = (1:nt)'*dt;                       % fpt_moving needs t(1) > 0
drift = mu - lambda_tot;
A     = cumsum(drift)*dt;                 % int_0^t (mu - lambda) ds
b     = theta - A;
bdot  = -drift;

[f, S, h] = fpt_moving(b, bdot, t);
end
