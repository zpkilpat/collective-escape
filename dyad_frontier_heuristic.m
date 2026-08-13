function [q_casc, T2, ET1, gap] = dyad_frontier_heuristic(theta, alpha, dt, Tfac, Nx)
% DYAD_FRONTIER_HEURISTIC  Closed-form dyad observables for the constant-alpha
% heuristic, via the survivor's killed position density.
%
%   This is the computation the paper describes, not a per-agent suppression
%   factor applied to the naive rate. Two things encode the heuristic:
%
%     (1) ANTI-DRIFT DURING SILENCE. While the neighbour is still silent both
%         agents drift at mu_H - alpha, so the first-departure density AND the
%         survivor's position law are evaluated at the discounted drift
%         mu_tilde_1 = 1 - alpha (H=1), mu_tilde_0 = -(1 + alpha) (H=0).
%     (2) GROWING KICK. At the first departure (time s) the survivor is kicked
%         by J(s) = theta + alpha*s, the constant-rate form of the exact
%         Bayesian update theta - Lambda^surv(s) under Lambda^surv ~ -alpha*s.
%         POST-KICK THE SURVIVOR IS BARE: the anti-drift switches off, so
%         completion runs at the plain drift mu_H.
%
%   The survivor's pre-kick position law is the killed density on x < theta,
%   by the method of images with the image weight exp(mu*theta) that makes the
%   combination vanish on the boundary (Eqs. S3-S4). It is UN-NORMALISED: its
%   integral is the survival probability, so the survival weight is carried
%   inside and must not be applied again.
%
%   At alpha = 0 this reduces exactly to the naive pulsatile rule with the
%   fixed kick kappa = theta -- NOT to independent agents. Independence is
%   kappa = 0, which gives q^2. run_selftest cross-checks the alpha = 0 case
%   against dyad_cascade_analytic (the closed-form G-integral of S2); the two
%   derivations share no code and agree to six significant figures.
%
%   Inputs
%     theta : threshold
%     alpha : constant discount rate (alpha = 0 -> naive with kappa = theta)
%     dt    : time quadrature step        (default 0.0015)
%     Tfac  : integrate t over [dt, Tfac*theta]  (default 20)
%     Nx    : position grid points        (default 5000)
%
%   Outputs
%     q_casc : P(both depart | H=0)
%     T2     : E[T_(2) | H=1], collective response time
%     ET1    : E[T_(1) | H=1], first-departure time
%     gap    : E[T_(2) - T_(1) | H=1]
%
%   ZPK 2026

if nargin < 3 || isempty(dt),   dt   = 0.0015; end
if nargin < 4 || isempty(Tfac), Tfac = 20;     end
if nargin < 5 || isempty(Nx),   Nx   = 5000;   end

sv = (dt : dt : Tfac*theta)';
xg = linspace(-25, theta, Nx)';

%% ===================== H = 0: false-alarm cascade =====================
mu  = -1;
mud = mu - alpha;                       % (1) discounted silence-phase drift
f   = igd(sv, mud, theta);              % trigger density at the discounted drift

q_casc = 0;
for i = 1:numel(sv)
    s    = sv(i);
    pa   = p_abs(xg, s, mud, theta);    % survivor position law, killed, un-normalised
    x0   = xg + theta + alpha*s;        % (2) growing kick J(s)
    comp = eventual_hit(x0, mu, theta); % post-kick survivor is BARE (drift mu)
    comp(x0 >= theta) = 1;              % already across at the kick
    q_casc = q_casc + 2*f(i) * trapz(xg, pa .* comp) * dt;
end

%% ===================== H = 1: collective response time =====================
mu  = 1;
mud = mu - alpha;
ET1 = trapz(sv, survcf(sv, mud, theta).^2);      % dyad first departure
f   = igd(sv, mud, theta);

gap_num = 0;  gap_den = 0;
for i = 1:numel(sv)
    s   = sv(i);
    pa  = p_abs(xg, s, mud, theta);
    x0  = xg + theta + alpha*s;
    rem = max(theta - x0, 0) / mu;      % bare mean remaining first-passage time
    gap_num = gap_num + 2*f(i) * trapz(xg, pa .* rem) * dt;
    gap_den = gap_den + 2*f(i) * trapz(xg, pa)        * dt;
end
gap = gap_num / max(gap_den, 1e-300);
T2  = ET1 + gap;
end

%% ---------------- local analytic primitives (sigma^2 = 2) ----------------

function f = igd(t, mu, th)
% Inverse-Gaussian first-passage density.
t = t(:);  f = zeros(size(t));  ok = t > 0;  tt = t(ok);
f(ok) = th ./ sqrt(4*pi*tt.^3) .* exp(-(th - mu*tt).^2 ./ (4*tt));
end

function S = survcf(t, mu, th)
% Survival function.
t = t(:);  S = ones(size(t));  ok = t > 0;  tt = t(ok);
S(ok) = normcdf((th - mu*tt)./sqrt(2*tt)) ...
      - exp(mu*th).*normcdf((-th - mu*tt)./sqrt(2*tt));
S = min(max(S,0),1);
end

function p = p_abs(x, s, mu, th)
% Killed position density on x < theta by the method of images. The image
% source sits at the reflection 2*theta of the origin across the barrier,
% weighted by exp(mu*theta) so the combination vanishes at x = theta.
% UN-NORMALISED: int p dx = S(s), the survival probability.
x    = x(:);
base = exp(-(x - mu*s).^2 ./ (4*s)) ./ sqrt(4*pi*s);
img  = exp(mu*th) .* exp(-((x - 2*th) - mu*s).^2 ./ (4*s)) ./ sqrt(4*pi*s);
p    = max(base - img, 0);
end

function p = eventual_hit(x0, mu, th)
% Probability of ever reaching th from x0 at constant drift mu, sigma^2 = 2.
% For mu < 0 this is exp(mu*(th-x0)); at x0 = 0, mu = -1 it recovers exp(-th).
x0 = x0(:);
if mu >= 0
    p = ones(size(x0));
else
    p = min(max(exp(mu*(th - x0)), 0), 1);
end
end
