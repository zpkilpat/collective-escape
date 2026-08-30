function [Mhat, out] = fit_pool_mle(t, opts)
%FIT_POOL_MLE  Pooling count by maximum likelihood, with a FREE TIME SCALE.
%              A CROSS-CHECK on the moment-ratio estimator, not the headline.
%
%  READ THIS BEFORE USING. The model carries no timescale of its own: units
%  are fixed by sigma^2 = 2 and mu = 1, so model time is not seconds. At
%  M = 13.5 the predicted mean latency is about 1.8 model units against an
%  observed 4.92 s. A likelihood written on raw seconds absorbs that
%  mismatch by dragging M down -- it returned 4.3 on the molly data, against
%  13.5 from the ratio. The fix is a free scale, t_obs = s * t_model,
%  profiled out at each M, which makes the estimator scale-invariant
%  (multiplying every latency by 3 leaves Mhat unchanged).
%
%  WITH THE SCALE FREE THIS ESTIMATOR HAS NO ADVANTAGE. Recovery at n = 125
%  is [10.1, 22.6] against [10.3, 22.5] for the moment ratio -- identical. A
%  scale-free model can only be identified by scale-free features of the
%  data, so the dimensionless ratio already extracts what there is.
%
%  M is not a free parameter sitting alongside the threshold. Once M is
%  chosen, the two response rates fix both the discounting rate and the
%  threshold,
%       alpha(M) = Eq. (13) at the M-pooled rates,
%       theta(M) = -log(1 - (1-q)^(1/M)),
%  so the model is a ONE-PARAMETER family and the likelihood can be profiled
%  directly over M. The observed latency is the minimum of M first passages,
%  with density
%       f_(1)(t) = M f(t) S(t)^(M-1),
%  at drift 1-alpha(M) against theta(M).
%
%  Inputs
%     t       vector of latencies, seconds. If omitted, Data S1 is read from
%             the working directory via read_pacher and the 125 valid
%             latencies are used.
%     opts    .Mmin, .Mmax   search bounds (default 2, 60)
%             .ngrid         profile grid (default 600)
%             .qN, .TPN      response rates (default the molly values)
%             .censor        if given, the latencies are treated as
%                            right-censored at this time and the likelihood
%                            uses S(c)^M for the unobserved fraction
%             .ncens         number of censored observations (with .censor)
%  Outputs
%     Mhat    the maximizer
%     out     .Mgrid, .loglik, .ci (likelihood-ratio 95%), .alpha, .theta
%
%  ZPK 2026

if nargin < 1 || isempty(t)
    D = read_pacher();                       % finds Data S1 in the working dir
    t = D.tdet;                              % attacks only -- see read_pacher
    fprintf('  loaded %d detection latencies from Data S1 (mean %.2f s)\n', ...
            numel(t), mean(t));
end
if nargin < 2, opts = struct; end
if ~isfield(opts,'Mmin'),  opts.Mmin  = 2;    end
if ~isfield(opts,'Mmax'),  opts.Mmax  = 60;   end
if ~isfield(opts,'ngrid'), opts.ngrid = 600;  end
if ~isfield(opts,'qN'),    opts.qN    = 0.3210; end
if ~isfield(opts,'TPN'),   opts.TPN   = 0.7175; end

t = t(:); t = t(isfinite(t) & t > 0);
Mg = linspace(opts.Mmin, opts.Mmax, opts.ngrid)';
ll = nan(size(Mg));  sh = nan(size(Mg));
for i = 1:numel(Mg)
    % profile the time scale out at each M
    f = @(ls) -loglik(t, Mg(i), exp(ls), opts);
    ls = fminbnd(f, -3, 4);
    sh(i) = exp(ls);  ll(i) = -f(ls);
end

[~, k] = max(ll);
Mhat = Mg(k);

% likelihood-ratio interval, chi2(1) at 95% is 3.841, so half of it
cut = ll(k) - 1.920;
in  = ll >= cut;
lo  = Mg(find(in, 1, 'first'));
hi  = Mg(find(in, 1, 'last'));

out = struct('Mgrid', Mg, 'loglik', ll, 'scale', sh, 'ci', [lo hi], ...
             'Mhat', Mhat, 'shat', sh(k), ...
             'alpha', alpha_of(Mhat, opts), 'theta', theta_of(Mhat, opts), ...
             'n', numel(t));

if nargout == 0 || nargin < 1
    fprintf(['  Mhat = %.2f  [%.2f, %.2f]   alpha = %.4f   theta = %.4f' ...
             '   scale = %.2f s  (n = %d)\n'], ...
            Mhat, lo, hi, out.alpha, out.theta, out.shat, numel(t));
    fprintf('  cross-check only -- headline estimator is calib_pool\n');
end
end

%% ---------------------------------------------------------------------
function L = loglik(t, M, s, opts)
%LOGLIK  t_obs = s * t_model, so the density picks up a 1/s Jacobian.
tm = t / s;
th = theta_of(M, opts);
nu = 1 - alpha_of(M, opts);
lf = log(max(igpdf(tm, nu, th), realmin));
lS = log(max(survf(tm, nu, th), realmin));
L  = numel(t)*log(M) - numel(t)*log(s) + sum(lf) + (M-1)*sum(lS);
if isfield(opts,'censor') && ~isempty(opts.censor) && isfield(opts,'ncens')
    L = L + opts.ncens * M * log(max(survf(opts.censor/s, nu, th), realmin));
end
end

function a = alpha_of(M, opts)
q = 1 - (1-opts.qN).^(1./M);  m = (1-opts.TPN).^(1./M);
r = log(q)./log(m);  a = (r-1)./(r+1);
end

function th = theta_of(M, opts)
th = -log(1 - (1-opts.qN).^(1./M));
end

function f = igpdf(t, nu, th)
f = th ./ sqrt(4*pi*t.^3) .* exp(-(th - nu*t).^2 ./ (4*t));
end

function S = survf(t, nu, th)
S = ncdf((th - nu*t)./sqrt(2*t)) - exp(nu*th).*ncdf((-th - nu*t)./sqrt(2*t));
S = min(max(S, 0), 1);
end

function p = ncdf(z), p = 0.5*(1 + erf(z/sqrt(2))); end

