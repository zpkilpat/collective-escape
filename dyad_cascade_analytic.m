function [q_casc, cascade_time] = dyad_cascade_analytic(theta, kappa, dt, Tmax)
% DYAD_CASCADE_ANALYTIC  Closed-form naive-dyad cascade observables.
%
%   Implements the order-statistic integrals of Supplementary S2 for the naive
%   pulsatile rule (fixed kick kappa, no anti-drift):
%
%     P(cascade | H=0)      = 2 int f(t|H=0) G(t;kappa,theta) dt      (S5)-(S7)
%     E[cascade time | H=1] = 2 int f(t|H=1) C(t;kappa,theta) dt      (S8)-(S9)
%
%   G is the survivor's eventual false-alarm probability after the kick,
%   integrated against the H=0 killed density; C is the expected residual
%   crossing time after the kick under H=1.  sigma^2 = 2 throughout.
%
%   At kappa = 0, q_casc collapses to q^2 = exp(-2 theta) (independent agents).
%
%   ZPK 2026

if nargin < 3 || isempty(dt),   dt   = 0.002; end
if nargin < 4 || isempty(Tmax), Tmax = 25*theta; end

t = (dt:dt:Tmax)';
s2 = sqrt(2*t);

% ---- G(t;kappa,theta), Eq. (S6) ----------------------------------------
G = normcdf((theta + t)./s2) - normcdf((theta - kappa + t)./s2) ...
  - exp(-theta) * ( normcdf((-theta + t)./s2) - normcdf((-theta - kappa + t)./s2) ) ...
  + exp(-(theta - kappa)) * normcdf((theta - kappa - t)./s2) ...
  - exp(kappa) * normcdf((-theta - kappa - t)./s2);

f0 = ddm_ig(t, -1, theta).f;
q_casc = 2 * trapz(t, f0 .* G);

% ---- C(t;kappa,theta), Eq. (S8) ----------------------------------------
a = (theta - kappa - t)./s2;
b = (-theta - kappa - t)./s2;
C = (theta - kappa - t).*normcdf(a) + s2.*normpdf(a) ...
  - exp(theta) * ( (-theta - kappa - t).*normcdf(b) + s2.*normpdf(b) );

f1 = ddm_ig(t, 1, theta).f;
cascade_time = 2 * trapz(t, f1 .* max(C,0));
end

