function out = ddm_ig(t, mu, theta)
% DDM_IG  Single-agent first-passage quantities for the symmetric DDM.
%
%   dxi = mu dt + sigma dW,  sigma^2 = 2  (log-likelihood-ratio convention),
%   absorbing threshold theta > 0, xi(0) = 0.
%
%   out.f  inverse-Gaussian first-passage density   (Eq. 13)
%   out.S  survival function                        (Eq. 14)
%   out.h  hazard f/S
%
%   Under H=1 (mu=+1) absorption is certain and E[T] = theta/mu.
%   Under H=0 (mu=-1) absorption occurs with probability q = exp(-theta).
%
%   ZPK 2026

SIG2 = 2;
t = t(:);
f = zeros(size(t)); S = ones(size(t));
ok = t > 0; tt = t(ok);

f(ok) = theta ./ sqrt(2*pi*SIG2*tt.^3) .* exp(-(theta - mu*tt).^2 ./ (2*SIG2*tt));

z1 = (theta - mu*tt) ./ sqrt(SIG2*tt);
z2 = (-theta - mu*tt) ./ sqrt(SIG2*tt);
% exp(2 mu theta/sigma^2)*Phi(z2) evaluated in logs: stable for large theta
S(ok) = normcdf(z1) - exp(2*mu*theta/SIG2 + log(max(normcdf(z2), realmin)));
S = min(max(S, 0), 1);

out.f = f;
out.S = S;
out.h = f ./ max(S, 1e-300);
out.q = exp(-theta);          % false-alarm probability at mu = -1
end
