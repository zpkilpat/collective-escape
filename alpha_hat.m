function a = alpha_hat(M, qN, TPN)
%ALPHA_HAT  Eq. (13): the discounting rate implied by two group response
%           rates pooled over M responders.
%
%   q1 = 1 - (1-qN)^(1/M)      per-unit false alarm
%   m1 = (1-TPN)^(1/M)         per-unit miss
%   r  = log(q1)/log(m1) = (1+alpha)/(1-alpha)
%   alpha = (r-1)/(r+1)
%
%   Threshold and noise scale cancel in the ratio, so no SNR assumption is
%   needed. Reproduces 0.8836 at M = 7 and 0.9511 at M = 14.
%
%   NOTE. Bootstrap replicates must call this with the RESAMPLED rates.
%   calib_pool('alpha', M) uses the cached point estimates and would collapse
%   the interval width to zero.
%
%   ZPK 2026
q = 1 - (1-qN).^(1./M);
m = (1-TPN).^(1./M);
r = log(q)./log(m);
a = (r - 1)./(r + 1);
end
