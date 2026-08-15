function L = linf_ceiling(N)
% LINF_CEILING  Closed-form late-time saturation of the total survival drift.
%
%   L_inf(N) = [ (k+2) - 2 sqrt(k+1) ] / k,   k = N-1        (Eq. 2)
%
%   The physical root of k(1-L)^2 = 4L.  L_inf(2) = 3 - 2*sqrt(2) ~ 0.1716;
%   L_inf -> 1 as N -> inf, and slowly: 1 - L_inf ~ 2/sqrt(N-1).
%
%   ZPK 2026
k = N - 1;
L = ((k + 2) - 2*sqrt(k + 1)) ./ max(k, eps);
L(N <= 1) = 0;
end