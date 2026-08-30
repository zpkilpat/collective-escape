function [lo, hi] = wilson(k, n, z)
%WILSON  Wilson score interval for a binomial proportion. No toolbox.
%   Default z = 1.959963985 (95%). Used for every binned rate in the paper.
if nargin < 3, z = 1.959963985; end
p = k/n;  d = 1 + z^2/n;  c = p + z^2/(2*n);
h = z*sqrt(p*(1-p)/n + z^2/(4*n^2));
lo = (c - h)/d;  hi = (c + h)/d;
end

