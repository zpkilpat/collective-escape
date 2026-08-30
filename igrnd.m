function x = igrnd(mu, lam, n)
%IGRND  Inverse-Gaussian random draws by the Michael-Schucany-Haas method.
%
%   x ~ IG(mean mu, shape lam), n draws. No Statistics Toolbox.
%   Used by calib_pool and the supplementary robustness figures.
%
%   ZPK 2026
y = randn(n,1).^2;
x = mu + mu^2*y/(2*lam) - mu/(2*lam)*sqrt(4*mu*lam*y + mu^2*y.^2);
u = rand(n,1);
sw = u > mu./(mu + x);
x(sw) = mu^2 ./ x(sw);
end

