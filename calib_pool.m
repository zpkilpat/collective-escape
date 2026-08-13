function varargout = calib_pool(varargin)
%CALIB_POOL  Self-consistent calibration of the pooling count M.
%
%  The scale-free latency statistic lambda/Tbar is NOT threshold-free at the
%  M of interest. Measured spread across theta in [1,6] is 52% at M=5, 19% at
%  M=10, 11% at M=14, 6% at M=20, 1.5% at M=40. Calibrating at an arbitrary
%  theta therefore moves the estimate: theta=1 returns Mhat=14.7, theta=3.6
%  returns 13.6. This is the source of the drift across earlier runs, and it
%  is not Monte Carlo noise.
%
%  The threshold is not free. Once M is chosen, the false alarm rate fixes it,
%      theta_1(M) = -log(1 - (1-q)^(1/M)),
%  so the calibration is run at theta_1(M) AT EACH GRID POINT and the loop
%  closes. Nothing is tuned. Grid and draw count are set here and nowhere
%  else, so every figure inverts against the same curve.
%
%  Usage:
%     M            = calib_pool('invert', stat)     % statistic -> M
%     [Mg, sg]     = calib_pool('grid')             % the calibration curve
%     th           = calib_pool('theta1', M)        % threshold at M
%     a            = calib_pool('alpha', M)         % alpha at M (Eq. 13)
%     calib_pool('reset')                           % force recompute
%     calib_pool('rates', q, TP)                    % set rates, resets
%
%  Reproduces: stat 3.93 -> Mhat 13.5, theta_1 3.57, alpha 0.949,
%              L_inf(M) 0.573, gap 0.376. Stable to +-0.07 across seeds.
%
%  ZPK 2026

persistent Mg sg QN TPN
if isempty(QN), QN = 0.3210; TPN = 0.7175; end

GRID  = [2:1:40, 45, 50, 60, 75, 90]';
NDRAW = 2e6;                       % seed spread +-0.04 at this count
SEED  = 20260811;

switch lower(varargin{1})
    case 'reset'
        Mg = []; sg = []; return
    case 'rates'
        QN = varargin{2}; TPN = varargin{3}; Mg = []; sg = []; return
    case 'theta1'
        varargout{1} = theta1(varargin{2}, QN); return
    case 'alpha'
        varargout{1} = alpha_hat(varargin{2}, QN, TPN); return
end

if isempty(Mg)
    st = rng; rng(SEED);
    Mg = GRID; sg = nan(size(Mg));
    for i = 1:numel(Mg)
        M  = Mg(i);
        th = theta1(M, QN);                    % self-consistent, not free
        a  = alpha_hat(M, QN, TPN);
        x  = igrnd(th/(1-a), th^2/2, NDRAW);
        k  = floor(NDRAW/M);
        sg(i) = shape_over_mean(min(reshape(x(1:k*M), k, M), [], 2));
    end
    rng(st);
    if ~all(diff(sg) > 0)
        error('calib_pool:monotone', ...
              'calibration curve is not monotone -- check rates or draw count');
    end
end

switch lower(varargin{1})
    case 'grid'
        varargout{1} = Mg; varargout{2} = sg;
    case 'invert'
        M = interp1(sg, Mg, varargin{2}, 'linear', 'extrap');
        varargout{1} = min(max(M, Mg(1)), Mg(end));
    otherwise
        error('calib_pool:usage','unknown mode "%s"', varargin{1});
end
end

%% ---------------------------------------------------------------------
function th = theta1(M, qN)
th = -log(1 - (1-qN).^(1./M));
end


function s = shape_over_mean(t)
m = mean(t);  s = (1/(mean(1./t) - 1/m)) / m;
end
