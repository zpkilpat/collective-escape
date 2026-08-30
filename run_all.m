%% RUN_ALL  Regenerate every figure panel, in dependency order.
%
%  Writes standalone PDF and PNG panels into output/. Multi-panel figures are
%  assembled externally; panel letters are not drawn here.
%
%  FIGURE SET. Four main figures and seven supplementary, matching the PNAS
%  submission. See PANELMAP below for which script writes which panel; several
%  scripts do not emit the panel letter their own name suggests.
%
%  DEPENDENCY ORDER MATTERS IN THREE PLACES.
%    1. solve_bayes must precede make_fig3 and make_figS2, and must run at
%       theta_1 = 3.5678. Panel 3A reads the converged correction from that
%       struct rather than rebuilding it.
%    2. calib_pool caches its calibration in a persistent. Changing the
%       response rates without calling calib_pool('reset') silently inverts
%       against the old curve.
%    3. make_fig2 writes figS2C as well as its own four panels, so it must run
%       before Fig. S2 is assembled. It does not need to precede make_figS2,
%       which writes figS2A and figS2B only.
%
%  The N = 100 solve is the slow step, twenty-five fixed-point iterations of
%  two O(nt^2) Volterra sweeps each. Set SKIP_SOLVER = true to load the saved
%  struct instead; the solve below writes it.
%
%  SCRIPTS ARE CALLED THROUGH run_script. make_fig1, make_fig2 and make_figS1
%  are scripts, not functions, and open with "clear". Called directly they
%  would wipe this workspace, taking CSV, OUTDIR, MATF, SKIP_SOLVER and t0
%  with them. Running them inside a local function confines the clear.
%
%  Data S1 and S2 are not redistributed; see README.md.
%
%  ZPK 2026
clear; close all;
CSV    = 'sciadv_adt8600_data_s1.csv';
XLSX   = 'sciadv_adt8600_data_s2.xlsx';
OUTDIR = 'output';
MATF   = 'bayes_theta3568.mat';
SKIP_SOLVER = false;          % true: load MATF instead of re-solving
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end
calib_pool('reset');          % start from a clean calibration
t0 = tic;

%% ---- theory panels needing no data -----------------------------------
% make_figS4 sits here rather than with the empirical figures: it is a
% simulation study of the calibration and reads no data file, so gating it on
% the CSV would skip it for no reason.
fprintf('\n===== Fig 1 (schematic B, C, survival LLR D) =====\n'); run_script('make_fig1');
fprintf('\n===== Fig 2 (dyad) + Fig S2C =====\n');                 run_script('make_fig2');
fprintf('\n===== Fig S1 (single-agent frontier) =====\n');         run_script('make_figS1');
fprintf('\n===== Fig S3 (ceiling and hump) =====\n');              run_script('make_figS3');
fprintf('\n===== Fig S4 (pooling-count robustness) =====\n');      out.figS4 = make_figS4(OUTDIR);

%% ---- the N-agent solve, then everything that reads it ----------------
% NOTE the anchor is defined twice: here from the observed statistic, and
% inside make_figS3 as a literal 3.5678. They agree, but if the rates ever
% move only this one follows.
th1 = calib_pool('theta1', calib_pool('invert', 3.93));
fprintf('\n===== N-agent fixed point at theta_1 = %.4f =====\n', th1);
if SKIP_SOLVER && exist(MATF,'file')
    fprintf('  loading %s\n', MATF);
    S = load(MATF);  bayes = S.bayes;
else
    bayes = solve_bayes([2 5 10 20 40 100], th1);
    save(MATF, 'bayes', 'th1');       % so SKIP_SOLVER has something to find
    fprintf('  wrote %s\n', MATF);
end

fprintf('\n===== Fig 3 =====\n');                        out.fig3  = make_fig3(bayes, OUTDIR);
fprintf('\n===== Fig S2A, S2B (heuristic vs fixed point) =====\n');
out.figS2 = make_figS2(bayes, OUTDIR);

%% ---- empirical panels -------------------------------------------------
if ~exist(CSV,'file')
    warning('run_all:nodata', ...
        ['%s not found -- skipping the empirical figures. Download Data S1 ' ...
         'and S2 from the Pacher et al. archive; see README.md.'], CSV);
else
    fprintf('\n===== Fig 4 =====\n');    out.fig4  = make_fig4(CSV, OUTDIR);
    fprintf('\n===== Fig S5 =====\n');   out.figS5 = make_figS5(CSV, OUTDIR);
    fprintf('\n===== Fig S6 =====\n');   out.figS6 = make_figS6(CSV, OUTDIR);
    if exist(XLSX,'file')
        fprintf('\n===== Fig S7 =====\n');
        out.figS7 = make_figS7(CSV, XLSX, OUTDIR);
    else
        warning('run_all:nos2', '%s not found -- Fig S7 panels B, C skipped.', XLSX);
        out.figS7 = make_figS7(CSV, '', OUTDIR);
    end
end

fprintf('\n===== done in %.1f min; panels in %s/ =====\n', toc(t0)/60, OUTDIR);

%% ---- local -------------------------------------------------------------
function run_script(name)
%RUN_SCRIPT  Execute a figure SCRIPT in an isolated workspace, so that a
%   "clear" inside it cannot reach the driver's variables. Local functions in
%   scripts require R2016b or later.
run(name);
end
