# Rerun order

`make_figS7.m` deleted (byte-identical to `make_figS6.m` apart from the
export names, a leftover from the renumbering).

## Gate

`read_pacher.m` was not in the upload, so it is the one file not verified
against the corrected cluster key. All three cluster-dependent scripts now
call `assert_key(D)` and warn loudly if the key does not give 73/47/26.
Check it once before anything else:

    D = read_pacher('sciadv_adt8600_data_s1.csv');
    numel(unique(D.clust))    % must be 73

If it is not 73, fix the loader first. Nothing downstream is worth running.

## Order

    % 1  anchor chain -- theta_1 now 3.5678 in make_fig1/3/S1
    d2 = solve_bayes(2, 3.5678);
    make_fig3;                              % panels A, D, E
    make_figS1;

    % 2  N sweep, feeds two figures
    bayes = solve_bayes([2 5 10 20 40], 3.5678);
    make_figS2(bayes, 'output');
    make_fig4(bayes, 'output');             % q1_CI now [8.825, 29.403]

    % 3  self-contained
    make_figS3('output');                   % panel B now at theta_1(M)

    % 4  cluster-dependent, after the gate
    make_figS4('sciadv_adt8600_data_s1.csv', 'output');
    make_figS5('sciadv_adt8600_data_s1.csv', 'output');
    make_figS6('sciadv_adt8600_data_s1.csv', 'sciadv_adt8600_data_s2.xlsx', 'output');

    % 5  last -- rerunning fixes the Mhat = 14 inset label on its own
    make_fig5('sciadv_adt8600_data_s1.csv');   % 2nd arg is an opts STRUCT

## Expected console lines

    make_figS5  panel C: 73 clusters, 47 carrying attacks, 26 carrying flybys
    make_figS6  panel A: n = 125 timed attacks over 42 clusters
                         log-log slope -0.429, clustered 95% CI [-0.841, 0.071]
    make_fig5   bootstrap: M in [8.82, 29.40], alpha in [0.906, 0.982]
                           4000/4000 valid, P = 1.0000, min excess = 0.2043

## Numbers these runs are meant to settle

  - make_figS3 panel A: relative spread of lambda/Tbar over theta at
    M = 5,10,14,20,40. The Fig. S3 caption says 54/18/12/8/5; the console
    has said 52/19/11/6/1.5. The M = 40 entry differs by 3x.
  - make_figS3 panel B: censoring sensitivity, now at the self-consistent
    threshold. The caption's "5% reads as 25, 20% as 40" was produced at a
    hardcoded theta = 2.0 and will move.
  - make_figS3 panel C: the 0.38 -> 0.07 correlation gap.
  - make_figS4 panel A: the wedge minimum, quoted as 0.26.
  - make_figS4 panel B: the bootstrap median M, quoted as 14.1.
  - make_figS4 panel C: heterogeneity triple. Caption says 0.451/0.406/0.342
    at CV = 0/0.6/1.0; the console has said 0.441/0.414/0.369.
  - make_figS2 panel B: detection-time error series, quoted as 1.2% to 6.9%.
  - solve_bayes(2,...): the dyad numbers still written at theta_1 = 3.317,
    i.e. "cuts q_dyad by roughly a third", the 2.7% / 8.4% pair, and the
    survival discrepancies 3.1e-2 / 2.2e-3.

## Not touched

`make_fig2.m` is expository at theta = 2.0 and needs no rerun.

## Still open in the repo

  - `make_fig5.m` carries its own `load_pacher`, duplicating `read_pacher.m`.
  - `make_fig5(csv, opts)` takes an options struct where `make_figS3`--`S6`
    take an output-directory string in the same position. Passing 'output'
    to `make_fig5` errors in `defaults()`.
