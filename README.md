# Nonparametric Partial Directed Coherence (iPDC): tutorial code

MATLAB code for *A Tutorial on Nonparametric Estimation of Partial Directed
Coherence and Network Hierarchy* (Takahashi and Zhang),
<https://github.com/takahashiyd/nppdc>. The pipeline estimates
information PDC (iPDC) without fitting a VAR model:

    multitaper cross-spectrum  ->  coherency-preconditioned Wilson factorization  ->  iPDC
    ->  closed-form significance threshold, confidence band, bias correction
    ->  Hodge (gradient) decomposition of the directed flow with jackknife inference

## Requirements

MATLAB R2021a or later, Signal Processing Toolbox (`dpss`; the two Allen LFP
demos also use `butter`, `filtfilt` and `iirnotch`) and Statistics and Machine
Learning Toolbox (`finv`, `norminv`, `normcdf`, `chi2cdf`). No other
dependencies. Tested on MATLAB R2021a (Windows).

## Layout

    run_all.m   regenerate every figure and table (or a subset) into figures/
    src/        the library (add this folder to your path to use the method)
    demos/      one script per figure or table of the tutorial
    figures/    the PNGs used in the manuscript
    data/       the extracted Allen LFP subset (see data/README.md)

## Quick start

```matlab
addpath src
fs = 200;                                           % sampling rate (Hz)
A  = zeros(3, 3, 2);                                % the tutorial's reference VAR(2)
A(:,:,1) = [0.55 0 0; 0.40 0.556 0; 0 0.50 0.40];
A(:,:,2) = [-0.25 0 0; 0 -0.81 0; 0 0 -0.20];
Sigma = diag([1.0, 0.8, 1.2]);
X  = sim_var(4096, A, Sigma);                       % or your own N x T data
[S, f, info] = multitaper_csd(X, fs, 4, 7, 512, 0.5, 512);   % NW=4, K=7, segLen=T/8, 50% overlap
[H, Sigma_hat] = wilson_factor_normalized(S);
I  = ipdc(H, Sigma_hat);                            % nfft x N x N, I(:,i,j) = j -> i
[thr, ci_lo, ci_hi] = ipdc_asymp(H, Sigma_hat, info.effDOF, 0.05);   % F-form threshold, Fisher-z CI
I_corr = ipdc_bias_correct(I, info.effDOF, size(X,1));
out = ipdc_hodge_jack(X, fs, 4, 7, 256, 0.5, 256, 0.05);        % ranking scores +/- CI
```

The sectioned walk-through is `demos/nppdc_tutorial.m`.

## Library (`src/`)

| Function | Purpose |
|---|---|
| `multitaper_csd` | segment-averaged multitaper CSD; returns the effective DOF `info.effDOF = (2K+1) M_indep - 2` |
| `wilson_factor`, `wilson_factor_normalized` | Wilson minimum-phase factorization, plain and coherency-preconditioned (recommended) |
| `scalar_min_phase_factor` | cepstral scalar minimum-phase factor used by the preconditioner |
| `ipdc` | iPDC from `(H, Sigma)`, full quadratic form in `Sigma^{-1}` |
| `ipdc_asymp` | F-form H0 threshold and Fisher-z confidence band |
| `ipdc_bias_correct` | closed-form bias correction, `E|i|^2 = |i|^2 + (1-|i|^2)^2/n`, `n = L - N + 1` |
| `ipdc_hodge`, `ipdc_hodge_jack` | Hodge gradient/ranking decomposition; delete-one-block jackknife inference |
| `ipdc_from_var`, `ipdc_hawkes_theory` | closed-form iPDC for VAR and exponential-kernel Hawkes ground truth |
| `sim_var`, `sim_hawkes_mv` | simulators (VAR(p); Hawkes by Ogata thinning) |
| `ipdc_validate_mc` | Monte Carlo calibration helper |
| `paper_palette`, `style_paper_figure` | figure styling |

## Reproducing the tutorial

`run_all` runs every demo and writes to `figures/`; the figure/table to script
map is in its header. The Monte Carlo scripts take from minutes to hours each.
Each run appends a pass/fail line per demo to `figures/run_all.log`.

## Citation

If you use this code, please cite the tutorial (see `CITATION.cff`). If you use
the Allen data subset, also cite the Allen Brain Observatory (see
`data/README.md`).

## License

The code is released under the MIT License (`LICENSE`). The Allen data in
`data/allen/` are not covered by this license; they remain subject to the Allen
Institute's Terms of Use (see `data/README.md`).
