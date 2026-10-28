%% hawkes_demo.m  --  iPDC of a multivariate Hawkes process
% A 3-component Hawkes with directed chain 1 -> 2 -> 3 (no back-edges).
% Sketch of the test:
%
%   1. Simulate the Hawkes by Ogata thinning.
%   2. Bin events into a regular time grid -> count series X.
%   3. Run the recommended nonparametric pipeline (multitaper ->
%      coh-precond Wilson -> iPDC) on the mean-centred counts.
%   4. Compare to the closed-form iPDC computed directly from the
%      Hawkes kernel.
%   5. Add the asymptotic H0 threshold and 95% CI band and check
%      that null edges stay below threshold while real edges are
%      detected.

clear; close all
root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(7)

%% Ground truth
N     = 3;
mu    = [0.5; 0.5; 0.5];
alpha = [0.30  0.00  0.00;
         0.40  0.30  0.00;
         0.00  0.40  0.30];          % chain 1 -> 2 -> 3, plus self-excitation
beta  = 1.0;

assert(max(abs(eig(alpha))) < 1, 'Hawkes is not stable');
lambda_stat = (eye(N) - alpha) \ mu;
fprintf('Stationary rates  : [%.3f %.3f %.3f] events/sec\n', lambda_stat);

%% Simulate
T_max  = 5000;                        % sec
fprintf('Simulating Hawkes for T = %g sec ...\n', T_max);
tic
events = sim_hawkes_mv(mu, alpha, beta, T_max);
fprintf('done (%.1f s).  events per channel: [%d %d %d]\n', ...
        toc, cellfun(@numel, events));

%% Bin events
dt_bin = 0.2;                         % sec  (fs = 5 Hz, Nyquist 2.5 Hz)
fs     = 1 / dt_bin;
edges  = 0:dt_bin:T_max;
T_bins = numel(edges) - 1;
X = zeros(N, T_bins);
for i = 1:N
    X(i, :) = histcounts(events{i}, edges);
end
fprintf('Binning: dt = %.2g sec, fs = %g Hz, %d bins\n', dt_bin, fs, T_bins);
fprintf('Mean counts per bin   : [%.3f %.3f %.3f]\n', mean(X, 2));
fprintf('Variance to mean ratio: [%.3f %.3f %.3f]   (1 = Poisson)\n', ...
        var(X, 0, 2)./mean(X, 2));

% Mean-centre
X = X - mean(X, 2);

%% Recommended pipeline
segLen      = 1024;
overlapFrac = 0.5;
nfft        = segLen;
NW          = 4;
K           = 2*NW - 1;

[S_est, fW, mt] = multitaper_csd(X, fs, NW, K, segLen, overlapFrac, nfft);
[H_est, Sigma_est, ~, info] = wilson_factor_normalized(S_est);
Iota_est      = ipdc(H_est, Sigma_est);

fprintf('\nMultitaper:  K=%d, segLen=%d, overlap=%.0f%%, nfft=%d, effDOF=%d\n', ...
        K, segLen, 100*overlapFrac, nfft, mt.effDOF);
fprintf('Wilson:      converged=%d, iter=%d, ||S - psi*psi''||/nfft = %.2e\n', ...
        info.converged, info.iter, info.factor_err);

%% Closed-form iPDC from the Hawkes kernel
[Iota_th, f_th] = ipdc_hawkes_theory(mu, alpha, beta, fs, nfft);

%% Asymptotic stats (F-form)
alpha_sig = 0.05;
[thr_mod, ci_lo, ci_hi] = ipdc_asymp(H_est, Sigma_est, mt.effDOF, alpha_sig);

%% Plot:  truth vs estimate, with threshold and CI band
keep = fW <= fs/2;
figure('Color', 'w', 'Position', [80 80 1100 820]);
tl = tiledlayout(N, N, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:N
    for j = 1:N
        nexttile
        xb = [fW(keep); flipud(fW(keep))];
        yb = [ci_lo(keep, i, j); flipud(ci_hi(keep, i, j))];
        fill(xb, yb, [1 0.75 0.75], 'EdgeColor', 'none'); hold on
        plot(fW(keep), thr_mod(keep, i, j), 'b--', 'LineWidth', 1.0)
        plot(fW(keep), Iota_est(keep, i, j), 'r-',  'LineWidth', 1.0)
        plot(fW(keep), Iota_th(keep, i, j),  'k-',  'LineWidth', 1.8)
        ylim([0 1]); xlim([0 fs/2]); grid on
        if i == j, title(sprintf('self %d', i))
        else,      title(sprintf('%d \\rightarrow %d', j, i))
        end
        if i == N, xlabel('Hz'); end
        if j == 1, ylabel('iPDC'); end
    end
end
title(tl, sprintf('Hawkes iPDC: theory (black), estimate (red), F-test threshold (blue), 95%% CI (band).  T=%g s, effDOF=%d', ...
      T_max, mt.effDOF), 'FontWeight', 'bold');

style_paper_figure(gcf);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'figures/hawkes_validation.png', '-dpng', '-r150');
fprintf('SAVED_FIGURE: figures/hawkes_validation.png\n');

%% Numerical summary
fprintf('\n%-12s | %-8s | %-8s | %-8s\n', 'edge', 'max_th', 'max_est', 'detect?');
fprintf('%s\n', repmat('-', 1, 50));
for i = 1:N
    for j = 1:N
        if i == j, lbl = sprintf('self %d', i); else, lbl = sprintf('%d->%d', j, i); end
        mt_th  = max(Iota_th(keep, i, j));
        mt_est = max(Iota_est(keep, i, j));
        mt_thr = max(thr_mod(keep, i, j));
        det = mt_est > mt_thr;
        fprintf('%-12s | %-8.3f | %-8.3f | %-8s\n', lbl, mt_th, mt_est, ...
                ternary(det, 'YES', 'no'));
    end
end


function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
