%% hawkes_mc.m  --  Monte Carlo calibration of the asymptotic stats on Hawkes
% Runs M independent Hawkes realizations, computes iPDC, the F-form
% threshold and the Fisher-z CI on each, and accumulates the per-frequency
% type-I rate on null edges and the CI coverage on real edges.

clear; close all
root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(1)

%% Ground truth (matches hawkes_demo.m)
N     = 3;
mu    = [0.5; 0.5; 0.5];
alpha = [0.30  0.00  0.00;
         0.40  0.30  0.00;
         0.00  0.40  0.30];
beta  = 1.0;

T_max       = 5000;
dt_bin      = 0.2;
fs          = 1 / dt_bin;
segLen      = 1024;
overlapFrac = 0.5;
nfft        = segLen;
NW = 4; K = 2*NW - 1;
alpha_sig = 0.05;

%% Theoretical iPDC and edge classification
[I_th, f] = ipdc_hawkes_theory(mu, alpha, beta, fs, nfft);
keep = f <= fs/2;

null_mask = false(N, N);
real_mask = false(N, N);
for i = 1:N
    for j = 1:N
        if i == j, continue; end      % self-terms (saturated near 1) are not assessed
        if max(I_th(keep, i, j)) < 1e-6
            null_mask(i, j) = true;
        else
            real_mask(i, j) = true;
        end
    end
end
fprintf('Edges: %d null, %d real\n', sum(null_mask(:)), sum(real_mask(:)));

%% MC loop
M = 100;
exceed = zeros(nfft, N, N);
cover  = zeros(nfft, N, N);

fprintf('\nMC: M=%d, T_max=%g s, dt_bin=%g s\n', M, T_max, dt_bin);
tic
for m = 1:M
    events = sim_hawkes_mv(mu, alpha, beta, T_max);
    edges = 0:dt_bin:T_max;
    Tb    = numel(edges) - 1;
    X = zeros(N, Tb);
    for i = 1:N
        X(i, :) = histcounts(events{i}, edges);
    end
    X = X - mean(X, 2);

    [Smat, fW, mt] = multitaper_csd(X, fs, NW, K, segLen, overlapFrac, nfft);
    [H, Sg]       = wilson_factor_normalized(Smat);
    I             = ipdc(H, Sg);
    [thr, lo, hi] = ipdc_asymp(H, Sg, mt.effDOF, alpha_sig);

    exceed = exceed + double(I > thr);
    cover  = cover  + double(I_th >= lo & I_th <= hi);

    if mod(m, 10) == 0
        fprintf('  %3d / %d   (elapsed %.0f s)\n', m, M, toc);
    end
end

t1 = exceed / M;
cv = cover  / M;

%% Per-frequency aggregation (drop DC and Nyquist)
fk = find(keep);
fk = fk(fk > 1 & fk < floor(nfft/2) + 1);

t1_per_freq = zeros(numel(fk), 1);
cv_per_freq = zeros(numel(fk), 1);
for kk = 1:numel(fk)
    k = fk(kk);
    t1_per_freq(kk) = mean(t1(k, null_mask), 'all');
    cv_per_freq(kk) = mean(cv(k, real_mask), 'all');
end

fprintf('\nMC calibration (target type-I = %.2f, coverage = %.2f)\n', ...
        alpha_sig, 1 - alpha_sig);
fprintf('  type-I = %.3f   coverage = %.3f\n', mean(t1_per_freq), mean(cv_per_freq));

%% Plot
figure('Color', 'w', 'Position', [80 80 1100 380])
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile
plot(f(fk), t1_per_freq, 'k-', 'LineWidth', 1.2); hold on
yline(alpha_sig, 'r--', sprintf('\\alpha=%.2f', alpha_sig));
ylim([0 max(0.2, 1.3 * max(t1_per_freq))])
xlim([0 fs/2]); grid on
xlabel('Hz'); ylabel('empirical type-I'); title('Null edges')

nexttile
plot(f(fk), cv_per_freq, 'k-', 'LineWidth', 1.2); hold on
yline(1 - alpha_sig, 'r--', sprintf('%.0f%%', 100*(1 - alpha_sig)));
ylim([0 1])
xlim([0 fs/2]); grid on
xlabel('Hz'); ylabel('empirical coverage'); title('Real edges')
sgtitle(sprintf('Hawkes asymptotic-stats calibration (M=%d, L=%d)', ...
                M, mt.effDOF))

style_paper_figure(gcf);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'figures/hawkes_calibration.png', '-dpng', '-r150')
fprintf('SAVED_FIGURE: figures/hawkes_calibration.png\n');
