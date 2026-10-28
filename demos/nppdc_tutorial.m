%% Nonparametric Partial Directed Coherence -- Tutorial
% Estimate the information PDC (iPDC) of a multivariate time series
% without fitting a VAR model. Pipeline:
%
%   multitaper CSD -> coherency-preconditioned Wilson factorization -> iPDC
%
% We validate the pipeline against parametric ground truth on a 3-channel
% VAR(2) with a directed chain 1 -> 2 -> 3 (no back-edges).

clear; close all
root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(7)

%% 1. Ground-truth VAR(2)
% A directed chain  1 -> 2 -> 3  with self-loops on every channel.
N = 3;
p = 2;
A_true = zeros(N, N, p);
% Channel 2 carries a sharp self-resonance (pole radius r=0.9 at 40 Hz:
% a1=2*r*cos(2*pi*40/200)=0.556, a2=-r^2=-0.81) so that the 2->3
% transmission is band-limited, while 1->2 stays broadband.
A_true(:,:,1) = [ 0.55,  0.00,  0.00;
                  0.40,  0.556, 0.00;
                  0.00,  0.50,  0.40];
A_true(:,:,2) = [-0.25,  0.00,  0.00;
                  0.00, -0.81,  0.00;
                  0.00,  0.00, -0.20];
Sigma_true = diag([1.0, 0.8, 1.2]);

% Stability check: eigenvalues of the companion matrix must lie inside
% the unit circle.
companion = [reshape(A_true, N, N*p);
             eye(N*(p-1)), zeros(N*(p-1), N)];
assert(max(abs(eig(companion))) < 1, 'VAR is not stable')

%% 2. Simulate
T  = 4096;
fs = 200;                       % Hz
X  = sim_var(T, A_true, Sigma_true);
fprintf('Simulated X: %d channels x %d samples at %g Hz\n', size(X,1), size(X,2), fs)

%% 3. Parametric (ground-truth) iPDC
segLen      = 512;             % 2.56 s windows at 200 Hz
overlapFrac = 0.5;
nfft        = segLen;
[I_param, f] = ipdc_from_var(A_true, Sigma_true, fs, nfft);

%% 4. Multitaper settings
NW = 4;                        % time-bandwidth on each segment
K  = 2*NW - 1;                 % Slepian count
fprintf('Multitaper:  NW=%g, K=%d, segLen=%d, overlap=%.0f%%, nfft=%d\n', ...
        NW, K, segLen, 100*overlapFrac, nfft)

%% 5. Estimate
keep = f <= fs/2;
[S_X, ~, mt]     = multitaper_csd(X, fs, NW, K, segLen, overlapFrac, nfft);
[H_est, Sigma_est, ~, info] = wilson_factor_normalized(S_X);
I_np             = ipdc(H_est, Sigma_est);

fprintf('Wilson:      converged=%d, iter=%d\n', info.converged, info.iter);

% Quick diagnostics
maxE21 = max(I_np(keep, 1, 2));
maxE31 = max(I_np(keep, 1, 3));
sigErr = norm(Sigma_est - Sigma_true, 'fro');
fprintf('Estimate:    max iPDC 2->1 = %.3f, max iPDC 3->1 = %.3f, ||Sigma err|| = %.3f\n', ...
        maxE21, maxE31, sigErr);

%% 6. Overlay parametric vs the nonparametric estimate
figure('Color', 'w', 'Position', [80 80 980 760])
tl = tiledlayout(N, N, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:N
    for j = 1:N
        nexttile
        plot(f(keep), I_param(keep, i, j), 'k-', 'LineWidth', 1.8); hold on
        plot(f(keep), I_np(   keep, i, j), 'r-', 'LineWidth', 1.0)
        ylim([0 1]); xlim([0 fs/2]); grid on
        if i == j
            title(sprintf('self %d', i))
        else
            title(sprintf('%d \\rightarrow %d', j, i))
        end
        if i == N, xlabel('Hz'); end
        if j == 1, ylabel('iPDC'); end
    end
end
title(tl, 'iPDC: theoretical truth (black) vs nonparametric estimate (red)', ...
      'FontWeight', 'normal')
legend({'theoretical', 'nonparametric'}, 'Location', 'best')

%% 7. Asymptotic statistics: per-frequency H_0 threshold and 1-alpha CI
alpha = 0.05;
dof   = mt.effDOF;
[thr_mod, ci_lo, ci_hi] = ipdc_asymp(H_est, Sigma_est, dof, alpha);

fprintf('\nAsymptotic stats:  alpha=%.2f, eff DOF=%d\n', alpha, dof);

figure('Color', 'w', 'Position', [80 80 1100 820])
tl = tiledlayout(N, N, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:N
    for j = 1:N
        nexttile
        xb = [f(keep);                flipud(f(keep))];
        yb = [ci_lo(keep, i, j);      flipud(ci_hi(keep, i, j))];
        fill(xb, yb, [1 0.75 0.75], 'EdgeColor', 'none'); hold on
        plot(f(keep), thr_mod(keep, i, j), 'b--', 'LineWidth', 1.0)
        plot(f(keep), I_np(keep, i, j),  'r-', 'LineWidth', 1.0)
        plot(f(keep), I_param(keep, i, j),'k-','LineWidth', 1.6)
        ylim([0 1]); xlim([0 fs/2]); grid on
        if i == j, title(sprintf('self %d', i))
        else,      title(sprintf('%d \\rightarrow %d', j, i)); end
        if i == N, xlabel('Hz'); end
        if j == 1, ylabel('iPDC'); end
    end
end
title(tl, sprintf('iPDC with asymptotic H_0 threshold (blue dashed) and %.0f%% Fisher-z CI (red band)', 100*(1-alpha)), ...
      'FontWeight', 'normal')
legend({'CI band','H_0 threshold','estimate','truth'}, 'Location', 'best')

%% 8. Monte Carlo validation of the asymptotic stats
M = 200;
fprintf('\nRunning Monte Carlo validation (M=%d).\n', M);
val = ipdc_validate_mc(M, T, A_true, Sigma_true, fs, NW, K, segLen, overlapFrac, nfft, alpha);

% Classify edges by ground truth
keep_idx  = find(keep);
null_mask = false(N, N);
real_mask = false(N, N);
for i = 1:N
    for j = 1:N
        if i == j, continue; end      % self-terms (|iota| ~ 1, saturated) are not assessed
        peakT = max(val.I_true(keep, i, j));
        if peakT < 1e-6, null_mask(i,j) = true;
        else            real_mask(i,j) = true;
        end
    end
end

t1_per_freq = zeros(numel(keep_idx), 1);
cov_per_freq = zeros(numel(keep_idx), 1);
n_null = sum(null_mask(:));
n_real = sum(real_mask(:));
for kk = 1:numel(keep_idx)
    k = keep_idx(kk);
    t1_per_freq(kk) = mean(val.type1(k, null_mask), 'all');
    cov_per_freq(kk) = mean(val.coverage(k, real_mask), 'all');
end

fprintf('\nMC calibration  (target type-I = %.2f, target coverage = %.2f)\n', ...
        alpha, 1-alpha);
fprintf('  Null edges  (%d pairs):   mean type-I  = %.3f   (across freq)\n', ...
        n_null, mean(t1_per_freq));
fprintf('  Real edges  (%d pairs):   mean coverage = %.3f   (across freq)\n', ...
        n_real, mean(cov_per_freq));

figure('Color','w','Position',[80 80 900 360])
tiledlayout(1, 2, 'Padding','compact','TileSpacing','compact');
nexttile
plot(f(keep), t1_per_freq, 'k-', 'LineWidth', 1.2); hold on
yline(alpha, 'r--', 'nominal \alpha')
ylim([0 max(0.3, 1.2*max(t1_per_freq))]); xlim([0 fs/2]); grid on
xlabel('Hz'); ylabel('empirical type-I'); title('null edges')

nexttile
plot(f(keep), cov_per_freq, 'k-', 'LineWidth', 1.2); hold on
yline(1-alpha, 'r--', sprintf('nominal %.0f%%', 100*(1-alpha)))
ylim([0 1]); xlim([0 fs/2]); grid on
xlabel('Hz'); ylabel('empirical coverage'); title('real edges')
sgtitle(sprintf('Asymptotic-stats calibration over M = %d MC replicates', M))
