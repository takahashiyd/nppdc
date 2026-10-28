%% validate_bias_correction.m -- MC validation of iPDC bias correction
% Sweep true |iota|^2 by varying a single coupling coefficient.
% Compare raw and bias-corrected |iPDC|^2 expectations to truth.

clear; close all; root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(0)

fs = 200;  T = 4096;  nfft = 512;
NW = 4;  K = 2*NW - 1;  segLen = 512;  overlap = 0.5;
M  = 300;        % MC replicates per coupling level
N  = 3;

% Reference VAR(2): sweep A_1(2,1) (the 1->2 coupling at lag 1)
A_base = zeros(3,3,2);
% Channel 2 resonant (r=0.9 @ 40 Hz) to match the reference DGP; the swept
% 1->2 edge (i=2,j=1) is independent of channel 2's self-dynamics, so the
% bias results here are identical to the non-resonant case.
A_base(:,:,1) = [0.55 0 0; 0 0.556 0; 0 0.50 0.40];
A_base(:,:,2) = [-0.25 0 0; 0 -0.81 0; 0 0 -0.20];
Sigma = diag([1.0, 0.8, 1.2]);

coup_levels = [0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60];

% For each coupling, compute the true |iota|^2 at the peak frequency from
% the parametric formula; then run MC and record E[|i_hat|^2] (raw and
% bias-corrected) at the SAME peak.

results = struct('coup', [], 'true_peak_i2', [], ...
                 'E_raw_at_peak', [], 'E_corr_at_peak', [], 'm_used', []);

for cl = 1:numel(coup_levels)
    A = A_base;
    A(2, 1, 1) = coup_levels(cl);

    % stability check (companion eigenvalues)
    comp = [reshape(A, N, N*2); eye(N), zeros(N, N)];
    if max(abs(eig(comp))) >= 0.999
        fprintf('  level %d: unstable, skipping\n', cl);
        continue
    end

    % parametric truth
    [I_true, f] = ipdc_from_var(A, Sigma, fs, nfft);
    % focus on the i=2, j=1 edge -- this is the one we modulate
    iTrue2_21 = squeeze(I_true(:, 2, 1).^2);
    [peak_val, peak_k] = max(iTrue2_21);
    fprintf('coup A21=%.2f  true |iota|^2 peak = %.4f at f=%.1f Hz\n', ...
            coup_levels(cl), peak_val, f(peak_k));

    raw_acc  = zeros(M, 1);
    corr_acc = zeros(M, 1);
    for m = 1:M
        X = sim_var(T, A, Sigma);
        [Smt, fW, mt]   = multitaper_csd(X, fs, NW, K, segLen, overlap, nfft);
        [H, Sg]         = wilson_factor_normalized(Smt);
        I_est           = ipdc(H, Sg);
        I_corr          = ipdc_bias_correct(I_est, mt.effDOF, N);

        raw_acc(m)      = I_est(peak_k, 2, 1).^2;
        corr_acc(m)     = I_corr(peak_k, 2, 1).^2;
    end

    results(end+1).coup            = coup_levels(cl);
    results(end).true_peak_i2      = peak_val;
    results(end).E_raw_at_peak     = mean(raw_acc);
    results(end).E_corr_at_peak    = mean(corr_acc);
    results(end).std_raw_at_peak   = std(raw_acc);
    results(end).std_corr_at_peak  = std(corr_acc);
    results(end).m_used            = mt.effDOF - N;
end
results = results(2:end);   % drop the empty seed

%% Plot
fig = figure('Color','w','Position',[60 60 1100 440]);
tl = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

true_i2 = [results.true_peak_i2];
E_raw   = [results.E_raw_at_peak];
E_corr  = [results.E_corr_at_peak];
m_used  = results(1).m_used;
predicted_raw = true_i2 + (1 - true_i2).^2 / (m_used + 1);   % Carter et al. (1973), n = m + 1 DOF

c_raw  = paper_palette('red');
c_corr = paper_palette('blue');
c_pred = paper_palette('gray');

nexttile
plot([0, 1], [0, 1], '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9); hold on
plot(true_i2, predicted_raw, '-',  'Color', c_pred, 'LineWidth', 1.2);
plot(true_i2, E_raw,  'o', 'MarkerSize', 7, 'LineWidth', 1.3, ...
     'MarkerEdgeColor', c_raw,  'MarkerFaceColor', 'w');
plot(true_i2, E_corr, 's', 'MarkerSize', 7, 'LineWidth', 1.3, ...
     'MarkerEdgeColor', c_corr, 'MarkerFaceColor', 'w');
xlim([0 max(true_i2)*1.08]); ylim([0 max(E_raw)*1.10])
xlabel('true |\iota|^2'); ylabel('mean estimated |\iota|^2')
L_used = m_used + N;
title(sprintf('Bias correction at the peak  (L=%d, m=%d, M=%d)', ...
              L_used, m_used, M));
legend({'y = x', 'predicted raw', 'raw \langle|\^\iota|^2\rangle', ...
        'corrected'}, 'Location','northwest', 'Interpreter','tex');

nexttile
plot([0, max(true_i2)*1.08], [0 0], '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9); hold on
plot(true_i2, E_raw  - true_i2, '-o', 'Color', c_raw,  'LineWidth', 1.3, ...
     'MarkerSize', 7, 'MarkerFaceColor', 'w');
plot(true_i2, E_corr - true_i2, '-s', 'Color', c_corr, 'LineWidth', 1.3, ...
     'MarkerSize', 7, 'MarkerFaceColor', 'w');
xlim([0 max(true_i2)*1.08])
xlabel('true |\iota|^2'); ylabel('bias  E[\langle|\^\iota|^2\rangle] - |\iota|^2')
title('Residual bias');
legend({'zero', 'raw', 'corrected'}, 'Location','southwest');

style_paper_figure(fig);

if ~exist('figures', 'dir'), mkdir('figures'); end
set(fig, 'PaperPositionMode','auto');
print(fig, 'figures/bias_correction.png', '-dpng', '-r150');
fprintf('SAVED figures/bias_correction.png\n');

fprintf('\nMC summary (m = %d):\n', m_used);
fprintf('  %-10s %-12s %-12s %-12s %-12s %-12s\n', ...
        'A_1(2,1)', 'true |i|^2', 'E[|i_hat|^2]', 'predicted', 'E[|i_corr|^2]', 'bias-corr');
for k = 1:numel(results)
    r = results(k);
    pred = r.true_peak_i2 + (1 - r.true_peak_i2)^2/(m_used + 1);   % Carter, n = m + 1
    fprintf('  %-10.2f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
            r.coup, r.true_peak_i2, r.E_raw_at_peak, pred, ...
            r.E_corr_at_peak, r.E_corr_at_peak - r.true_peak_i2);
end
