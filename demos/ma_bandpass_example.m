%% ma_bandpass_example.m
% Invertible-MA DGP with bandpass coupling 1 -> 2, illustrating
%   (a) frequency-localized iPDC structure from an MA-only DGP,
%   (b) the bias correction across frequency (averaged over MC).
%
% DGP (3-channel):
%   x_1(t) = e_1(t)
%   x_2(t) = sum_k b_k e_1(t-k) + e_2(t)
%   x_3(t) = e_3(t)
% b chosen so |B(omega)|^2 peaks in a clear bandpass band.

clear; close all; root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(1)

fs   = 200;
T    = 4096;
NW   = 4; K = 2*NW - 1; segLen = 512; overlap = 0.5; nfft = 512;
N    = 3;
M_mc = 100;       % MC reps for bias-vs-correction averaging

% Sharp bandpass MA: B(z) = z^{-1} * (1 - r^2 z^{-2})^p,
% i.e., p stacked zero-pairs at z = +-r inside the unit circle.
% The squared magnitude is |1 - r^2 e^{-2i omega}|^{2p}, peaking at
% omega = pi/2 (50 Hz at fs=200) with transition narrowness set by p.
% b_0 = 0 (leading delay) keeps the MA polynomial in the Wilson canonical
% form so the recovered iPDC matches the analytic iPDC of B.
r = 0.97; p = 5;
b_kernel = 1;
for j = 1:p, b_kernel = conv(b_kernel, [1, 0, -r^2]); end   % (1 - r^2 z^-2)^p
b = [0, b_kernel];   % prepend delay so b_0 = 0
Sigma_e = diag([1.0, 0.8, 1.2]);

% --- Parametric truth: spectrum and iPDC from the MA polynomial
f = (0:nfft-1)' * fs/nfft;
omega = 2*pi*f/fs;
B_w = zeros(nfft, 1);
for k = 1:numel(b)
    B_w = B_w + b(k) * exp(-1i*(k-1)*omega);
end
S_true = zeros(nfft, N, N);
for k = 1:nfft
    Mk = [ 1,        0, 0;
           B_w(k),   1, 0;
           0,        0, 1 ];
    S_true(k,:,:) = Mk * Sigma_e * Mk';
end
[H_true, Sg_true] = wilson_factor_normalized(S_true);
I_true = ipdc(H_true, Sg_true);

% --- MC sweep: raw and bias-corrected iPDC at every frequency
I_raw_acc  = zeros(M_mc, nfft, N, N);
I_corr_acc = zeros(M_mc, nfft, N, N);
dof_eff = NaN;
for m = 1:M_mc
    nb = numel(b);
    T_sim = T + nb;
    E = randn(T_sim, N) * chol(Sigma_e);
    xs = zeros(T_sim, N);
    for t = nb:T_sim
        xs(t, 1) = E(t, 1);
        xs(t, 2) = E(t, 2);
        for k = 1:nb
            xs(t, 2) = xs(t, 2) + b(k) * E(t-k+1, 1);
        end
        xs(t, 3) = E(t, 3);
    end
    X = xs(nb+1:end, :)';
    [Sh, ~, mt] = multitaper_csd(X, fs, NW, K, segLen, overlap, nfft);
    [H_hat, Sg_hat] = wilson_factor_normalized(Sh);
    I_raw  = ipdc(H_hat, Sg_hat);
    [~, I2_corr] = ipdc_bias_correct(I_raw, mt.effDOF, N);
    I_raw_acc(m,:,:,:)  = I_raw.^2;      % squared modulus
    I_corr_acc(m,:,:,:) = I2_corr;
    dof_eff = mt.effDOF;
end
I_raw_m2  = squeeze(mean(I_raw_acc,  1));   % MC mean of |iota_hat|^2
I_corr_m2 = squeeze(mean(I_corr_acc, 1));   % MC mean of corrected |iota|^2

alpha = 0.05;
[thr_mod, ~, ~] = ipdc_asymp(H_hat, Sg_hat, dof_eff, alpha);
thr_val = thr_mod(1,1,1);

% --- Plot: bandpass edge (2 <- 1) + a null edge (1 <- 2)
keep = 1:floor(nfft/2);
fkeep = f(keep);

fig = figure('Color','w','Position',[60 60 1100 440]);
tl = tiledlayout(1, 2, 'Padding','compact','TileSpacing','compact');

c_true = paper_palette('black');
c_raw  = paper_palette('red');
c_corr = paper_palette('blue');
c_thr  = [0.45 0.45 0.45];

nexttile
plot(fkeep, I_true(keep, 2, 1).^2, '-',  'Color', c_true, 'LineWidth', 2.2); hold on
plot(fkeep, I_raw_m2(keep,  2, 1), '-', 'Color', c_raw,  'LineWidth', 1.5);
plot(fkeep, I_corr_m2(keep, 2, 1), '-', 'Color', c_corr, 'LineWidth', 1.5);
yline(thr_val^2, '--', 'Color', c_thr, 'LineWidth', 1.0);
xlabel('frequency (Hz)'); ylabel('mean |\iota|^2 over MC')
title('Edge 2 \leftarrow 1: bandpass coupling', 'Interpreter','tex')
legend({'truth', 'raw mean', 'bias-corrected mean', sprintf('thr  \\alpha=%.2f', alpha)}, ...
       'Location','south', 'Interpreter','tex')
ymax = max([I_true(keep,2,1).^2; I_raw_m2(keep,2,1)])*1.10;
ylim([0 ymax])

nexttile
plot(fkeep, I_true(keep, 1, 2).^2, '-',  'Color', c_true, 'LineWidth', 2.2); hold on
plot(fkeep, I_raw_m2(keep,  1, 2), '-', 'Color', c_raw,  'LineWidth', 1.5);
plot(fkeep, I_corr_m2(keep, 1, 2), '-', 'Color', c_corr, 'LineWidth', 1.5);
yline(thr_val^2, '--', 'Color', c_thr, 'LineWidth', 1.0);
xlabel('frequency (Hz)'); ylabel('mean |\iota|^2 over MC')
title('Null edge 1 \leftarrow 2 (no true coupling)', 'Interpreter','tex')
legend({'truth (= 0)', 'raw mean', 'bias-corrected mean', sprintf('thr  \\alpha=%.2f', alpha)}, ...
       'Location','northeast', 'Interpreter','tex')
ymax2 = max(I_raw_m2(keep,1,2)) * 1.4;
ylim([0 ymax2])

sgtitle(sprintf('Bandpass MA DGP  (T=%d, L=%d, m=%d, MC=%d)', ...
                T, dof_eff, dof_eff - N, M_mc));

style_paper_figure(fig);

if ~exist('figures', 'dir'), mkdir('figures'); end
set(fig, 'PaperPositionMode','auto');
print(fig, 'figures/ma_bandpass.png', '-dpng', '-r150');
fprintf('SAVED figures/ma_bandpass.png\n');

% --- Numeric summary
fprintf('\n=== Bandpass MA example  (T=%d, m=%d, MC=%d) ===\n', T, dof_eff - N, M_mc);
[~, kp] = max(I_true(keep,2,1));
fprintf('Edge 2 <- 1 (bandpass):\n');
fprintf('  peak f               : %.1f Hz\n', fkeep(kp));
fprintf('  truth at peak f      : %.4f\n', I_true(kp,2,1)^2);
fprintf('  raw mean at peak f   : %.4f\n', I_raw_m2(kp,2,1));
fprintf('  corrected at peak f  : %.4f\n', I_corr_m2(kp,2,1));
fprintf('Null edge 1 <- 2:\n');
kin = keep(2:end);   % exclude DC (real-valued there, different null law)
fprintf('  raw mean (avg over 0<f<fs/2): %.4f  (theory null mean 1/(m+1) = %.4f)\n', ...
        mean(I_raw_m2(kin,1,2)), 1/(dof_eff - N + 1));
fprintf('  corrected mean              : %.4f\n', mean(I_corr_m2(kin,1,2)));
