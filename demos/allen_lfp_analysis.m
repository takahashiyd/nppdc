%% allen_lfp_analysis.m -- iPDC on V1<->LM LFP from Allen Neuropixels
% Loads data/allen/subset_movie.mat, runs the multitaper -> Wilson -> iPDC
% pipeline on the 16-channel V1+LM LFP, contrasts with bivariate analysis,
% and saves a figure to figures/allen_lfp.png.

clear; close all; root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root

S = load('data/allen/subset_movie.mat');

% Build 16-channel matrix (V1 on top, LM on bottom). The two probes have
% slightly different sample counts in the window; truncate to the smaller.
Tmin = min(size(S.V1_lfp, 2), size(S.LM_lfp, 2));
X    = [single(S.V1_lfp(:, 1:Tmin)); single(S.LM_lfp(:, 1:Tmin))];   % 16 x Tmin
fs   = S.fs_lfp_v1;
[N, T] = size(X);
NV1 = size(S.V1_lfp, 1); NLM = size(S.LM_lfp, 1);
fprintf('LFP: N=%d (%d V1 + %d LM), T=%d at fs=%.2f Hz (%.1f s)\n', ...
        N, NV1, NLM, T, fs, T/fs);

% Highpass at 1 Hz to remove DC drift
[bH, aH] = butter(4, 1/(fs/2), 'high');
X = filtfilt(bH, aH, double(X.')).';
% Notch at 60 Hz (US line) -- a couple of channels show a small peak
[bN, aN] = iirnotch(60/(fs/2), 60/(fs/2)/30);
X = filtfilt(bN, aN, X.').';
% Center
X = X - mean(X, 2);

%% Multitaper -> Wilson -> iPDC
segLen      = 2048;       % 1.64 s segments at 1250 Hz
overlapFrac = 0.5;
NW          = 4;
K           = 2*NW - 1;
nfft        = segLen;

fprintf('Multitaper segLen=%d, K=%d ...\n', segLen, K);
[S_mt, fW, mt] = multitaper_csd(X, fs, NW, K, segLen, overlapFrac, nfft);
fprintf('  effDOF = %d\n', mt.effDOF);

fprintf('Wilson factorization ...\n');
[H, Sg, ~, info] = wilson_factor_normalized(S_mt);
fprintf('  converged=%d, iter=%d, ||S - psi*psi''||/nfft = %.2e\n', ...
        info.converged, info.iter, info.factor_err);

fprintf('iPDC ...\n');
I = ipdc(H, Sg);

%% Asymptotic threshold + CI (use F-form, alpha=0.05)
alpha_sig = 0.05;
[thr_mod, ci_lo, ci_hi] = ipdc_asymp(H, Sg, mt.effDOF, alpha_sig);

% Restrict to the analyzable band (drop DC and Nyquist tiny edges)
keep = (fW >= 2) & (fW <= 120);

% Aggregate iPDC by source/target area:
%   V1->LM = mean over (i in LM rows, j in V1 cols)
%   LM->V1 = mean over (i in V1 rows, j in LM cols)
v1_idx = 1:NV1; lm_idx = NV1+1 : NV1+NLM;

% Squared iPDC, easier to average
I2 = I.^2;

% (i, j) = target i, source j  in our convention
v1_to_lm = squeeze(mean(mean(I2(:, lm_idx, v1_idx), 2), 3));   % nfft x 1
lm_to_v1 = squeeze(mean(mean(I2(:, v1_idx, lm_idx), 2), 3));

% Also: within-area mean (excluding self-loops)
within_v1 = zeros(nfft, 1); within_lm = zeros(nfft, 1);
for k = 1:nfft
    Ik = squeeze(I2(k, :, :));
    Av = Ik(v1_idx, v1_idx); Av(eye(NV1)==1) = NaN; within_v1(k) = nanmean(Av(:));
    Al = Ik(lm_idx, lm_idx); Al(eye(NLM)==1) = NaN; within_lm(k) = nanmean(Al(:));
end

thr2 = (thr_mod(:,1,1).^2);   % constant across (i,j,f); just take one column

%% Bivariate comparison on a representative pair: midlayer V1 channel vs midlayer LM
midV1 = round(NV1/2); midLM = round(NLM/2);
X_pair = X([midV1, NV1 + midLM], :);
[Sp, fp, mtp] = multitaper_csd(X_pair, fs, NW, K, segLen, overlapFrac, nfft);
[Hp, Sgp]     = wilson_factor_normalized(Sp);
Ip            = ipdc(Hp, Sgp);
biv_v1_to_lm  = squeeze(Ip(:, 2, 1).^2);     % source=col 1 = V1
biv_lm_to_v1  = squeeze(Ip(:, 1, 2).^2);
coh_bivar     = abs(Sp(:,1,2)).^2 ./ max(real(Sp(:,1,1)) .* real(Sp(:,2,2)), eps);

% bivariate H0 threshold (N = 2 channels, same L)
[thr_p, ~, ~] = ipdc_asymp(Hp, Sgp, mtp.effDOF, alpha_sig);
thr2_biv = thr_p(:,1,1).^2;

%% Plot -- all four panels share one logarithmic y-axis so that the
% multivariate (~1e-3), within-area (~2e-2) and bivariate (~1e-1) values
% can be compared at a glance
yl = [1e-4 1];
fig = figure('Color','w','Position',[60 60 1200 700]);
tl  = tiledlayout(2, 2, 'TileSpacing','compact','Padding','compact');

nexttile
semilogy(fW(keep), v1_to_lm(keep), 'r', 'LineWidth', 1.6); hold on
semilogy(fW(keep), lm_to_v1(keep), 'b', 'LineWidth', 1.6);
semilogy(fW(keep), thr2(keep),     'k--');
grid on; xlim([2 120]); ylim(yl)
xlabel('Hz'); ylabel('|iPDC|^2')
title('Multivariate iPDC: cross-area mean')
legend({'V1\rightarrowLM','LM\rightarrowV1','H_0 thr (\alpha=0.05)'}, 'Location','northeast')

nexttile
semilogy(fW(keep), within_v1(keep), 'r', 'LineWidth', 1.4); hold on
semilogy(fW(keep), within_lm(keep), 'b', 'LineWidth', 1.4);
semilogy(fW(keep), thr2(keep),      'k--');
grid on; xlim([2 120]); ylim(yl)
xlabel('Hz'); ylabel('|iPDC|^2')
title('Multivariate iPDC: within-area mean')
legend({'V1\rightarrowV1','LM\rightarrowLM','H_0 thr'}, 'Location','northeast')

nexttile
semilogy(fp(keep), biv_v1_to_lm(keep), 'r', 'LineWidth', 1.6); hold on
semilogy(fp(keep), biv_lm_to_v1(keep), 'b', 'LineWidth', 1.6);
semilogy(fp(keep), thr2_biv(keep),     'k--');
grid on; xlim([2 120]); ylim(yl)
xlabel('Hz'); ylabel('|iPDC|^2'); title('Bivariate iPDC: one mid-layer pair')
legend({'V1\rightarrowLM','LM\rightarrowV1','H_0 thr'}, 'Location','northeast')

nexttile
semilogy(fp(keep), coh_bivar(keep), 'k', 'LineWidth', 1.4); hold on
grid on; xlim([2 120]); ylim(yl)
xlabel('Hz'); ylabel('|C_{12}(f)|^2')
title('Ordinary coherence: same pair')

if ~exist('figures', 'dir'), mkdir('figures'); end
style_paper_figure(fig);
set(fig, 'PaperPositionMode','auto');
print(fig, 'figures/allen_lfp.png', '-dpng', '-r150');
fprintf('SAVED figures/allen_lfp.png\n');

%% Console summary at canonical bands
gamma_band = fW >= 30 & fW <= 80;
beta_band  = fW >= 12 & fW <= 30;
alpha_band = fW >= 8  & fW <= 12;
fprintf('\nMean |iPDC|^2 by band (V1->LM / LM->V1):\n');
fprintf('  alpha (8-12 Hz):  %.4f / %.4f   ratio = %.2f\n', mean(v1_to_lm(alpha_band)), mean(lm_to_v1(alpha_band)), mean(v1_to_lm(alpha_band))/mean(lm_to_v1(alpha_band)));
fprintf('  beta  (12-30 Hz): %.4f / %.4f   ratio = %.2f\n', mean(v1_to_lm(beta_band)),  mean(lm_to_v1(beta_band)),  mean(v1_to_lm(beta_band))/mean(lm_to_v1(beta_band)));
fprintf('  gamma (30-80 Hz): %.4f / %.4f   ratio = %.2f\n', mean(v1_to_lm(gamma_band)), mean(lm_to_v1(gamma_band)), mean(v1_to_lm(gamma_band))/mean(lm_to_v1(gamma_band)));
fprintf('  threshold |iPDC|^2 = %.4f\n', thr2(round(end/2)));
