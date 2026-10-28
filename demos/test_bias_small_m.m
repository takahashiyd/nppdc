%% test_bias_small_m.m -- raw vs bias-corrected iPDC at small m (Table tab:bias_small_m)
% Spectrogram-sized windows; chain 1 -> 2 -> 3 with A_1(2,1) = 0.60 (same 1 -> 2
% iPDC as the top point of validate_bias_correction.m).
root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(0)

fs = 200; nfft = 256;
NW = 4; K = 2*NW - 1; overlap = 0.5;
N = 3;
M_mc = 300;
A21 = 0.60;

A = zeros(3,3,2);
A(:,:,1) = [0.55 0 0; 0 0.30 0; 0 0.50 0.40]; A(2,1,1) = A21;
A(:,:,2) = [-0.25 0 0; 0 -0.10 0; 0 0 -0.20];
Sigma = diag([1.0, 0.8, 1.2]);

[I_true, f] = ipdc_from_var(A, Sigma, fs, nfft);
[peak_val, peak_k] = max(squeeze(I_true(:,2,1).^2));

% Spectrogram-window-sized configurations.
cfgs = struct('T',      {256, 512, 1024, 2048}, ...
              'segLen', {128, 128, 128,  256});   % segLen must be <= nfft (256)

fprintf('truth |i|^2 at peak = %.4f  (f = %.1f Hz)\n', peak_val, f(peak_k));
fprintf('\n%-5s %-7s %-3s %-3s %-8s %-8s %-9s %-9s\n', ...
        'T', 'segLen', 'L', 'm', 'raw', 'closed', 'RAWbias', 'CLbias');
for c = 1:numel(cfgs)
    T = cfgs(c).T; segLen = cfgs(c).segLen;
    raw = zeros(M_mc,1); cf = zeros(M_mc,1);
    skipped = 0;
    for r = 1:M_mc
        X = sim_var(T, A, Sigma);
        try
            [Smt, ~, mt] = multitaper_csd(X, fs, NW, K, segLen, overlap, nfft);
            [H, Sg] = wilson_factor_normalized(Smt);
            I_est = ipdc(H, Sg);
            I_cf = ipdc_bias_correct(I_est, mt.effDOF, N);
            raw(r) = I_est(peak_k,2,1)^2;
            cf(r)  = I_cf (peak_k,2,1)^2;
        catch
            skipped = skipped + 1;
            raw(r) = NaN; cf(r) = NaN;
        end
    end
    raw = raw(~isnan(raw)); cf = cf(~isnan(cf));
    L = mt.effDOF; m = L - N;
    fprintf('%-5d %-7d %-3d %-3d %-8.4f %-8.4f %+.4f   %+.4f\n', ...
            T, segLen, L, m, mean(raw), mean(cf), ...
            mean(raw) - peak_val, mean(cf) - peak_val);
    if skipped > 0
        fprintf('  (%d/%d reps skipped due to Wilson/MT failure)\n', skipped, M_mc);
    end
end
