%% more_checks.m  --  deeper validation of the asymptotic iPDC statistics
% Three checks of the asymptotic iPDC statistics on the reference DGP:
%
%   1. alpha-sweep of the F-form type-I rate, alpha in {1e-3, ..., 0.2}
%   2. power curve on the 1 -> 2 edge at alpha = 0.05
%   3. type-I and coverage under t_5 (heavy-tailed) innovations

clear; close all
root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(42)

%% Common setup
N = 3; fs = 200; T = 4096;
segLen = 512; overlapFrac = 0.5; nfft = segLen;
NW = 4; K = 2*NW - 1;

% Channel 2 resonant (r=0.9 @ 40 Hz) -> band-limited 2->3; 1->2 broadband.
A1 = [ 0.55  0.00  0.00;  0.40  0.556 0.00;  0.00  0.50  0.40];
A2 = [-0.25  0.00  0.00;  0.00 -0.81  0.00;  0.00  0.00 -0.20];
A_DGP = cat(3, A1, A2);
Sigma_DGP = diag([1.0, 0.8, 1.2]);

alpha0 = 0.05;

%% Check 1: alpha-sweep calibration (null |iota|^2 collected by MC)
M_qa = 300;
fprintf('=== Null collection for the alpha-sweep (M=%d) ===\n', M_qa);
[null_ipdc2, dof_eff, ~] = mc_null_collect(M_qa, T, A_DGP, Sigma_DGP, fs, NW, K, segLen, overlapFrac, nfft);
alpha_levels = [0.001 0.005 0.01 0.02 0.05 0.10 0.20];
alpha_sweep(null_ipdc2, dof_eff, nfft, alpha_levels, 'figures/check_alpha_sweep.png');

%% Check 2: power curve
couplings = [0.0 0.05 0.10 0.20 0.30 0.40 0.50];
M_pow = 100;
fprintf('\n=== Power curve sweep (M=%d per coupling) ===\n', M_pow);
power_curve(couplings, M_pow, T, A_DGP, Sigma_DGP, fs, NW, K, segLen, overlapFrac, nfft, alpha0, 'figures/check_power.png');

%% Check 3: non-Gaussian innovations (t_5)
M_t = 200;
nu_t = 5;
fprintf('\n=== Non-Gaussian innovations t_%d (M=%d) ===\n', nu_t, M_t);
nongauss_check(M_t, T, A_DGP, Sigma_DGP, nu_t, fs, NW, K, segLen, overlapFrac, nfft, alpha0, 'figures/check_nongauss.png');


%% ---- local helpers -------------------------------------------------------

function [vals, dof_eff, f] = mc_null_collect(M, T, A, Sigma, fs, NW, K, segLen, overlap, nfft)
    N = size(A, 1);
    [I_true, f] = ipdc_from_var(A, Sigma, fs, nfft);
    null_mask = false(N, N);
    for i = 1:N, for j = 1:N
        if max(I_true(:, i, j)) < 1e-6, null_mask(i, j) = true; end
    end; end
    [iIdx, jIdx] = find(null_mask);
    nE = numel(iIdx);

    vals = zeros(M, nfft, nE);
    dof_eff = 0;
    for m = 1:M
        X = sim_var(T, A, Sigma);
        [Smat, fW, mt] = multitaper_csd(X, fs, NW, K, segLen, overlap, nfft);
        [H, Sg]       = wilson_factor_normalized(Smat);
        I = ipdc(H, Sg);
        I2 = I.^2;
        for e = 1:nE
            vals(m, :, e) = I2(:, iIdx(e), jIdx(e));
        end
        dof_eff = mt.effDOF;
        if mod(m, 50) == 0, fprintf('  null collect: %d / %d\n', m, M); end
    end
end

function alpha_sweep(null_ipdc2, dof, nfft, alpha_levels, outfile)
    f_keep = 2:floor(nfft/2);
    T_stat = dof * null_ipdc2(:, f_keep, :);
    flat = T_stat(:);

    % F-form threshold, m = L - N
    N_ch = 3;  m = dof - N_ch;
    emp = zeros(numel(alpha_levels), 1);
    for k = 1:numel(alpha_levels)
        Fc  = finv(1 - alpha_levels(k), 2, 2*m);
        emp(k) = mean(flat > dof * Fc / (m + Fc));
    end

    fprintf('\nAlpha-sweep (DC and Nyquist excluded, %d samples)\n', numel(flat));
    fprintf('  %-10s %-12s %-10s\n', 'nominal', 'empirical', 'ratio');
    for k = 1:numel(alpha_levels)
        fprintf('  %-10.4f %-12.4f %-10.3f\n', alpha_levels(k), emp(k), emp(k)/alpha_levels(k));
    end

    figure('Color','w','Position',[100 100 900 380]);
    tiledlayout(1, 2, 'Padding','compact','TileSpacing','compact');
    nexttile
    loglog(alpha_levels, emp, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor','k'); hold on
    lim = [min(alpha_levels)/2  max(alpha_levels)*2];
    loglog(lim, lim, 'r--', 'LineWidth', 1.0)
    xlabel('nominal \alpha'); ylabel('empirical type-I'); grid on
    title('Calibration (log-log)')
    nexttile
    semilogx(alpha_levels, emp ./ alpha_levels, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor','k'); hold on
    yline(1, 'r--');
    xlabel('nominal \alpha'); ylabel('empirical / nominal'); grid on
    title('Ratio  (1 = perfect)')
    sgtitle('Asymptotic threshold calibration across \alpha levels')
    style_paper_figure(gcf);
    set(gcf,'PaperPositionMode','auto'); print(gcf, outfile, '-dpng', '-r150');
    fprintf('SAVED_FIGURE: %s\n', outfile);
end

function power_curve(couplings, M, T, A_base, Sigma, fs, NW, K, segLen, overlap, nfft, alpha, outfile)
    f_keep = 2:floor(nfft/2);
    nC = numel(couplings);
    power_max = zeros(nC, 1);
    power_avg = zeros(nC, 1);
    true_max  = zeros(nC, 1);

    for c = 1:nC
        A = A_base;
        A(2, 1, 1) = couplings(c);             % the 1 -> 2 coupling
        N = size(A, 1);
        comp = [reshape(A, N, N*2); eye(N), zeros(N, N)];
        if max(abs(eig(comp))) >= 1
            warning('Unstable VAR at coupling=%.2f; skipping', couplings(c)); continue
        end
        [I_true, ~] = ipdc_from_var(A, Sigma, fs, nfft);
        true_max(c) = max(I_true(f_keep, 2, 1));

        n_reject = zeros(numel(f_keep), 1);
        for m = 1:M
            X = sim_var(T, A, Sigma);
            [Smat, fW, mt] = multitaper_csd(X, fs, NW, K, segLen, overlap, nfft);
            [H, Sg]        = wilson_factor_normalized(Smat);
            I             = ipdc(H, Sg);
            [thr, ~, ~]   = ipdc_asymp(H, Sg, mt.effDOF, alpha);
            rej = I(f_keep, 2, 1) > thr(f_keep, 2, 1);
            n_reject = n_reject + rej(:);
        end
        pw = n_reject / M;
        power_max(c) = max(pw);
        power_avg(c) = mean(pw);
        fprintf('  coupling=%.2f  true|iota_21|_max=%.3f  power_max=%.3f  power_avg=%.3f\n', ...
                couplings(c), true_max(c), power_max(c), power_avg(c));
    end

    figure('Color','w','Position',[100 100 760 520]);
    plot(true_max, power_max, 'bo-', 'LineWidth', 1.6, 'MarkerFaceColor','b'); hold on
    plot(true_max, power_avg, 'rs-', 'LineWidth', 1.6, 'MarkerFaceColor','r');
    yline(alpha, 'k--', sprintf('\\alpha=%.2f', alpha), 'LabelHorizontalAlignment','left');
    xlabel('true |\iota_{21}|_{max}'); ylabel('detection probability'); grid on
    ylim([0 1.05])
    title(sprintf('Power curve on 1\\rightarrow2 edge  (M=%d, \\alpha=%.2f, L=%d)', ...
                  M, alpha, mt.effDOF))
    legend({'max over freq','avg over freq','\alpha'}, 'Location','best')
    style_paper_figure(gcf);
    set(gcf,'PaperPositionMode','auto'); print(gcf, outfile, '-dpng', '-r150');
    fprintf('SAVED_FIGURE: %s\n', outfile);
end

function nongauss_check(M, T, A, Sigma, nu_t, fs, NW, K, segLen, overlap, nfft, alpha, outfile)
    N = size(A, 1);
    [I_true, f] = ipdc_from_var(A, Sigma, fs, nfft);
    null_mask = false(N, N); real_mask = false(N, N);
    for i = 1:N, for j = 1:N
        if i == j, continue; end      % self-terms (saturated near 1) are not assessed
        if max(I_true(:, i, j)) < 1e-6, null_mask(i,j) = true;
        else                            real_mask(i,j) = true; end
    end; end

    exceed = zeros(nfft, N, N);
    cover  = zeros(nfft, N, N);
    L = chol(Sigma, 'lower');
    burn = 500;

    for m = 1:M
        Ttot = T + burn;
        % standardized t_{nu_t} innovations with unit variance per channel
        Z = randn(N, Ttot);
        Vch = zeros(N, Ttot);
        for q = 1:nu_t, Vch = Vch + randn(N, Ttot).^2; end
        Ec = Z .* sqrt(nu_t ./ Vch) / sqrt(nu_t/(nu_t-2));
        E = L * Ec;

        X = zeros(N, Ttot);
        for t = 3:Ttot
            X(:, t) = E(:, t) + A(:,:,1)*X(:, t-1) + A(:,:,2)*X(:, t-2);
        end
        X = X(:, burn+1:end);

        [Smat, fW, mt] = multitaper_csd(X, fs, NW, K, segLen, overlap, nfft);
        [H, Sg]       = wilson_factor_normalized(Smat);
        I = ipdc(H, Sg);
        [thr, ci_lo, ci_hi] = ipdc_asymp(H, Sg, mt.effDOF, alpha);

        exceed = exceed + double(I    > thr);
        cover  = cover  + double(I_true >= ci_lo & I_true <= ci_hi);
        if mod(m, 50) == 0, fprintf('  t-MC: %d / %d\n', m, M); end
    end
    type1 = exceed / M;
    cov   = cover  / M;

    f_keep = 2:floor(nfft/2);
    t1 = zeros(numel(f_keep), 1);
    cv = zeros(numel(f_keep), 1);
    for kk = 1:numel(f_keep)
        k = f_keep(kk);
        t1(kk) = mean(type1(k, null_mask), 'all');
        cv(kk) = mean(cov(k, real_mask), 'all');
    end
    fprintf('  empirical type-I (null, mean over freq) = %.3f\n', mean(t1));
    fprintf('  empirical coverage (real, mean over freq) = %.3f\n', mean(cv));

    figure('Color','w','Position',[100 100 1000 360]);
    tiledlayout(1, 2, 'Padding','compact','TileSpacing','compact');
    nexttile
    plot(f(f_keep), t1, 'k-', 'LineWidth', 1.2); hold on
    yline(alpha, 'r--', sprintf('\\alpha=%.2f', alpha));
    ylim([0 max(0.2, 1.3*max(t1))]); xlim([0 fs/2]); grid on
    xlabel('Hz'); ylabel('empirical type-I'); title('null edges (t innovations)')
    nexttile
    plot(f(f_keep), cv, 'k-', 'LineWidth', 1.2); hold on
    yline(1 - alpha, 'r--', sprintf('%.0f%%', 100*(1 - alpha)));
    ylim([0 1.05]); xlim([0 fs/2]); grid on
    xlabel('Hz'); ylabel('empirical coverage'); title('real edges (t innovations)')
    sgtitle(sprintf('Robustness to t_{%d} innovations  (M=%d)', nu_t, M))
    style_paper_figure(gcf);
    set(gcf,'PaperPositionMode','auto'); print(gcf, outfile, '-dpng', '-r150');
    fprintf('SAVED_FIGURE: %s\n', outfile);
end
