%% Preconditioning stress test
% Two stress cases that exercise the coherency-matrix preconditioner
% (wilson_factor_normalized). Each case runs baseline vs. coh-precond on
% M_mc independent realizations and reports, averaged over replicates:
%   type-I  : rejection rate of the F-form threshold on the null edges
%             (per frequency, 0 < f < fs/2; target alpha = 0.05)
%   RMSE    : iPDC RMSE against the parametric truth (all edges)
%   SigErr  : relative Frobenius error of Sigma_hat
%   rho     : factorization residual ||S - psi psi'||_F / ||S||_F
% The overlay figures use the first realization only.
%
%   Case A: channel variances span 8 orders of magnitude
%           -> stresses Wilson's matrix conditioning
%   Case B: channel-1 has an AR(2) self-pole at r = 0.97 near f0 = 15 Hz
%           -> sharp spectral peak, narrower than the multitaper bandwidth

clear; close all
root = fileparts(fileparts(mfilename('fullpath')));   % project root (code/)
addpath(fullfile(root, 'src'), fullfile(root, 'demos')); cd(root);   % data/ and figures/ are relative to root
rng(7)

N = 3;
fs = 200;
T  = 4096;
segLen      = 512;
overlapFrac = 0.5;
nfft        = segLen;
NW = 4;  K = 2*NW - 1;
M_mc = 100; alpha = 0.05;

%% Case A. Extreme channel-scale mismatch
A1 = [ 0.55  0.00  0.00;   0.40  0.30  0.00;   0.00  0.50  0.40];
A2 = [-0.25  0.00  0.00;   0.00 -0.10  0.00;   0.00  0.00 -0.20];
A_A = cat(3, A1, A2);
Sigma_A = diag([1e-4, 1, 1e4]);                 % var ratio 1e8

fprintf('=== Case A: channel-scale mismatch  Sigma_diag = [%.0e, %.0e, %.0e] ===\n', ...
        diag(Sigma_A));
[I_param, f] = ipdc_from_var(A_A, Sigma_A, fs, nfft);

resA = mc_metrics(A_A, Sigma_A, I_param, f, T, fs, NW, K, segLen, overlapFrac, nfft, M_mc, alpha);
save_overlay(resA, I_param, f, fs, 'figures/stress_A_scale_mismatch.png', ...
             'Case A: channel-scale mismatch (Var ratio 1e8)');

%% Case B. Near-unit-circle self-pole on channel 1
r  = 0.97;
f0 = 15;                                        % Hz
ph = 2*pi*f0/fs;
A1_B = [ 2*r*cos(ph)  0.00  0.00;   0.40  0.30  0.00;   0.00  0.50  0.40];
A2_B = [-r^2          0.00  0.00;   0.00 -0.10  0.00;   0.00  0.00 -0.20];
A_B  = cat(3, A1_B, A2_B);
Sigma_B = diag([1.0, 0.8, 1.2]);

companion = [reshape(A_B, N, N*2);   eye(N), zeros(N, N)];
assert(max(abs(eig(companion))) < 1, 'Case B VAR is not stable')

fprintf('\n=== Case B: sharp peak  pole r=%.2f at f0=%g Hz on ch1 ===\n', r, f0);
fprintf('     peak width ~ f0/Q ~ %.2f Hz   vs multitaper bw = %.2f Hz\n', ...
        f0/(1/(2*(1-r))), NW * fs / segLen);
[I_param, f] = ipdc_from_var(A_B, Sigma_B, fs, nfft);

resB = mc_metrics(A_B, Sigma_B, I_param, f, T, fs, NW, K, segLen, overlapFrac, nfft, M_mc, alpha);
save_overlay(resB, I_param, f, fs, 'figures/stress_B_peaky.png', ...
             sprintf('Case B: sharp pole r=%.2f at %g Hz', r, f0));

%% --- Local helpers ------------------------------------------------------
function res1 = mc_metrics(A, Sig, I_param, f, T, fs, NW, K, segLen, overlapFrac, nfft, M_mc, alpha)
    N = size(A, 1);
    fin = f > 0 & f < fs/2;
    nullm = false(N); for i=1:N, for j=1:N, nullm(i,j) = i~=j && max(I_param(:,i,j)) < 1e-8; end, end
    acc = zeros(2, 4);   % [typeI rmse sigErr rho] per variant
    for r = 1:M_mc
        X = sim_var(T, A, Sig);
        [res, dof] = run_variants(X, fs, NW, K, segLen, overlapFrac, nfft);
        if r == 1, res1 = res; end
        m = dof - N; Fc = finv(1 - alpha, 2, 2*m); thr = sqrt(Fc / (m + Fc));
        for v = 1:2
            I = res(v).I;
            rej = I(fin, :, :) > thr;  rej = rej(:, nullm);
            d = I(fin,:,:) - I_param(fin,:,:);
            acc(v,:) = acc(v,:) + [mean(rej(:)), sqrt(mean(d(:).^2)), ...
                norm(res(v).Sigma - Sig, 'fro') / norm(Sig, 'fro'), res(v).rho];
        end
    end
    acc = acc / M_mc;
    fprintf('\n%-14s | %-7s | %-9s | %-8s | %-9s   (M=%d, %d null edges)\n', ...
        'variant', 'type-I', 'iPDC RMSE', 'SigErr', 'rho', M_mc, nnz(nullm));
    fprintf('%s\n', repmat('-', 1, 62));
    for v = 1:2
        fprintf('%-14s | %7.3f | %9.3f | %8.3f | %9.2e\n', res1(v).name, acc(v,:));
    end
end
function [res, dof] = run_variants(X, fs, NW, K, segLen, overlapFrac, nfft)
    % (A) baseline
    [S, ~, mt] = multitaper_csd(X, fs, NW, K, segLen, overlapFrac, nfft);
    dof = mt.effDOF;
    nf = size(S, 1); Snrm = 0;
    for k = 1:nf, Snrm = Snrm + norm(squeeze(S(k,:,:)), 'fro'); end
    Snrm = Snrm / nf;
    [Ha, Sa, ~, ia] = wilson_factor(S);
    res(1).name = 'baseline';        res(1).I = ipdc(Ha, Sa);
    res(1).Sigma = Sa;               res(1).info = ia;   res(1).rho = ia.factor_err / Snrm;

    % (B) +coh-precond
    [Hc, Sc, ~, ic] = wilson_factor_normalized(S);
    res(2).name = '+coh-precond';    res(2).I = ipdc(Hc, Sc);
    res(2).Sigma = Sc;               res(2).info = ic;   res(2).rho = ic.factor_err / Snrm;
end

function save_overlay(res, I_param, f, fs, outfile, ttl)
    N = size(I_param, 2);
    keep = f <= fs/2;
    figure('Color','w','Position',[80 80 1200 820]);
    tl = tiledlayout(N, N, 'Padding','compact','TileSpacing','compact');
    cols = {'r','b','g','m'};
    names = {res.name};
    for i = 1:N
        for j = 1:N
            nexttile
            plot(f(keep), I_param(keep, i, j), 'k-', 'LineWidth', 2.0); hold on
            for v = 1:numel(res)
                plot(f(keep), res(v).I(keep, i, j), cols{v}, 'LineWidth', 0.9);
            end
            ylim([0 1]); xlim([0 fs/2]); grid on
            if i == j, title(sprintf('self %d', i))
            else,      title(sprintf('%d \\rightarrow %d', j, i))
            end
            if i == N, xlabel('Hz'); end
            if j == 1, ylabel('iPDC'); end
        end
    end
    title(tl, ttl, 'FontWeight','normal');
    legend(['truth', names], 'Location','northoutside', 'Orientation','horizontal');
    set(gcf, 'PaperPositionMode','auto');
    print(gcf, outfile, '-dpng', '-r110');
    fprintf('SAVED_FIGURE: %s\n', outfile);
end
