function val = ipdc_validate_mc(M, T, A, Sigma_true, fs, NW, K, segLen, overlapFrac, nfft, alpha)
%IPDC_VALIDATE_MC  Monte Carlo calibration of the asymptotic iPDC stats.
%
%   val = ipdc_validate_mc(M, T, A, Sigma, fs, NW, K, segLen, overlapFrac, nfft, alpha)
%
%   Simulates M independent realizations of the VAR(p) defined by (A, Sigma),
%   runs the recommended nonparametric pipeline on each (multitaper ->
%   coherency-preconditioned Wilson -> iPDC), and accumulates the empirical:
%
%     - per-frequency rejection rate of H_0 (should be ~ alpha on null edges)
%     - per-frequency coverage of the 1-alpha CI of the true |iota_ij(f)|
%
%   The rejection rate is the type-I rate on null edges and the power on
%   real edges. Coverage is informative on real edges only: on null edges
%   the CI lower bound is clipped at 0.
%
%   Inputs match nppdc_tutorial.m semantics.
%
%   Output struct fields
%     val.f               nfft x 1 freq grid
%     val.I_true          nfft x N x N parametric ground-truth iPDC
%     val.type1           nfft x N x N empirical rejection rate
%     val.coverage        nfft x N x N empirical CI coverage
%     val.mean_estimate   nfft x N x N mean of the M estimates
%     val.nominal_alpha   scalar
%     val.M               scalar (number of replicates)

    [I_true, f] = ipdc_from_var(A, Sigma_true, fs, nfft);
    N           = size(A, 1);

    accept_rej = zeros(nfft, N, N);
    cover      = zeros(nfft, N, N);
    sum_est    = zeros(nfft, N, N);

    for m = 1:M
        X = sim_var(T, A, Sigma_true);

        % Recommended pipeline: multitaper -> coh-precond Wilson -> iPDC
        [S, fW, mt]    = multitaper_csd(X, fs, NW, K, segLen, overlapFrac, nfft);
        [H, Sg]        = wilson_factor_normalized(S);
        I              = ipdc(H, Sg);

        dof    = mt.effDOF;
        [thr_mod, ci_lo, ci_hi] = ipdc_asymp(H, Sg, dof, alpha);

        accept_rej = accept_rej + double(I > thr_mod);
        cover      = cover      + double(I_true >= ci_lo & I_true <= ci_hi);
        sum_est    = sum_est    + I;

        if mod(m, 25) == 0
            fprintf('  MC progress: %d / %d\n', m, M);
        end
    end

    val.f             = f;
    val.I_true        = I_true;
    val.type1         = accept_rej / M;        % per (f, i, j); meaningful on null edges
    val.coverage      = cover     / M;        % per (f, i, j); meaningful on real edges
    val.mean_estimate = sum_est   / M;
    val.nominal_alpha = alpha;
    val.M             = M;
end
