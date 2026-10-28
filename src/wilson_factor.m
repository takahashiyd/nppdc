function [H, Sigma, psi, info] = wilson_factor(S, tol, maxiter)
%WILSON_FACTOR  Wilson minimum-phase spectral matrix factorization.
%
%   [H, Sigma, psi, info] = wilson_factor(S, tol, maxiter)
%
%   Given a Hermitian positive-definite cross-spectral density matrix
%   sampled on a uniform two-sided frequency grid, this routine finds the
%   minimum-phase factor psi(f) such that
%
%       S(f) = psi(f) * psi(f)'    (Hermitian transpose).
%
%   It then converts to the (H, Sigma) parameterization
%
%       S(f) = H(f) * Sigma * H(f)',
%
%   which is the form needed for partial directed coherence.
%
%   Inputs
%     S       : nfft x N x N two-sided CSD tensor (nfft preferably even).
%     tol     : relative tolerance on ||psi^(k+1) - psi^(k)||_F (default 1e-9)
%     maxiter : maximum iterations (default 500)
%
%   Outputs
%     H     : nfft x N x N transfer function
%     Sigma : N x N real symmetric positive-definite innovation covariance
%     psi   : nfft x N x N minimum-phase factor (S = psi * psi')
%     info  : struct with fields .iter, .converged, .resid (final relative
%             change of psi), .factor_err (mean Frobenius residual of
%             S - psi*psi' over frequencies).
%
%   References
%     Wilson, G. T. (1972). The factorization of matricial spectral
%       densities. SIAM J. Appl. Math. 23, 420-426.
%     Dhamala, Rangarajan, Ding (2008). Estimating Granger causality from
%       Fourier and wavelet transforms of time series data. PRL 100, 018701.
%     FieldTrip's sfactorization_wilson.m (reference implementation).

    if nargin < 2 || isempty(tol),     tol = 1e-9;  end
    if nargin < 3 || isempty(maxiter), maxiter = 500; end

    [nfft, N, N2] = size(S);
    assert(N == N2, 'S must be square in the last two dimensions');

    % --- Initial guess: psi(f) = chol(Gamma_0) for all f ---------------
    Gamma0 = squeeze(mean(S, 1));
    Gamma0 = real((Gamma0 + Gamma0') / 2);
    [L, pflag] = chol(Gamma0, 'lower');
    if pflag ~= 0
        % Tiny ridge if numerics push us slightly indefinite.
        Gamma0 = Gamma0 + 1e-12 * trace(Gamma0) * eye(N);
        L = chol(Gamma0, 'lower');
    end
    psi = repmat(reshape(L, [1, N, N]), [nfft, 1, 1]);

    I_N = eye(N);
    converged = false;
    resid = NaN;

    % --- Wilson iteration ----------------------------------------------
    for iter = 1:maxiter
        psi_prev = psi;

        % g(f) = psi(f)^{-1} S(f) psi(f)^{-H} + I
        g = zeros(nfft, N, N);
        for k = 1:nfft
            psik = squeeze(psi(k,:,:));
            Sk   = squeeze(S(k,:,:));
            g(k,:,:) = psik \ Sk / psik' + I_N;
        end

        % Causal projection in the time domain
        gt = ifft(g, [], 1);
        g0 = squeeze(gt(1,:,:));
        g0 = real(g0);
        gt(1,:,:) = tril(g0, -1) + 0.5 * diag(diag(g0));
        gt(floor(nfft/2)+1:nfft, :, :) = 0;
        gplus = fft(gt, [], 1);

        % psi <- psi * gplus
        for k = 1:nfft
            psi(k,:,:) = squeeze(psi(k,:,:)) * squeeze(gplus(k,:,:));
        end

        delta = psi(:) - psi_prev(:);
        resid = norm(delta) / max(norm(psi_prev(:)), eps);
        if resid < tol
            converged = true;
            break
        end
    end

    if ~converged
        warning('wilson_factor:noConvergence', ...
            'Wilson factorization did not converge in %d iterations (resid=%.3g).', ...
            maxiter, resid);
    end

    % --- Extract Sigma and H from psi ----------------------------------
    % psi(tau=0) is the lag-0 time-domain coefficient of psi(f),
    % equal to (1/nfft) * sum_k psi(f_k) = mean(psi, 1).
    psi_t0 = squeeze(mean(psi, 1));
    psi_t0 = real(psi_t0);                 % theory: real lower-triangular
    Sigma = psi_t0 * psi_t0.';
    Sigma = (Sigma + Sigma.') / 2;          % symmetrize

    H = zeros(nfft, N, N);
    invPsi0 = inv(psi_t0);
    for k = 1:nfft
        H(k,:,:) = squeeze(psi(k,:,:)) * invPsi0;
    end

    % Diagnostic: residual of the factorization
    fe = 0;
    for k = 1:nfft
        psik = squeeze(psi(k,:,:));
        fe = fe + norm(squeeze(S(k,:,:)) - psik * psik', 'fro');
    end
    info.iter        = iter;
    info.converged   = converged;
    info.resid       = resid;
    info.factor_err  = fe / nfft;
end
