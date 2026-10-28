function [H, Sigma, psi, info] = wilson_factor_normalized(S, tol, maxiter)
%WILSON_FACTOR_NORMALIZED  Wilson factorization with coherency preconditioning.
%
%   [H, Sigma, psi, info] = wilson_factor_normalized(S, tol, maxiter)
%
%   Factorizes S(f) = H(f) * Sigma * H(f)' the same way as WILSON_FACTOR,
%   but first normalizes each diagonal of S(f) to 1 by dividing rows and
%   columns by the scalar minimum-phase factor d_i(f) (|d_i|^2 = S_ii):
%
%       R(f) = D(f)^{-1} * S(f) * D(f)^{-H}       (unit diagonal, |R_ij| <= 1)
%       R(f) = psi_R(f) * psi_R(f)'               (Wilson on R)
%       psi_S(f) = D(f) * psi_R(f)                (minimum-phase product)
%
%   This is a standard preconditioner for matricial spectral factorization
%   (e.g., Sayed & Kailath 2001). R(f) is far better conditioned than S(f)
%   when channel variances or spectral peaks differ by orders of magnitude,
%   and Wilson on R typically converges in a handful of iterations.
%
%   Inputs and outputs match WILSON_FACTOR. info.normalized == true.

    if nargin < 2 || isempty(tol),     tol = 1e-9;  end
    if nargin < 3 || isempty(maxiter), maxiter = 500; end

    [nfft, N, N2] = size(S);
    assert(N == N2, 'S must be square in last two dims');

    % --- Scalar min-phase factor per channel ---------------------------
    d = zeros(nfft, N);
    for n = 1:N
        d(:, n) = scalar_min_phase_factor(real(S(:, n, n)));
    end

    % --- Form coherency matrix R(f) = D^{-1} S D^{-H} ------------------
    R = zeros(nfft, N, N);
    for k = 1:nfft
        invDk = diag(1 ./ d(k, :));
        R(k,:,:) = invDk * squeeze(S(k,:,:)) * invDk';
    end

    % --- Wilson on R ---------------------------------------------------
    [~, ~, psi_R, info] = wilson_factor(R, tol, maxiter);

    % --- Recompose psi_S(f) = D(f) * psi_R(f) --------------------------
    psi = zeros(nfft, N, N);
    for k = 1:nfft
        psi(k,:,:) = diag(d(k, :)) * squeeze(psi_R(k,:,:));
    end

    % --- Extract Sigma and H from psi (same logic as wilson_factor) ----
    psi_t0 = real(squeeze(mean(psi, 1)));
    Sigma = psi_t0 * psi_t0.';
    Sigma = (Sigma + Sigma.') / 2;

    H = zeros(nfft, N, N);
    invPsi0 = inv(psi_t0);
    for k = 1:nfft
        H(k,:,:) = squeeze(psi(k,:,:)) * invPsi0;
    end

    info.normalized = true;
end
