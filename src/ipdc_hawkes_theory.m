function [I, f, H, S] = ipdc_hawkes_theory(mu, alpha, beta, fs, nfft)
%IPDC_HAWKES_THEORY  Closed-form iPDC for an exponential-kernel Hawkes process.
%
%   [I, f, H, S] = ipdc_hawkes_theory(mu, alpha, beta, fs, nfft)
%
%   The Bartlett spectrum of the multivariate Hawkes process with
%   exponential kernel phi_{ij}(t) = alpha(i,j) * beta * exp(-beta * t) is
%
%       f_Bartlett(omega) = (1 / 2 pi) * H(omega) * D * H(omega)^H,
%
%   where D = diag(lambda) with lambda = (I - alpha)^{-1} mu the
%   stationary rates, and H(omega) = (I - Phi_hat(omega))^{-1} with
%   Phi_hat_{ij}(omega) = alpha(i,j) * beta / (beta + i omega). This is the
%   density per unit angular frequency; the per-Hz density used in the
%   tutorial is 2 pi times larger, which leaves iPDC unchanged.
%
%   This has the same algebraic form as the VAR spectral density,
%   S = H Sigma H^H, with Sigma = D. The iPDC (Takahashi-Baccala-
%   Sameshima 2010, eq. 6) is therefore
%
%       iota_{ij}(omega) = |Abar_{ij}(omega)| / sqrt(lambda_i)
%                          / sqrt( abar_j(omega)^H * D^{-1} * abar_j(omega) ),
%
%   with Abar = I - Phi_hat. Because D = diag(lambda) is diagonal, the
%   quadratic form in the denominator collapses to
%   sum_n |Abar_{nj}|^2 / lambda_n, but the general definition (and the
%   implementation in src/ipdc.m, called below) uses the full form.
%
%   Inputs
%     mu     : N x 1 baseline rates
%     alpha  : N x N excitation matrix (rho(alpha) < 1 for stability)
%     beta   : scalar decay rate
%     fs     : sampling rate of the Hz grid (Hz). Use the binned-data
%              sampling rate to match the nonparametric estimate.
%     nfft   : number of frequency points (two-sided grid)
%
%   Outputs
%     I : nfft x N x N theoretical iPDC
%     f : nfft x 1 Hz grid (0:nfft-1)*fs/nfft
%     H : nfft x N x N transfer function
%     S : nfft x N x N Bartlett spectrum

    N = numel(mu);
    lambda = (eye(N) - alpha) \ mu(:);

    f = (0:nfft-1).' * (fs / nfft);
    omega = 2*pi * f;
    H = zeros(nfft, N, N);
    S = zeros(nfft, N, N);
    Dm = diag(lambda);
    I_N = eye(N);

    for k = 1:nfft
        Phi = alpha * (beta / (beta + 1i * omega(k)));
        A   = I_N - Phi;
        Hk  = A \ I_N;
        H(k,:,:) = Hk;
        S(k,:,:) = Hk * Dm * Hk' / (2*pi);
    end

    I = ipdc(H, Dm);
end
