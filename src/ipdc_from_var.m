function [I, f, H, S] = ipdc_from_var(A, Sigma, fs, nfft)
%IPDC_FROM_VAR  Parametric iPDC from known VAR coefficients (ground truth).
%
%   [I, f, H, S] = ipdc_from_var(A, Sigma, fs, nfft)
%
%   Computes the iPDC directly from the VAR coefficient tensor:
%       Abar(f) = I - sum_r A_r exp(-i 2 pi f r / fs)
%       H(f)    = Abar(f)^{-1}
%       S(f)    = H(f) Sigma H(f)'
%   then calls IPDC(H, Sigma). Used to validate the nonparametric pipeline.
%
%   Inputs
%     A     : N x N x p VAR coefficient tensor
%     Sigma : N x N innovation covariance
%     fs    : sampling rate (Hz)
%     nfft  : number of frequency points (two-sided grid)
%
%   Outputs
%     I : nfft x N x N parametric iPDC
%     f : nfft x 1 frequency vector (Hz)
%     H : nfft x N x N transfer function
%     S : nfft x N x N spectral density

    [N, ~, p] = size(A);
    f = (0:nfft-1).' * (fs / nfft);
    omega = 2*pi*f / fs;
    I_N = eye(N);

    H = zeros(nfft, N, N);
    S = zeros(nfft, N, N);
    for k = 1:nfft
        Abar = I_N;
        for r = 1:p
            Abar = Abar - A(:,:,r) * exp(-1i * omega(k) * r);
        end
        Hk = Abar \ I_N;
        H(k,:,:) = Hk;
        S(k,:,:) = Hk * Sigma * Hk';
    end

    I = ipdc(H, Sigma);
end
