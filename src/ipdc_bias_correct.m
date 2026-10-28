function [I_corr, I2_corr] = ipdc_bias_correct(I, dof, N)
%IPDC_BIAS_CORRECT  Finite-sample bias correction for the iPDC modulus.
%
%   [I_corr, I2_corr] = ipdc_bias_correct(I, dof, N)
%
%   By Takahashi et al. (2010, Thm 1), iPDC equals an ordinary coherence
%   between constructed scalar processes; in sample |i_hat|^2 behaves as a
%   partial coherence with q = N - 1 partials from a CSD estimate with L complex
%   DOF, i.e. a sample coherence with n = L - q = m + 1 DOF (m = L - N,
%   the Beta(1, m) null parameter). It inherits the classical
%   squared-coherence bias (Carter, Knapp & Nuttall 1973; Jarvis & Mitra
%   2001):
%
%       E[|i_hat|^2] ~= |i|^2 + (1 - |i|^2)^2 / n,      n = m + 1 = L - N + 1.
%
%   At the null (|i|^2 = 0) this is exactly the Beta(1, m) mean 1/(m+1);
%   the bias shrinks quadratically in (1 - |i|^2) as the true value grows.
%
%   Inverting the quadratic u^2/n - u + (1 - |i_hat|^2) = 0 in
%   u = 1 - |i|^2 (taking the root that tends to 1 - |i_hat|^2 as
%   n -> inf) gives the bias-corrected estimator
%
%       |i_corr|^2 = 1 - ( n - sqrt(n^2 - 4 n (1 - |i_hat|^2)) ) / 2,
%
%   clipped to [0, 1]. Observed values below the null floor 1/n map to 0.
%
%   Inputs
%     I    : nfft x N x N iPDC modulus (output of ipdc.m)
%     dof  : effective complex DOF L from multitaper_csd info.effDOF
%     N    : number of channels
%
%   Outputs
%     I_corr  : nfft x N x N corrected iPDC modulus (clipped to [0, 1])
%     I2_corr : nfft x N x N corrected |iPDC|^2 (clipped to [0, 1])

    m = dof - N;
    if m <= 4
        error('ipdc_bias_correct:lowDOF', ...
              'Need m = dof - N > 4 (got dof=%d, N=%d)', dof, N);
    end
    n       = m + 1;
    I2      = I.^2;
    disc    = max(n^2 - 4 * n * (1 - I2), 0);
    I2_corr = 1 - (n - sqrt(disc)) / 2;
    I2_corr = min(max(I2_corr, 0), 1);
    I_corr  = sqrt(I2_corr);
end
