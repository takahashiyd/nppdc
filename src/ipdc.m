function I = ipdc(H, Sigma)
%IPDC  Information partial directed coherence from (H, Sigma).
%
%   I = ipdc(H, Sigma)
%
%   Computes the information PDC of Takahashi, Baccalá & Sameshima (2010)
%   from a transfer function H(f) and innovation covariance Sigma. The
%   directed influence is from source j (column) to target i (row).
%
%   Definition (Takahashi-Baccala-Sameshima 2010, eq. 6)
%       Abar(f)         = H(f)^{-1}                  (freq-domain VAR coeffs)
%       iota_{ij}(f)    = |Abar_{ij}(f)| / sqrt(sigma_{ii})
%                         ----------------------------------------
%                         sqrt( abar_j(f)^H * Sigma^{-1} * abar_j(f) )
%
%   The denominator is the FULL quadratic form in Sigma^{-1}, not the
%   diagonal-only sum sum_n |Abar_{nj}|^2 / sigma_{nn}. The two agree iff
%   Sigma is diagonal; for general Sigma the diagonal-only form is gPDC,
%   not iPDC.
%
%   Inputs
%     H     : nfft x N x N transfer function
%     Sigma : N x N innovation covariance (Hermitian positive definite)
%
%   Output
%     I     : nfft x N x N iPDC values in [0, 1]

    [nfft, N, ~] = size(H);
    sig = diag(Sigma);
    if any(sig <= 0)
        error('ipdc:nonPosDiag', 'Sigma must have positive diagonal entries');
    end
    inv_sqrt_sig = 1 ./ sqrt(sig);       % N x 1
    Sigma_inv    = inv(Sigma);
    Sigma_inv    = (Sigma_inv + Sigma_inv') / 2;  % enforce Hermitian

    I = zeros(nfft, N, N);
    for k = 1:nfft
        Hk   = squeeze(H(k,:,:));
        Abar = Hk \ eye(N);              % freq-domain VAR coefficient matrix
        num  = abs(Abar) .* inv_sqrt_sig;            % N x N, row-scaled
        den  = zeros(1, N);
        for j = 1:N
            aj     = Abar(:, j);
            qj     = real(aj' * Sigma_inv * aj);
            den(j) = sqrt(max(qj, eps));
        end
        I(k,:,:) = num ./ den;
    end
end
