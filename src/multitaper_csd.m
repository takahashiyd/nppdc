function [S, f, info] = multitaper_csd(X, fs, NW, K, segLen, overlapFrac, nfft)
%MULTITAPER_CSD  Segment-averaged two-sided multitaper cross-spectral density.
%
%   [S, f] = multitaper_csd(X, fs, NW, K)
%   [S, f] = multitaper_csd(X, fs, NW, K, segLen, overlapFrac, nfft)
%   [S, f, info] = multitaper_csd(...)
%
%   Cuts X into overlapping segments of length segLen, computes a K-taper
%   multitaper CSD on each segment, and averages. Segment averaging is
%   essential here: a single multitaper estimate (K~7) has ~1/sqrt(K) ~ 38%
%   relative SE per frequency, and downstream Wilson factorization is
%   nonlinear enough that this variance creates spurious off-diagonal
%   structure in the inferred transfer function.
%
%   Inputs
%     X           : N x T multivariate time series (real-valued)
%     fs          : sampling rate (Hz)
%     NW          : time-bandwidth product on each segment (e.g., 3 or 4)
%     K           : number of tapers (typically 2*NW - 1)
%     segLen      : segment length in samples       (default: T, single seg)
%     overlapFrac : segment overlap fraction in [0,1) (default: 0.5)
%     nfft        : FFT length, must be >= segLen   (default: segLen)
%
%   Outputs
%     S    : nfft x N x N two-sided complex CSD on f = (0:nfft-1)*fs/nfft
%     f    : nfft x 1 frequency vector (Hz)
%     info : struct with fields .segLen, .nseg, .step, .nfft, .M_indep and
%            .effDOF = (2K+1)*M_indep - 2
%
%   Normalization. With orthonormal DPSS tapers (sum_t v_k(t)^2 = 1), the
%   unweighted taper-average yields a CSD whose IFFT recovers biased
%   covariances at small lags: (1/nfft) * sum_k S(f_k, n, n) ~= var(x_n).
%   Wilson's algorithm expects this convention.

    [N, T] = size(X);
    if nargin < 5 || isempty(segLen),      segLen = T;          end
    if nargin < 6 || isempty(overlapFrac), overlapFrac = 0.5;   end
    if nargin < 7 || isempty(nfft),        nfft = segLen;       end
    if nfft < segLen
        error('multitaper_csd:badNfft', 'nfft (%d) must be >= segLen (%d)', nfft, segLen);
    end
    if segLen > T
        error('multitaper_csd:badSegLen', 'segLen (%d) cannot exceed T (%d)', segLen, T);
    end

    step = max(1, round(segLen * (1 - overlapFrac)));
    nseg = max(1, floor((T - segLen) / step) + 1);

    tapers = dpss(segLen, NW, K);             % segLen x K, orthonormal

    S = zeros(nfft, N, N);
    for s = 1:nseg
        idx = (s-1)*step + (1:segLen);
        if idx(end) > T, break; end
        Xs = X(:, idx);                       % N x segLen

        XF = zeros(nfft, K, N);               % taperedDFT(f, k, n)
        for n = 1:N
            XF(:, :, n) = fft(tapers .* Xs(n, :).', nfft, 1);
        end

        for i = 1:N
            for j = 1:N
                S(:, i, j) = S(:, i, j) + mean(XF(:,:,i) .* conj(XF(:,:,j)), 2);
            end
        end
    end
    S = S / nseg;
    f = (0:nfft-1).' * (fs / nfft);

    info.segLen = segLen;
    info.nseg   = nseg;
    info.step   = step;
    info.nfft   = nfft;
    % Effective complex DOF L of the iPDC estimator, an empirical Monte
    % Carlo calibration (calibrate_dof.m). A single-frequency plug-in
    % coherence would have K*M_indep complex DOF; the roughly twofold larger
    % value reflects the pooling across neighbouring frequencies in Wilson
    % factorization of a smoothed spectrum.
    M_indep        = floor(T / segLen);
    info.effDOF    = (2 * K + 1) * M_indep - 2;
    info.M_indep   = M_indep;
end
