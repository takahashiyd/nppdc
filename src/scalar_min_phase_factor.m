function d = scalar_min_phase_factor(Sd)
%SCALAR_MIN_PHASE_FACTOR  Min-phase factor of a scalar autospectrum (cepstrum).
%
%   d = scalar_min_phase_factor(Sd)
%
%   Given a real positive autospectrum Sd(f) sampled on a uniform two-sided
%   frequency grid (length nfft), returns a complex minimum-phase factor
%   d(f) such that |d(f)|^2 == Sd(f) at every grid point.
%
%   Algorithm (cepstral / "homomorphic" factorization)
%     c(tau)   = ifft(log Sd(f))       % real cepstrum, even-symmetric in tau
%     c+(0)    = 0.5 * c(0)
%     c+(tau)  = c(tau)   for tau > 0
%     c+(tau)  = 0        for tau < 0
%     log d(f) = fft(c+(tau))          % causal -> minimum phase
%     d(f)     = exp(log d(f))
%
%   Used by wilson_factor_normalized to precondition the cross-spectrum
%   into the coherency matrix R(f) = D(f)^{-1} S(f) D(f)^{-H}, which has
%   unit diagonal and bounded off-diagonals.

    Sd   = max(real(Sd(:)), eps);          % numerical floor
    nfft = numel(Sd);
    ceps = real(ifft(log(Sd)));
    ceps(1) = 0.5 * ceps(1);
    ceps(floor(nfft/2)+1:nfft) = 0;        % zero out negative lags + Nyquist
    d = exp(fft(ceps));
end
