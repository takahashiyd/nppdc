function [thr_mod, ci_lo, ci_hi, info] = ipdc_asymp(H, Sigma, dof, alpha)
%IPDC_ASYMP  Asymptotic significance threshold and CI band for iPDC.
%
%   [thr, ci_lo, ci_hi, info] = ipdc_asymp(H, Sigma, dof, alpha)
%
%   Returns the H_0 significance threshold and the (1 - alpha) confidence
%   interval for the iPDC modulus |iota_ij(f)|.
%
%   By Takahashi et al. (2010, Thm 1), iPDC is the coherence between the
%   innovation w_i and the partialized series eta_j (channel j given the
%   other N-1 channels),
%
%       iota_{ij}(f) = C_{w_i, eta_j}(f).
%
%   Both series are estimated from all N channels, so in sample the estimate
%   behaves as a partial coherence with q = N - 1 partials (Brillinger 1981,
%   Thm 8.5.1 and p. 306):
%
%       H_0:  |iota_hat|^2 ~ Beta(1, m),  m = L - N,
%             equivalently m |iota_hat|^2 / (1 - |iota_hat|^2) ~ F_{2, 2m},
%       H_1:  Var(atanh |iota_hat|) ~= 1/(2 n),  n = L - N + 1
%             (Brillinger 1981, eqs. 8.7.16-8.7.21).
%
%   L is the effective complex DOF returned by MULTITAPER_CSD,
%   L = (2 K + 1) M_indep - 2. The threshold is constant across (i, j, f);
%   the band has constant half-width in the atanh domain,
%
%       tanh( atanh|iota_hat| -/+ z_{1-alpha/2} / sqrt(2 n) ).
%
%   Inputs
%     H     : nfft x N x N transfer function (Wilson)
%     Sigma : N x N innovation covariance (Wilson)
%     dof   : effective complex DOF L (multitaper_csd info.effDOF)
%     alpha : significance level (default 0.05)
%
%   Outputs
%     thr_mod : nfft x N x N threshold on |iota_ij(f)|; reject H_0 if
%               the estimate exceeds it
%     ci_lo   : nfft x N x N lower 1-alpha CI on |iota_ij(f)|, clipped to 0
%     ci_hi   : nfft x N x N upper 1-alpha CI on |iota_ij(f)|, clipped to 1
%     info    : .thr_sq, .half_width (atanh domain), .alpha, .dof

    if nargin < 4 || isempty(alpha), alpha = 0.05; end
    [nfft, N, ~] = size(H);

    % threshold: Beta(1, m) null in F form
    m = dof - N;
    assert(m > 0, 'ipdc_asymp requires dof > N (got dof=%d, N=%d)', dof, N);
    F_crit  = finv(1 - alpha, 2, 2*m);
    thr_sq  = F_crit / (m + F_crit);
    thr_mod = sqrt(thr_sq) * ones(nfft, N, N);

    % CI: Fisher-z band, back-transformed
    z          = norminv(1 - alpha/2);
    n          = dof - N + 1;
    half_width = z / sqrt(2 * n);
    I          = ipdc(H, Sigma);
    Iz         = atanh(min(I, 1 - 1e-12));
    ci_lo      = max(tanh(Iz - half_width), 0);
    ci_hi      = min(tanh(Iz + half_width), 1);

    info.thr_sq     = thr_sq;
    info.half_width = half_width;
    info.alpha      = alpha;
    info.dof        = dof;
end
