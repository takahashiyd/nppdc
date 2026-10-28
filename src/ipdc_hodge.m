function out = ipdc_hodge(I2, m, alpha)
%IPDC_HODGE  Hodge (gradient/ranking) decomposition of the iPDC flow field.
%
%   out = ipdc_hodge(I2)            decomposition only
%   out = ipdc_hodge(I2, m, alpha)  + asymptotic hierarchy test
%
%   Treats the per-frequency squared-iPDC matrix as a directed edge flow on
%   the complete graph K_N and extracts the curl-free (gradient/ranking)
%   component -- the feedforward/feedback hierarchy of the channels. The
%   cyclic (curl) part is reported only as residual energy; it is not
%   decomposed (it depends on the choice of 2-cells/triangles).
%
%   Inputs
%     I2    : nfft x N x N squared iPDC, I2(:,i,j) = |iPDC_{i<-j}|^2 (j -> i)
%     m     : partial-coherence DOF (m = L - N) for the asymptotic test
%     alpha : test level (default 0.05)
%
%   Outputs (per frequency unless noted)
%     out.scores      : nfft x N   ranking potential s (high = sink/receiver)
%     out.div         : nfft x N   divergence (net inflow) of the flow
%     out.grad_energy : nfft x 1   gradient (hierarchy) energy E_g
%     out.tot_energy  : nfft x 1   total edge-flow energy
%     out.grad_frac   : nfft x 1   E_g / tot_energy in [0,1]
%   If m is supplied:
%     out.T           : nfft x 1   chi^2 statistic, ~ chi2_{N-1} under H0
%     out.pval        : nfft x 1   p-value of the global hierarchy test
%     out.reject      : nfft x 1   logical, p < alpha
%     out.df          : N - 1
%     out.score_seH0  : 1 x 1      H0 standard error of each score s_i

    [nf, N, N2] = assertsize(I2);
    if nargin < 3 || isempty(alpha), alpha = 0.05; end

    P = eye(N) - ones(N)/N;          % L0^+ = (1/N) P on K_N (uniform weights)
    iu = triu(true(N), 1);

    out.scores = zeros(nf, N);
    out.div    = zeros(nf, N);
    out.grad_energy = zeros(nf, 1);
    out.tot_energy  = zeros(nf, 1);

    for k = 1:nf
        R   = squeeze(I2(k, :, :));  % R(i,j) = rho_{ij} (j -> i)
        F   = R - R.';               % skew net flow F_ij = rho_ij - rho_ji
        div = sum(F, 2);             % net inflow at each node
        s   = (P * div) / N;         % ranking potential  (sum-zero)
        out.div(k, :)        = div.';
        out.scores(k, :)     = s.';
        out.grad_energy(k)   = sum(div.^2) / N;     % = sum_{i<j} grad_ij^2
        out.tot_energy(k)    = sum(F(iu).^2);
    end
    out.grad_frac = out.grad_energy ./ max(out.tot_energy, realmin);

    if nargin >= 2 && ~isempty(m)
        % Leading-order null that treats the edge flows as independent.
        % Flows sharing a node are correlated, so this test is
        % anti-conservative; IPDC_HODGE_JACK gives calibrated inference.
        varRho = m / ((m+1)^2 * (m+2));   % Var of |iota|^2 under H0 (Beta(1,m))
        sigF2  = 2 * varRho;              % Var of net flow per edge under H0
        out.T       = out.grad_energy / sigF2;       % approx ~ chi2_{N-1}
        out.df      = N - 1;
        out.pval    = 1 - chi2cdf(out.T, N - 1);
        out.reject  = out.pval < alpha;
        out.score_seH0 = sqrt(sigF2 * (N - 1) / N^2);
    end
end

function [nf, N, N2] = assertsize(I2)
    sz = size(I2);
    assert(numel(sz) == 3 && sz(2) == sz(3), ...
        'I2 must be nfft x N x N');
    nf = sz(1); N = sz(2); N2 = sz(3);
end
