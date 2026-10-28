function out = ipdc_hodge_jack(X, fs, NW, K, segLen, overlapFrac, nfft, alpha, ngroup)
%IPDC_HODGE_JACK  Hodge (gradient/ranking) decomposition of iPDC with
%                 delete-one-block jackknife inference.
%
%   out = ipdc_hodge_jack(X, fs, NW, K, segLen, overlapFrac, nfft, alpha, ngroup)
%
%   Decomposes the per-frequency directed iPDC flow into its curl-free
%   (gradient / ranking) component. The divergence covariance Sigma_div is
%   estimated by a delete-one-block jackknife over non-overlapping blocks of
%   length segLen, which accounts for the correlation between edge flows
%   that share a node. Two tests are returned: a per-node Wald test on the
%   ranking scores (the primary output) and a global Satterthwaite test for
%   the presence of any hierarchy. About 16 blocks (segLen ~ T/16) give
%   calibrated inference (hodge_jack_validate.m).
%
%   ngroup (default min(number of blocks, 40)) groups contiguous blocks into
%   the jackknife units, which bounds the cost on long records. overlapFrac
%   is not used: the jackknife needs independent blocks.
%
%   Outputs (per frequency)
%     out.f            : nfft x 1 frequency vector
%     out.scores       : nfft x N ranking potential (high = sink/receiver)
%     out.grad_energy  : nfft x 1 gradient energy E_g = div' L0^+ div
%     out.grad_frac    : nfft x 1 E_g / total edge-flow energy
%     GLOBAL hierarchy test (H0: gradient component = 0), via Satterthwaite
%       out.T_glob     : nfft x 1 E_g
%       out.p_glob     : nfft x 1 p-value (Satterthwaite-scaled chi^2)
%       out.reject_glob: nfft x 1 logical, p_glob < alpha
%     PER-NODE score test (H0: s_i = 0), Wald with jackknife SE
%       out.score_se   : nfft x N
%       out.score_p    : nfft x N (two-sided)
%       out.score_ci   : nfft x N x 2 (lower, upper) at level alpha
%       out.score_reject : nfft x N logical
%     out.nseg, out.Sigma_div (nfft x N x N)

    if nargin < 8 || isempty(alpha), alpha = 0.05; end
    [N, T] = size(X);
    % non-overlapping blocks: deleting one of several overlapping segments
    % leaves near-duplicate replicates and collapses the variance estimate
    step = segLen;
    nseg = max(2, floor(T / segLen));
    % about 16 blocks are needed for a stable Sigma_div
    if nseg < 12
        warning('ipdc_hodge_jack:fewBlocks', ...
            ['only %d non-overlapping blocks (T/segLen); inference will be ' ...
             'conservative with low power -- use segLen <= T/16.'], nseg);
    end
    tapers = dpss(segLen, NW, K);
    f = (0:nfft-1).' * (fs / nfft);

    % group the blocks into njk contiguous leave-one-out units
    if nargin < 9 || isempty(ngroup), ngroup = min(nseg, 40); end
    njk = min(max(2, ngroup), nseg);

    % ---- accumulate per-group and full cross-spectra (streamed) -------
    Gsum = zeros(nfft, N, N, njk);  Gcnt = zeros(njk, 1);
    Fsum = zeros(nfft, N, N);       nblk = 0;
    for s = 1:nseg
        idx = (s-1)*step + (1:segLen);
        if idx(end) > T, break; end
        Xs = X(:, idx);
        XF = zeros(nfft, K, N);
        for n = 1:N, XF(:, :, n) = fft(tapers .* Xs(n, :).', nfft, 1); end
        Ss = zeros(nfft, N, N);
        for i = 1:N
            for j = 1:N
                Ss(:, i, j) = mean(XF(:, :, i) .* conj(XF(:, :, j)), 2);
            end
        end
        g = floor((s-1) * njk / nseg) + 1;
        Gsum(:, :, :, g) = Gsum(:, :, :, g) + Ss;  Gcnt(g) = Gcnt(g) + 1;
        Fsum = Fsum + Ss;                          nblk = nblk + 1;
    end
    Sfull = Fsum / nblk;

    % ---- full-sample decomposition -----------------------------------
    Hf = ipdc_hodge(local_sq(Sfull));

    % ---- delete-one-group jackknife of the divergence ----------------
    divJ = zeros(njk, nfft, N);
    for g = 1:njk
        S_g = (Fsum - Gsum(:, :, :, g)) / (nblk - Gcnt(g));   % leave group g out
        h   = ipdc_hodge(local_sq(S_g));
        divJ(g, :, :) = h.div;
    end

    % ---- per-frequency inference -------------------------------------
    L0     = N*eye(N) - ones(N);
    L0plus = pinv(L0);                       % = (1/N)(I - J/N)
    z = norminv(1 - alpha/2);

    out.f = f;
    out.scores = Hf.scores;
    out.grad_energy = Hf.grad_energy;
    out.grad_frac   = Hf.grad_frac;
    out.T_glob = Hf.grad_energy;
    out.p_glob = ones(nfft, 1);
    out.score_se = zeros(nfft, N);
    out.score_p  = ones(nfft, N);
    out.score_ci = zeros(nfft, N, 2);
    out.Sigma_div = zeros(nfft, N, N);

    for k = 1:nfft
        dk   = squeeze(divJ(:, k, :));               % njk x N
        dbar = mean(dk, 1);
        Cs   = ((njk-1)/njk) * ((dk - dbar).' * (dk - dbar));     % jackknife cov
        out.Sigma_div(k, :, :) = Cs;

        % global Satterthwaite test on E_g = div' L0+ div
        Mmat = L0plus * Cs;
        tr1 = trace(Mmat); tr2 = trace(Mmat*Mmat);
        if tr1 > 0 && tr2 > 0
            gsat = tr2/tr1; df = tr1^2/tr2;      % Satterthwaite scale, dof
            out.p_glob(k) = 1 - chi2cdf(Hf.grad_energy(k)/gsat, df);
        end

        % per-node score CIs / Wald test
        Covs = L0plus * Cs * L0plus;
        se = sqrt(max(diag(Covs), 0)).';
        s  = Hf.scores(k, :);
        out.score_se(k, :) = se;
        zz = s ./ max(se, realmin);
        out.score_p(k, :)  = 2*(1 - normcdf(abs(zz)));
        out.score_ci(k, :, 1) = s - z*se;
        out.score_ci(k, :, 2) = s + z*se;
    end
    out.reject_glob   = out.p_glob   < alpha;
    out.score_reject  = out.score_p  < alpha;
    out.nseg  = nseg;
    out.ngrp  = njk;
    out.divJ  = divJ;          % njk x nfft x N jackknife divergence replicates
end

function I2 = local_sq(S)
    [H, Sg] = wilson_factor_normalized(S);
    I2 = ipdc(H, Sg).^2;
end
